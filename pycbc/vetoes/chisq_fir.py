"""Signal-consistency estimators and filter spectra for hierarchical / ratio FIR
matched filtering.
"""

import types

import numpy as np

from pycbc.types import FrequencySeries
from pycbc.vetoes.chisq import (
    SingleDetPowerChisq,
    power_chisq_bin_weights,
    power_chisq_bins_from_sigmasq_series,
)

try:
    import matchedfilter as _mf
except ImportError:
    _mf = None


def _require_matchedfilter():
    if _mf is None:
        raise ImportError(
            "matchedfilter is required for ratio/FIR filtering but is not installed"
        )


def _psd_cache_key(stilde):
    psd = getattr(stilde, "psd", None)
    return None if psd is None else getattr(psd, "_cache_key", id(psd))


def _stilde_cache_key(stilde):
    return getattr(stilde, "_cache_key", id(stilde))


class ChisqFilters:
    """A fine group's chisq filters: each template's spectrum on the chisq block
    (length block_length, conjugated, centre tap at lag 0), made from its taps
    when first read -- most templates never have a trigger loud enough for a chisq.

    The chisq is computed on a block of its own, centred on each trigger, so its
    resolution is pycbc's choice and independent of how matchedfilter blocks the
    series for filtering. The block must hold every filter with room for the lags
    evaluated around the trigger (reach either side, for the auto-chisq); otherwise
    this raises rather than lengthen the block, since its length is part of the
    statistic. valid: (T, 2) [start, stop) of each filter's valid lags in its block.
    """

    def __init__(
        self,
        taps,
        counts,
        block_length,
        reach=0,
        tap_sample_rate=2048,
        data_sample_rate=2048,
    ):
        _require_matchedfilter()
        n = self.block_length = int(block_length)
        self._taps = taps
        self._counts = np.asarray(counts, dtype=np.int64)
        self._rates = dict(
            tap_sample_rate=tap_sample_rate, data_sample_rate=data_sample_rate
        )
        self.valid = _mf.template_valid(self._counts, n, **self._rates)
        centre = n // 2
        short = np.flatnonzero(
            (centre - reach < self.valid[:, 0]) | (centre + reach >= self.valid[:, 1])
        )
        if len(short):
            lo, hi = self.valid[short[0]]
            raise ValueError(
                "--chisq-block-length %d is too short for this template group: a "
                "ratio filter has valid lags [%d, %d), not the block centre "
                "+-%d" % (n, lo, hi, reach)
            )
        self._rows = {}

    def __len__(self):
        return len(self._counts)

    def prefetch(self, indices):
        """Batch-compute and cache spectra for any uncached row indices in `indices`."""
        n = len(self._counts)
        missing = [
            i
            for i in dict.fromkeys(int(x) % n for x in np.asarray(indices).ravel())
            if 0 <= i < n and i not in self._rows
        ]
        if not missing:
            return
        idx = np.asarray(missing, dtype=np.int64)
        spectra, _ = _mf.template_spectra(
            self._taps[idx], self._counts[idx], self.block_length, **self._rates
        )
        for k, row_idx in enumerate(missing):
            self._rows[row_idx] = spectra[k]

    def __getitem__(self, i):
        if not isinstance(i, (int, np.integer)):
            raise TypeError(
                f"ChisqFilters index must be an integer, got {type(i).__name__}"
            )
        i = int(i)
        n = len(self._counts)
        if i < -n or i >= n:
            raise IndexError(f"ChisqFilters index {i} out of range for {n} templates")
        if i < 0:
            i += n
        row = self._rows.get(i)
        if row is None:
            spectra, _ = _mf.template_spectra(
                self._taps[i : i + 1],
                self._counts[i : i + 1],
                self.block_length,
                **self._rates,
            )
            row = self._rows[i] = spectra[0]
        return row


def chisq_middle_spectra(taps, counts, ids, block_length):
    """Each middle filter's spectrum for the power chisq, keyed by middle id.

    On the smallest power-of-two multiple of the chisq block that holds the filter,
    so decimating by an integer gives its exact DTFT samples on the chisq grid; in
    the convention of TimeDomainFilterBank.filters_f (conjugate, centre tap at lag 0).
    """
    _require_matchedfilter()
    counts = np.asarray(counts, dtype=np.int64)
    taps = np.asarray(taps, dtype=np.float32)
    ns = np.array(
        [max(int(block_length), 1 << int(np.ceil(np.log2(int(c) + 1)))) for c in counts]
    )
    out = {}
    for n in np.unique(ns):
        rows = np.flatnonzero(ns == n)
        m = int(counts[rows].max())
        spec, _ = _mf.template_spectra(taps[rows, :m], counts[rows], int(n))
        for r, row in enumerate(rows):
            out[int(ids[row])] = spec[r]
    return out


class RatioPowerChisq(SingleDetPowerChisq):
    """Fast power chisq and block template/correlation reconstruction for
    hierarchical ratio/FIR filtering.
    """

    def __init__(
        self,
        num_bins,
        sample_rate,
        flow=20.0,
        f_final=800.0,
        top_template=None,
        upper_f_by_coarse=None,
        bank=None,
        snr_threshold=None,
    ):
        super().__init__(num_bins=num_bins, snr_threshold=snr_threshold)
        if self.do and self.num_bins is not None:
            self.num_bins = int(self.num_bins)
        self.sample_rate = float(sample_rate)
        self.flow = float(flow)
        self.f_final = float(f_final)
        self.top_template = top_template
        self.upper_f_by_coarse = (
            upper_f_by_coarse if upper_f_by_coarse is not None else {}
        )
        self.bank = bank
        self._batch = None

        self.template_cache = {}
        self.psd_cache = {}
        self.coarse_w_cache = {}
        self.ref_f_cache = {}

    def clear(self):
        self.template_cache.clear()
        self.psd_cache.clear()
        self.coarse_w_cache.clear()
        self.ref_f_cache.clear()
        self._batch = None

    def begin_group(self):
        """Release the block spectra of the previous (segment, middle group) pass."""
        self.ref_f_cache.clear()

    def set_batch(
        self, stilde, coarse_idx, ids_local, gids, t_idxs, ref_snr, filters_f
    ):
        """Bind the current trigger batch for template(j) and correlation(j, reach)."""
        self.ref_f_cache.clear()
        self._batch = (
            stilde,
            coarse_idx,
            ids_local,
            gids,
            t_idxs,
            ref_snr,
            filters_f,
        )

    def coarse_weight(self, stilde, coarse_idx, block_len):
        """(psd key, whitened top x upper spectrum of a middle group) on a block's
        frequency grid: (top_dec * upper_f) / sqrt(psd), None without a top template.
        """
        df = self.sample_rate / block_len
        n_half = block_len // 2 + 1
        psd_key = _psd_cache_key(stilde)

        if psd_key not in self.psd_cache:
            dec = max(1, int(round(df / stilde.delta_f)))
            self.psd_cache[psd_key] = stilde.psd[::dec][:n_half].numpy()
        psd_pad = self.psd_cache[psd_key]

        coarse_key = (psd_key, coarse_idx)
        if coarse_key not in self.coarse_w_cache:
            if self.top_template is not None:
                dec_top = max(1, int(round(df / self.top_template.delta_f)))
                top_pts = self.top_template[::dec_top]
                m_top = min(len(top_pts), n_half)
                top_dec = np.zeros(n_half, dtype=np.complex64)
                top_dec[:m_top] = top_pts[:m_top]
                upper_f = self.upper_f_by_coarse.get(coarse_idx)
                if upper_f is not None:
                    step_u = max(1, len(upper_f) // block_len)
                    coarse_dec = top_dec * np.conj(upper_f[::step_u][:n_half])
                else:
                    coarse_dec = top_dec
                kmin = int(self.flow / df)
                kmax = min(int(self.f_final / df), n_half - 1)
                cw = (coarse_dec / np.sqrt(np.maximum(psd_pad, 1e-30))).astype(
                    np.complex64
                )
                cw[:kmin] = 0.0
                cw[kmax:] = 0.0
                self.coarse_w_cache[coarse_key] = cw
            else:
                self.coarse_w_cache[coarse_key] = None
        return psd_key, self.coarse_w_cache[coarse_key]

    def ref_spectrum(self, stilde, coarse_idx, ref_snr, tstart, block_len):
        """Band-limited FFT of the block of the middle series starting at tstart,
        or None if the series is too short.
        """
        cache_key = (_stilde_cache_key(stilde), coarse_idx, tstart, block_len)
        if cache_key not in self.ref_f_cache:
            ref = ref_snr[tstart : tstart + block_len]
            if ref is None or len(ref) < block_len:
                return None
            df = self.sample_rate / block_len
            kmin = int(self.flow / df)
            kmax = min(int(self.f_final / df), block_len // 2)
            spec = np.fft.fft(ref).astype(np.complex64)
            spec[:kmin] = 0.0
            spec[kmax:] = 0.0
            self.ref_f_cache[cache_key] = spec
        return self.ref_f_cache[cache_key]

    def template(self, j):
        """Return the whitened fine template for batch trigger j on the block grid
        with its power-chisq bins and SG-chisq metadata attached, cached per
        (psd, coarse_idx, local_id).
        """
        stilde, coarse_idx, ids_local, gids, _t_idxs, _ref_snr, filters_f = self._batch
        local_id = int(ids_local[j])
        gid = int(gids[j])
        bank = self.bank
        block_len = filters_f.block_length
        df = self.sample_rate / block_len
        n_half = block_len // 2 + 1
        psd_key, coarse_w = self.coarse_weight(stilde, coarse_idx, block_len)
        key = (psd_key, coarse_idx, local_id)
        if key not in self.template_cache:
            fine = np.zeros(block_len, dtype=np.complex64)
            bins = None
            chisq_plan = None
            if coarse_w is not None:
                fine_w = coarse_w * np.conj(filters_f[local_id][:n_half])
                fine[:n_half] = fine_w
                kmin = int(self.flow / df)
                kmax = min(int(self.f_final / df), n_half - 1)
                s_vec = np.cumsum(np.abs(fine_w, dtype=np.float64) ** 2)
                if s_vec[kmax - 1] > 0 and self.do:
                    raw_bins = power_chisq_bins_from_sigmasq_series(
                        s_vec, self.num_bins, kmin, kmax
                    )
                    bins, weights = power_chisq_bin_weights(s_vec, raw_bins)
                    bins = bins.astype(np.int32)
                    p = len(bins) - 1
                    if p >= 2:
                        inv_wb = weights.astype(np.float32)
                        k_freqs = np.arange(bins[0], bins[-1], dtype=np.float32)
                        center_phases = np.exp(
                            (2j * np.pi * (block_len // 2) / block_len) * k_freqs
                        ).astype(np.complex64)
                        bin_starts = bins[:-1] - bins[0]
                        chisq_plan = (
                            bins,
                            p,
                            inv_wb,
                            k_freqs,
                            center_phases,
                            bin_starts,
                        )
            tmplt = FrequencySeries(fine, delta_f=df)
            tmplt.bins = bins
            tmplt._chisq_plan = chisq_plan
            if bank is not None and getattr(bank, "table", None) is not None:
                tmplt.f_lower = float(bank.table["f_lower"][gid])
                tmplt.params = types.SimpleNamespace(
                    template_hash=bank.table["template_hash"][gid]
                )
            else:
                tmplt.f_lower = float(self.flow)
                tmplt.params = types.SimpleNamespace(template_hash=gid)
            self.template_cache[key] = tmplt
        return self.template_cache[key]

    def correlation(self, j, reach=0):
        """Return (corr, index_in_block) for batch trigger j on the chisq block
        centred on t_idx, or (None, None) if reach exceeds the valid FIR lags.
        """
        stilde, coarse_idx, ids_local, _gids, t_idxs, ref_snr, filters_f = self._batch
        local_id = int(ids_local[j])
        t_idx = int(t_idxs[j])
        block_len = filters_f.block_length
        tstart = min(max(t_idx - block_len // 2, 0), len(ref_snr) - block_len)
        if reach > 0 and hasattr(filters_f, "valid"):
            valid_start, valid_stop = filters_f.valid[local_id]
            if (
                t_idx - reach - tstart < valid_start
                or t_idx + reach - tstart >= valid_stop
            ):
                return None, None
        ref_f = self.ref_spectrum(stilde, coarse_idx, ref_snr, tstart, block_len)
        if ref_f is None:
            return None, None
        return (filters_f[local_id] * ref_f) / block_len, t_idx - tstart

    def values(self, corr, snrv, snr_norm, psd, indices, template):
        """Power chisq with unequal-weight bins, for SignalConsistencyChain."""
        plan = getattr(template, "_chisq_plan", None)
        if corr is None or plan is None:
            return None, None
        bins, p, inv_wb, k_freqs, center_phases, bin_starts = plan
        dt = int(indices[0])
        block_len = len(corr)
        if dt == block_len // 2:
            phases = center_phases
        else:
            phases = np.exp((2j * np.pi * dt / block_len) * k_freqs).astype(
                np.complex64
            )
        phased_corr = corr[bins[0] : bins[-1]] * phases
        z_b = np.add.reduceat(phased_corr, bin_starts)
        rho_sum = np.sum(z_b)
        chi_val = max(
            0.0, float(np.sum((np.abs(z_b) ** 2) * inv_wb) - abs(rho_sum) ** 2)
        )
        return [float(chi_val)], [2 * p - 2]
