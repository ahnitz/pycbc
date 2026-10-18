"""Signal-consistency (power, auto, and sine-Gaussian chisq) estimators for
hierarchical / ratio FIR matched filtering.
"""

import time
import types
import numpy as np

import pycbc.pnutils
from pycbc.types import FrequencySeries, float32
from pycbc.vetoes.autochisq import SingleDetAutoChisq, autochisq_from_precomputed
from pycbc.vetoes.chain import SignalConsistencyChain
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
    if psd is None:
        return None
    return getattr(psd, "_cache_key", id(psd))


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
            (centre - reach < self.valid[:, 0])
            | (centre + reach >= self.valid[:, 1])
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
        missing = []
        seen = set()
        for raw_i in np.asarray(indices).ravel():
            i = int(raw_i)
            if i < 0:
                i += n
            if 0 <= i < n and i not in self._rows and i not in seen:
                seen.add(i)
                missing.append(i)
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
            raise IndexError(
                f"ChisqFilters index {i} out of range for {n} templates"
            )
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
        [
            max(int(block_length), 1 << int(np.ceil(np.log2(int(c) + 1))))
            for c in counts
        ]
    )
    out = {}
    for n in np.unique(ns):
        rows = np.flatnonzero(ns == n)
        m = int(counts[rows].max())
        spec, _ = _mf.template_spectra(taps[rows, :m], counts[rows], int(n))
        for r, row in enumerate(rows):
            out[int(ids[row])] = spec[r]
    return out


class RatioBlockInputs:
    """Inputs for SignalConsistencyChain on one middle group's triggers.
    Each trigger's chisq block is centred on it (clipped to the series).
    """

    def __init__(
        self,
        power_chisq,
        ids_local,
        t_idxs,
        snrs,
        fine_indices,
        stilde,
        ref_snr,
        filters_f,
        coarse_idx,
        bank=None,
        shifts=None,
    ):
        self.power_chisq = power_chisq
        self.ids_local = np.asarray(ids_local)
        self.t_idxs = np.asarray(t_idxs)
        self.snrs = np.asarray(snrs)
        self.fine_indices = np.asarray(fine_indices)
        self.stilde = stilde
        self.ref_snr = ref_snr
        self.filters_f = filters_f
        self.coarse_idx = coarse_idx
        self.bank = bank
        local_shifts = shifts[self.ids_local] if shifts is not None else 0
        self.phys_t = self.t_idxs - local_shifts
        self.template_id = self.fine_indices[self.ids_local]
        self.block_len = filters_f.block_length
        self.power_chisq.coarse_weight(stilde, coarse_idx, self.block_len)
        self.templates = {}

    def template(self, j):
        """The whitened fine template on the block's grid, with its power
        chisq bins (None if unusable) and template metadata for SG chisq.
        """
        local_id = int(self.ids_local[j])
        if local_id not in self.templates:
            n_half = self.block_len // 2 + 1
            coarse_w = self.power_chisq.coarse_weight(
                self.stilde, self.coarse_idx, self.block_len
            )[1]
            fine = np.zeros(self.block_len, dtype=np.complex64)
            if coarse_w is not None:
                fine[:n_half] = coarse_w * np.conj(
                    self.filters_f[local_id][:n_half]
                )
            tmplt = FrequencySeries(
                fine, delta_f=self.power_chisq.sample_rate / self.block_len
            )
            tmplt.bins = self.power_chisq.bins(
                self.stilde,
                self.coarse_idx,
                local_id,
                self.filters_f,
                self.block_len,
            )
            gid = int(self.template_id[j])
            if self.bank is not None and getattr(self.bank, "table", None) is not None:
                row = self.bank.table[gid]
                tmplt.f_lower = float(self.bank.table["f_lower"][gid])
                tmplt.params = types.SimpleNamespace(
                    template_hash=self.bank.table["template_hash"][gid]
                )
                if tmplt.bins is not None:
                    tmplt.fpeak = float(tmplt.bins[-1])
                else:
                    tmplt.fpeak = float(
                        pycbc.pnutils.named_frequency_cutoffs.get(
                            "IMRPhenomDPeak", lambda r: 100.0
                        )(row)
                    )
            else:
                tmplt.f_lower = float(self.power_chisq.flow)
                tmplt.params = types.SimpleNamespace(template_hash=gid)
                if tmplt.bins is not None:
                    tmplt.fpeak = float(tmplt.bins[-1])
            self.templates[local_id] = tmplt
        return self.templates[local_id]

    def correlation(self, j, reach=0):
        local_id = int(self.ids_local[j])
        t_idx = int(self.t_idxs[j])
        tstart = min(
            max(t_idx - self.block_len // 2, 0),
            len(self.ref_snr) - self.block_len,
        )
        if reach > 0 and hasattr(self.filters_f, "valid"):
            valid_start, valid_stop = self.filters_f.valid[local_id]
            if (
                t_idx - reach - tstart < valid_start
                or t_idx + reach - tstart >= valid_stop
            ):
                return None, None
        ref_f = self.power_chisq.ref_spectrum(
            self.stilde, self.coarse_idx, self.ref_snr, tstart, self.block_len
        )
        if ref_f is None:
            return None, None
        corr = (self.filters_f[local_id] * ref_f) / self.block_len
        return corr, t_idx - tstart


class RatioPowerChisq(SingleDetPowerChisq):
    """Calculates fast power chisq for ratio/FIR filtering using reconstructed
    frequency-domain template spectra and exact template-power bin weighting.
    """

    extra_columns = {"fpeak": float32}

    def __init__(
        self,
        num_bins,
        snr_threshold,
        sample_rate,
        flow=20.0,
        f_final=800.0,
        chisq_window=1.0,
        chisq_iterations=None,
        top_template=None,
        upper_f_by_coarse=None,
    ):
        super().__init__(num_bins=num_bins, snr_threshold=snr_threshold)
        if self.do and self.num_bins is not None:
            self.num_bins = int(self.num_bins)
        if self.snr_threshold is not None:
            self.snr_threshold = float(self.snr_threshold)
        self.sample_rate = float(sample_rate)
        self.flow = float(flow)
        self.f_final = float(f_final)
        self.chisq_window_samp = int(chisq_window * self.sample_rate)
        self.max_iterations = chisq_iterations
        self.top_template = top_template
        self.upper_f_by_coarse = (
            upper_f_by_coarse if upper_f_by_coarse is not None else {}
        )

        self.bins_cache = {}
        self.psd_cache = {}
        self.coarse_w_cache = {}
        self.ref_f_cache = {}
        self.chains = {}
        self.last_extra = {"fpeak": 0.0}

    def clear(self):
        self.bins_cache.clear()
        self.psd_cache.clear()
        self.coarse_w_cache.clear()
        self.ref_f_cache.clear()
        self.chains.clear()

    def begin_group(self):
        """Release the block spectra of the previous (segment, middle group) pass."""
        self.ref_f_cache.clear()

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
                    upper_dec = np.conj(upper_f[::step_u][:n_half])
                    coarse_dec = top_dec * upper_dec
                else:
                    coarse_dec = top_dec
                kmax_w = min(int(self.f_final / df), n_half - 1)
                coarse_dec[kmax_w:] = 0.0
                psd_safe = np.maximum(psd_pad, 1e-30)
                cw = (coarse_dec / np.sqrt(psd_safe)).astype(np.complex64)
                cw[kmax_w:] = 0.0
                self.coarse_w_cache[coarse_key] = cw
            else:
                self.coarse_w_cache[coarse_key] = None
        return psd_key, self.coarse_w_cache[coarse_key]

    def ref_spectrum(self, stilde, coarse_idx, ref_snr, tstart, block_len):
        """FFT of the block of the middle series starting at tstart, or None
        if the series is too short.
        """
        cache_key = (_stilde_cache_key(stilde), coarse_idx, tstart, block_len)
        if cache_key not in self.ref_f_cache:
            ref = ref_snr[tstart : tstart + block_len]
            if ref is None or len(ref) < block_len:
                return None
            self.ref_f_cache[cache_key] = np.fft.fft(ref).astype(np.complex64)
        return self.ref_f_cache[cache_key]

    def bins(self, stilde, coarse_idx, local_id, filters_f, block_len):
        """(bins, p, inv_wb, k_freqs, center_phases, bin_starts, fpeak) of a
        template, or None if it has no usable bins.
        """
        df = self.sample_rate / block_len
        n_half = block_len // 2 + 1
        psd_key, coarse_w = self.coarse_weight(stilde, coarse_idx, block_len)
        bins_key = (psd_key, coarse_idx, local_id)
        if bins_key not in self.bins_cache:
            if coarse_w is None:
                return None
            ratio = filters_f[local_id]
            ratio_dec = np.conj(ratio[:n_half])
            fine_w = coarse_w * ratio_dec

            kmin = int(self.flow / df)
            kmax = min(int(self.f_final / df), n_half - 1)
            power = np.abs(fine_w, dtype=np.float64) ** 2
            power[:kmin] = 0.0
            power[kmax:] = 0.0

            s_vec = np.zeros(n_half, dtype=np.float64)
            s_vec[kmin:kmax] = np.cumsum(power[kmin:kmax], dtype=np.float64)
            s_vec[kmax:] = s_vec[kmax - 1]
            if s_vec[kmax - 1] <= 0:
                self.bins_cache[bins_key] = None
                return None

            raw_bins = power_chisq_bins_from_sigmasq_series(
                s_vec, self.num_bins, kmin, kmax
            )
            bins, weights = power_chisq_bin_weights(s_vec, raw_bins)
            bins = bins.astype(np.int32)
            p = len(bins) - 1
            if p < 2:
                self.bins_cache[bins_key] = None
                return None
            inv_wb = weights.astype(np.float32)
            k_freqs = np.arange(bins[0], bins[-1], dtype=np.float32)
            center_phases = np.exp(
                (2j * np.pi * (block_len // 2) / block_len) * k_freqs
            ).astype(np.complex64)
            bin_starts = bins[:-1] - bins[0]
            fstep = float(bins[-2] - bins[-3]) * df
            fpeak = float(bins[-2]) * df + fstep
            self.bins_cache[bins_key] = (
                bins,
                p,
                inv_wb,
                k_freqs,
                center_phases,
                bin_starts,
                fpeak,
            )
        return self.bins_cache.get(bins_key)

    def values(self, corr, snrv, snr_norm, psd, indices, template):
        """Power chisq with unequal-weight bins, for SignalConsistencyChain."""
        if corr is None or getattr(template, "bins", None) is None:
            return None, None
        bins, p, inv_wb, k_freqs, center_phases, bin_starts, fpeak = template.bins
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
        self.last_extra = {"fpeak": float(fpeak)}
        return [float(chi_val)], [2 * p - 2]

    def chain(self, autochisq=None, max_rchi=None, sg_chisq=None, steps=None):
        """Build or retrieve the SignalConsistencyChain for this configuration."""
        auto_on = autochisq is not None and getattr(autochisq, "do", False)
        sg_on = sg_chisq is not None and getattr(sg_chisq, "do", False)
        steps_key = tuple(steps) if steps is not None else None
        key = (auto_on, max_rchi, sg_on, steps_key)
        if key not in self.chains:
            num_b = int(self.num_bins) if self.num_bins is not None else 64
            defaults = {
                "chisq": 2 * num_b - 2,
                "chisq_dof": 2 * num_b - 2,
                "fpeak": 0.0,
            }
            if auto_on:
                defaults["cont_chisq"] = autochisq.dof
            if sg_on:
                defaults["sg_chisq"] = 1.0
            if steps is None:
                chain_steps = []
                if self.do:
                    if self.snr_threshold:
                        chain_steps.append(f"snr:{self.snr_threshold!r}:lower:skip")
                    chain_steps.append("power")
                if auto_on:
                    if max_rchi is not None and max_rchi > 0:
                        chain_steps.append(f"traditional:{max_rchi!r}:upper_inc:skip")
                    chain_steps.append("auto")
                if sg_on:
                    chain_steps.append(
                        f"newsnr:{sg_chisq.snr_threshold!r}:lower_inc:skip"
                    )
                    chain_steps.append("sg")
            else:
                chain_steps = list(steps)
            self.chains[key] = SignalConsistencyChain(
                self if self.do else None,
                autochisq if auto_on else None,
                sg_chisq if sg_on else None,
                steps=chain_steps,
                window=self.chisq_window_samp,
                max_iterations=self.max_iterations,
                defaults=defaults,
            )
        return self.chains[key]

    def compute(
        self,
        ids_local,
        t_idxs,
        snrs,
        fine_indices,
        stilde,
        ref_snr,
        filters_f,
        coarse_idx,
        args,
        autochisq=None,
        sg_chisq=None,
        bank=None,
        shifts=None,
    ):
        """Evaluate power chisq (and optional auto-chisq and SG chisq) through
        SignalConsistencyChain.

        Returns (chisq_vals, chisq_dof_vals, auto_vals, fpeak_vals) when called
        without sg_chisq, or (keep, vals, stage_times) when `sg_chisq` is passed
        or `return_chain_dict` is used.
        """
        snrs_arr = np.asarray(snrs)
        snr_amp = np.abs(snrs_arr)
        num_trigs = len(snr_amp)
        num_b = int(self.num_bins) if self.num_bins is not None else 64
        auto_on = autochisq is not None and getattr(autochisq, "do", False)

        cthresh = self.snr_threshold or 0.0
        above = snr_amp >= cthresh
        self.last_kept = np.ones(num_trigs, dtype=bool)
        if not np.any(above):
            chisq_vals = np.full(num_trigs, 2 * num_b - 2, dtype=float32)
            chisq_dof_vals = np.full(num_trigs, 2 * num_b - 2, dtype=np.int32)
            auto_vals = (
                np.full(num_trigs, autochisq.dof, dtype=float32)
                if auto_on
                else None
            )
            fpeak_vals = np.zeros(num_trigs, dtype=float32)
            if sg_chisq is not None and getattr(sg_chisq, "do", False):
                self.last_sg_vals = np.ones(num_trigs, dtype=float32)
            return chisq_vals, chisq_dof_vals, auto_vals, fpeak_vals

        if hasattr(filters_f, "prefetch"):
            filters_f.prefetch(np.asarray(ids_local)[above])

        max_rchi = getattr(
            args,
            "autochi_max_rchisq",
            getattr(args, "autochisq_max_rchisq", 5.0),
        )
        custom_steps = getattr(args, "signal_consistency_chain", None)
        chain = self.chain(
            autochisq=autochisq if auto_on else None,
            max_rchi=max_rchi,
            sg_chisq=sg_chisq,
            steps=custom_steps,
        )
        inputs = RatioBlockInputs(
            self,
            ids_local,
            t_idxs,
            snrs_arr,
            fine_indices,
            stilde,
            ref_snr,
            filters_f,
            coarse_idx,
            bank=bank,
            shifts=shifts,
        )
        try:
            kept, vals = chain.values(
                inputs.snrs,
                inputs.phys_t,
                inputs.template_id,
                inputs.template,
                inputs.correlation,
                stilde=stilde,
            )
        finally:
            self.ref_f_cache.clear()

        self.last_kept = kept
        chisq_vals = vals["chisq"].astype(float32, copy=False)
        chisq_dof_vals = vals["chisq_dof"].astype(np.int32, copy=False)
        auto_vals = (
            vals["cont_chisq"].astype(float32, copy=False) if auto_on else None
        )
        fpeak_vals = vals.get(
            "fpeak", np.zeros(num_trigs, dtype=float32)
        ).astype(float32, copy=False)
        if sg_chisq is not None and getattr(sg_chisq, "do", False):
            self.last_sg_vals = vals.get(
                "sg_chisq", np.ones(num_trigs, dtype=float32)
            ).astype(float32, copy=False)
        return chisq_vals, chisq_dof_vals, auto_vals, fpeak_vals


class RatioAutoChisq(SingleDetAutoChisq):
    """pycbc's auto-chisq on the SNR series the ratio filter reconstructs.

    Options, statistic and degrees of freedom are SingleDetAutoChisq's; only its
    two inputs are built differently, from the block spectra the power chisq
    uses. Over an N-sample block of the middle series with FFT D, the template's
    SNR is IFFT(filters_f * D), and its autocorrelation is the IFFT of
    its whitened power |coarse_w * conj(filters_f)|^2. Each trigger gets a block
    centred on it, so its lags stay inside the block's valid samples.
    """

    def __init__(self, power_chisq, args):
        if getattr(args, "autochi_stride", 0) > 0 and getattr(
            args, "autochi_reverse_template", False
        ):
            raise ValueError(
                "--autochi-reverse-template is not supported by pycbc_inspiral_fir"
            )
        super().__init__(
            args.autochi_stride,
            args.autochi_number_points,
            onesided=args.autochi_onesided,
            twophase=args.autochi_two_phase,
            take_maximum_value=args.autochi_max_valued,
            maximal_value_dof=args.autochi_max_valued_dof,
        )
        if self.do and not power_chisq.do:
            raise ValueError(
                "auto-chisq is computed with the power chisq: it needs --chisq-bins"
            )
        self.power_chisq = power_chisq
        self.autocors = {}
        self.elapsed = 0.0

    @staticmethod
    def insert_option_group(parser):
        SingleDetAutoChisq.insert_option_group(parser)
        parser.add_argument(
            "--autochi-max-rchisq",
            "--autochisq-max-rchisq",
            type=float,
            default=5.0,
            metavar="RCHISQ",
            help="Skip the auto-chisq when a trigger's reduced power chisq "
            "(chisq / dof) is above RCHISQ; skipped triggers keep the "
            "default auto-chisq (dof, reduced value 1).",
        )

    def clear(self):
        self.autocors.clear()

    def autocorrelation(self, stilde, coarse_idx, local_id, filters_f, band):
        N = len(filters_f[local_id])
        psd_key, coarse_w = self.power_chisq.coarse_weight(stilde, coarse_idx, N)
        key = (psd_key, coarse_idx, local_id)
        if key not in self.autocors:
            power = np.zeros(N)
            power[band] = (
                np.abs(coarse_w[band] * np.conj(filters_f[local_id][band])) ** 2
            )
            autocor = np.fft.ifft(power)
            reach = self.stride * self.num_points
            lags = np.r_[0 : reach + 1, N - reach : N]
            self.autocors[key] = autocor[lags] / autocor[0]
        return self.autocors[key]

    def value(self, stilde, coarse_idx, ref_snr, filters_f, local_id, t):
        """Auto-chisq of one trigger, or None if its lags would leave the series."""
        t0 = time.perf_counter()
        reach = self.stride * self.num_points
        N = len(filters_f[local_id])
        start = min(max(int(t) - N // 2, 0), len(ref_snr) - N)
        valid_start, valid_stop = filters_f.valid[local_id]
        if t - reach - start < valid_start or t + reach - start >= valid_stop:
            return None
        ref_f = self.power_chisq.ref_spectrum(
            stilde, coarse_idx, ref_snr, start, N
        )
        if ref_f is None:
            return None
        df = self.power_chisq.sample_rate / N
        band = slice(
            int(self.power_chisq.flow / df),
            min(int(self.power_chisq.f_final / df), N // 2),
        )
        corr = np.zeros(N, dtype=np.complex128)
        corr[band] = filters_f[local_id][band] * ref_f[band]
        snr = np.fft.ifft(corr)[t - start - reach : t - start + reach + 1]
        _, achisq = autochisq_from_precomputed(
            snr,
            snr,
            self.autocorrelation(
                stilde, coarse_idx, local_id, filters_f, band
            ),
            [reach],
            stride=self.stride,
            num_points=self.num_points,
            oneside=self.one_sided,
            twophase=self.two_phase,
            maxvalued=self.take_maximum_value,
        )
        self.elapsed += time.perf_counter() - t0
        return float(achisq[0])

    def evaluate_trigger(self, batch, j):
        """Hook called by _AutoStage in SignalConsistencyChain."""
        # Retrieve underlying RatioBlockInputs from bound method batch.template
        inputs = getattr(batch.template, "__self__", None)
        if inputs is None:
            return None
        local_id = int(inputs.ids_local[j])
        t_idx = int(inputs.t_idxs[j])
        return self.value(
            inputs.stilde,
            inputs.coarse_idx,
            inputs.ref_snr,
            inputs.filters_f,
            local_id,
            t_idx,
        )
