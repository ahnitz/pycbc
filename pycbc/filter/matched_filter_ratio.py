import os
import numpy as np

from pycbc.fft import IFFT
from pycbc.filter.matchedfilter import (
    correlate, get_cutoff_indices, matched_filter_core,
)
from pycbc.types import Array, complex64, zeros

# matchedfilter back end is required for ratio/FIR filtering.
try:
    import matchedfilter as _mf
except ImportError:
    _mf = None


def matchedfilter_available():
    """True if the matchedfilter back end can be used."""
    return _mf is not None


class MatchedFilterRatioControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    Acts as the PyCBC adapter to the C++/AVX-optimized matchedfilter.TimeDomainFilterBank.
    """

    def __init__(self, snr_threshold, delta_f,
                 high_frequency_cutoff=None, batch_size=64,
                 tap_sample_rate=2048, engine_sample_rate=2048,
                 engine='matchedfilter-hierarchical', false_dismissal=1e-3,
                 coarse_band_hz=0, first_stage_snr=0, **kwargs):
        if _mf is None:
            raise ImportError(
                "matchedfilter is required for ratio/FIR filtering but is not installed"
            )

        self.delta_f = delta_f
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff

        self.tap_sr = int(tap_sample_rate)
        self.engine_sr = int(engine_sample_rate)

        self.threshold_sq = float(snr_threshold**2)
        self.batch_size = batch_size

        # Taps are generated at tap_sample_rate but filtered against data at
        # engine_sample_rate; the ratio must be an exact integer.
        exact_ratio = self.tap_sr / self.engine_sr
        self.decimation_factor = int(np.round(exact_ratio))
        if abs(exact_ratio - self.decimation_factor) > 1e-5 or self.decimation_factor < 1:
            raise ValueError(
                f"Multi-rate Error: The bank sample rate ({self.tap_sr} Hz) must be "
                f"an exact integer multiple of the engine sample "
                f"rate ({self.engine_sr} Hz).\n"
                f"Calculated ratio was {exact_ratio:.4f}. Please use an engine "
                f"sample rate that evenly divides the bank sample rate."
            )

        self._ref_plan = None
        self._ref_direct = os.environ.get('PYCBC_RATIO_REFDIRECT', '1') != '0'

        # matchedfilter backend mode: 'hier' or 'flat'
        mode = {'matchedfilter': 'flat',
                'matchedfilter-hierarchical': 'hier'}.get(engine, engine)
        self._engine_mode = os.environ.get('PYCBC_RATIO_ENGINE', mode).lower()
        self._engine_fd = float(os.environ.get(
            'PYCBC_RATIO_ENGINE_FD', false_dismissal))
        self._coarse_band_hz = float(os.environ.get('PYCBC_RATIO_BAND_HZ', coarse_band_hz))
        self._device = os.environ.get('PYCBC_RATIO_DEVICE') or None
        self._first_stage_snr = float(os.environ.get(
            'PYCBC_RATIO_FIRST_STAGE_SNR', first_stage_snr))
        self._td_bank = None
        self._ap_ref_key = None

    def prepare_filters(self, fir_taps, tap_counts):
        """
        Prepare frequency-domain filters for a batch of taps using TimeDomainFilterBank.
        """
        self._td_bank = _mf.TimeDomainFilterBank(
            fir_taps, tap_counts,
            tap_sample_rate=self.tap_sr,
            data_sample_rate=self.engine_sr,
            engine=self._engine_mode,
            threshold=self.snr_threshold,
            false_dismissal=self._engine_fd,
            first_stage_snr=self._first_stage_snr,
            coarse_band_hz=self._coarse_band_hz,
            device=self._device,
            max_batch_size=self.batch_size,
        )
        n_taps_max = int(np.max(tap_counts))
        return self._td_bank.filters_f, n_taps_max

    def make_reference_series(self, stilde, psd, ref_template):
        """Return the scaled reference SNR series and its template sigmasq."""
        h_norm = ref_template.sigmasq(psd)
        if self._ref_direct:
            snr, norm = self._reference_snr(stilde, psd, ref_template, h_norm)
        else:
            snr, _, norm = matched_filter_core(
                ref_template, stilde, psd=psd,
                low_frequency_cutoff=ref_template.f_lower,
                high_frequency_cutoff=self.f_high, h_norm=h_norm)
        scale = (norm * stilde.delta_t) / self.decimation_factor
        # The cached IFFT output is overwritten by the next reference; upper
        # series must survive while all of its middle children are processed.
        return np.asarray(snr.numpy(), dtype=np.complex64).copy() * scale, h_norm

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, indices,
                        valid_slice=None, reference_series=None,
                        profile_template=None):
        """
        Process a single data segment through TimeDomainFilterBank.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        if isinstance(ref_template, (int, float, np.floating)):
            h_norm = float(ref_template)
        elif getattr(self, '_cached_hnorm_key', None) == (id(ref_template), id(psd)):
            h_norm = self._cached_hnorm
        else:
            h_norm = self._cached_hnorm = ref_template.sigmasq(psd)
            self._cached_hnorm_key = (id(ref_template), id(psd))

        profile = profile_template if profile_template is not None else ref_template
        if self._engine_mode in ('hier', 'check') and getattr(self, '_ap_ref_key', None) != (id(profile), id(psd)):
            if self._td_bank is not None:
                self._td_bank.set_reference_from_template(stilde, psd, profile, f_high=self.f_high)
                self._ap_ref_key = (id(profile), id(psd))

        if reference_series is None:
            if self._ref_direct:
                snr, norm = self._reference_snr(stilde, psd, ref_template, h_norm)
            else:
                snr, _, norm = matched_filter_core(
                    ref_template, stilde, psd=psd,
                    low_frequency_cutoff=ref_template.f_lower,
                    high_frequency_cutoff=self.f_high, h_norm=h_norm)
            self.ref_snr = snr.numpy()
            self.ref_snr *= (norm * stilde.delta_t) / self.decimation_factor
        else:
            self.ref_snr = reference_series if isinstance(reference_series, np.ndarray) and reference_series.dtype == np.complex64 else np.asarray(reference_series, dtype=np.complex64)
            if len(self.ref_snr) != (len(stilde) - 1) * 2:
                raise ValueError('reference series length does not match data segment')

        res = self._td_bank.filter_series(self.ref_snr, valid_slice=valid_slice)
        local_idxs = res.template_indices
        t_idxs = res.sample_indices
        snr_vals = res.snr
        tstarts = res.block_starts

        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, tstarts, h_norm
        else:
            return [], [], [], tstarts, h_norm

    def _reference_snr(self, stilde, psd, ref_template, h_norm):
        """The reference SNR series, driving the class-based IFFT directly.

        matched_filter_core goes through pycbc.fft's *function* API, which on
        the MKL backend builds a DFTI descriptor, uses it once and frees it on
        every call -- 2 ms of a 3 ms transform at 2^20.  The class API caches
        the descriptor, which is what it is for, and the plan cache this object
        already keeps for the block FFTs serves here too.  The buffers are
        cached with it, so the four per-call allocations go as well.

        Same arithmetic as matched_filter_core, in the same order.
        """
        N = (len(stilde) - 1) * 2
        kmin, kmax = get_cutoff_indices(
            ref_template.f_lower, self.f_high, stilde.delta_f, N)
        plan, qt, q = self._get_ref_plan(N, kmax)
        # Only the low strip is cleared per call.  correlate writes solely
        # [kmin:kmax] and execute() reads qt and writes q, so nothing ever
        # dirties the region above kmax -- and kmax derives from self.f_high,
        # an engine attribute, so it is fixed for the whole run.  Clearing it
        # once with the buffer removes 4 MB of stores per segment.  kmin comes
        # from the template's f_lower and does vary, so that strip stays.
        qt[:kmin] = 0
        # stilde arrives overwhitened (see pycbc_inspiral_fir), so there is
        # no PSD division on this path at all.
        correlate(ref_template[kmin:kmax], stilde[kmin:kmax], qt[kmin:kmax])
        plan.execute()
        norm = (4.0 * stilde.delta_f) / np.sqrt(h_norm)
        return q, norm

    def _get_ref_plan(self, size, kmax):
        """Cached IFFT plan and its buffers for the reference filter.

        Keyed on kmax as well as size: zeros() hands back a cleared buffer,
        which is what lets _reference_snr skip clearing above kmax on every
        call.  If the upper cutoff ever did change, this reallocates rather
        than leaving a stale band in the transform.
        """
        cached = self._ref_plan
        if cached is None or cached[0] != (size, kmax):
            qt = zeros(size, dtype=complex64)
            q = zeros(size, dtype=complex64)
            plan = IFFT(qt, q)
            cached = self._ref_plan = ((size, kmax), plan, qt, q)
        _, plan, qt, q = cached
        return plan, Array(qt, copy=False), q
