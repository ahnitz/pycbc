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
                 high_frequency_cutoff=None,
                 tap_sample_rate=2048, engine_sample_rate=2048,
                 engine='matchedfilter-hierarchical', false_dismissal=1e-3,
                 coarse_band_hz=0, first_stage_snr=0, peak_window=None, **kwargs):
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
        # Granularity of the primary peak search, in seconds: one peak is reported per span this
        # long (pycbc passes its per-template cluster window), so how matchedfilter blocks the
        # series is its own choice and does not change which peaks survive clustering.
        self.peak_binsize = int(round(float(peak_window) * self.engine_sr)) if peak_window else None

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
        env_band = os.environ.get('PYCBC_RATIO_BAND_HZ')
        if env_band:
            parts = [float(x.strip()) for x in env_band.split(',') if x.strip()]
            self._coarse_band_hz = tuple(int(x) for x in parts) if len(parts) > 1 else parts[0]
        elif isinstance(coarse_band_hz, (tuple, list)):
            self._coarse_band_hz = tuple(int(x) for x in coarse_band_hz)
        else:
            self._coarse_band_hz = float(coarse_band_hz)
        self._device = os.environ.get('PYCBC_RATIO_DEVICE') or None
        self._first_stage_snr = float(os.environ.get(
            'PYCBC_RATIO_FIRST_STAGE_SNR', first_stage_snr))
        self._td_bank = None
        self._ap_ref_key = None

    def prepare_filters(self, fir_taps, tap_counts):
        """Build the fine-stage TimeDomainFilterBank for a batch of taps.

        Nothing is returned: with the peak granularity stated, the bank settles its
        block sizes at its first reference, so its spectra do not exist yet (the
        chisq has its own, from prepare_chisq_filters).
        """
        fft_lengths = None
        if os.environ.get('PYCBC_RATIO_FFT_LENGTH'):
            fft_lengths = [int(x.strip()) for x in os.environ['PYCBC_RATIO_FFT_LENGTH'].split(',') if x.strip()]
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
            fft_lengths=fft_lengths,
            # the peak granularity: with it stated, the block size is matchedfilter's choice
            binsize=self.peak_binsize,
        )

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
        # scale is a float64 scalar: multiplying into a complex64 output keeps the
        # double-precision product and its rounding, without a complex128 temporary.
        out = np.empty(len(snr), dtype=np.complex64)
        np.multiply(snr.numpy(), scale, out=out, casting='same_kind')
        return out, h_norm

    def prepare_chisq_filters(self, fir_taps, tap_counts, block_length, reach=0):
        """Template spectra on the chisq block: length block_length, the convention of filters_f.

        The chisq is computed on a block of its own, centred on each trigger, so its
        resolution is pycbc's choice and independent of how matchedfilter blocks the
        series for filtering. The block must hold every filter with room for the lags
        evaluated around the trigger (`reach` either side, for the auto-chisq);
        otherwise this raises rather than lengthen the block, since its length is part
        of the statistic.
        """
        n = int(block_length)
        bank = _mf.TimeDomainFilterBank(
            fir_taps, tap_counts=np.asarray(tap_counts, dtype=np.int64),
            tap_sample_rate=self.tap_sr, data_sample_rate=self.engine_sr,
            engine='matchedfilter', fft_lengths=[n])
        for g in bank.groups:
            centre = n // 2
            lo, hi = g['c_bad'], g['c_bad'] + g['n_valid']
            if g['n'] != n or centre - reach < lo or centre + reach >= hi:
                raise ValueError(
                    "--chisq-block-length %d is too short for this template group: its "
                    "longest ratio filter needs valid lags [%d, %d) around the block centre "
                    "+-%d" % (n, lo, hi, reach))
        self._chisq_bank = bank
        return bank.filters_f

    @property
    def chisq_bank(self):
        """The flat bank holding the chisq spectra (block layout and valid margins)."""
        return self._chisq_bank

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
        if self._engine_mode in ('hier', 'check') and getattr(self, '_ap_ref_key', None) != id(profile):
            if self._td_bank is not None:
                self._td_bank.set_reference_from_template(stilde, psd, profile, f_high=self.f_high)
                self._ap_ref_key = id(profile)

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

        res = self._td_bank.filter_series(self.ref_snr, windows=valid_slice, binsize=self.peak_binsize)
        self.ref_snr = None
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


def get_dynamic_snr_renorm_factor(series, dt=None, window_duration=8.0, hollow_duration=0.5, variance_floor=1.0, scale=None, max_boost_factor=None):
    """Compute rolling window hollow dynamic SNR renormalization envelope.

    Downweights non-stationary noise excursions and glitch tails while preserving
    short gravitational-wave signal peaks by excluding a central hollow window.
    Supports two-sided renormalization when variance_floor < 1.0, allowing
    effective SNR enhancement during pristine quiet periods up to max_boost_factor.

    Parameters
    ----------
    series : ndarray or Array or TimeSeries
        Complex SNR time series (or scaled reference series if scale is provided).
    dt : float, optional
        Sampling interval in seconds of the series. If None, inferred from series.delta_t
        if available, else defaults to 1.0.
    window_duration : float, optional
        Full width of the outer rolling window (seconds, default: 8.0s).
    hollow_duration : float, optional
        Half-width of the inner hollow exclusion window (seconds, default: 0.5s,
        giving an excluded region of +/- 0.5s around each sample).
    variance_floor : float, optional
        Minimum variance threshold to renormalize (default: 1.0). When set < 1.0
        (e.g. 0.5 or 0.25), allows effective SNR enhancement during quiet periods.
    scale : float, optional
        Amplitude scaling factor if series is an internally scaled buffer where
        z(t) = series(t) / scale. If None, series is assumed to already be a
        complex SNR series where Var(Re(z)) = Var(Im(z)) = 1.0 in stationary Gaussian noise.
    max_boost_factor : float, optional
        Safety ceiling for the maximum allowable SNR boost factor (default: None,
        meaning bounded only by 1.0 / sqrt(variance_floor)). Recommended: 1.25–1.50.

    Returns
    -------
    renorm_factor : ndarray (float32 or float64)
        1D array of renormalization scale factors in (0.0, max_boost_factor or 1/sqrt(variance_floor)].
    """
    factor = DynamicSNRRenormFactor(
        series, dt=dt, window_duration=window_duration, hollow_duration=hollow_duration,
        variance_floor=variance_floor, scale=scale, max_boost_factor=max_boost_factor)
    return factor[np.arange(len(factor))]


class DynamicSNRRenormFactor(object):
    """The envelope of get_dynamic_snr_renorm_factor, evaluated where it is read.

    A search reads the envelope only at its triggers, a few hundred samples of a
    series of ~10^6; building it everywhere was most of the reference stage's
    time. This keeps the one O(n) part (the running sum of the local power) and
    evaluates the rest at the requested indices with the same arithmetic, so the
    values are identical. Supports len() and indexing by integer arrays.
    """

    def __init__(self, series, dt=None, window_duration=8.0, hollow_duration=0.5,
                 variance_floor=1.0, scale=None, max_boost_factor=None):
        self._n = n = len(series)
        self._ones = (os.environ.get('PYCBC_DISABLE_DYNAMIC_SNR_RENORM', '0') == '1' or n == 0)
        if dt is None:
            dt = float(getattr(series, 'delta_t', 1.0))
        else:
            dt = float(dt)
        if dt <= 0:
            self._ones = True
        if self._ones:
            self._dtype = np.float32
            return

        if window_duration <= 2.0 * hollow_duration:
            raise ValueError(
                f"window_duration ({window_duration}s) must be greater than 2 * hollow_duration ({2.0 * hollow_duration}s)"
            )

        arr = np.asarray(series)
        # Local power P(t) = 0.5 * |z(t)|^2
        power = arr.real * arr.real
        power += arr.imag * arr.imag
        if scale is not None and float(scale) > 0:
            inv_scale = 1.0 / float(scale)
            power *= (0.5 * inv_scale * inv_scale)
        else:
            power *= 0.5

        np.nan_to_num(power, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

        self._w_outer = int(round((0.5 * window_duration) / dt))
        self._w_inner = int(round(hollow_duration / dt))
        self._cumsum = np.empty(n + 1, dtype=np.float64)
        self._cumsum[0] = 0.0
        np.cumsum(power, dtype=np.float64, out=self._cumsum[1:])
        self._floor = max(float(variance_floor), 1e-12)
        self._max_boost = max_boost_factor
        self._dtype = np.float64 if arr.dtype == np.complex128 else np.float32

    def __len__(self):
        return self._n

    def __getitem__(self, idx):
        idx = np.asarray(idx)
        if self._ones:
            return np.ones(idx.shape, dtype=self._dtype)
        n, cumsum = self._n, self._cumsum
        l_out = np.clip(idx - self._w_outer, 0, n)
        r_out = np.clip(idx + self._w_outer + 1, 0, n)
        l_in = np.clip(idx - self._w_inner, 0, n)
        r_in = np.clip(idx + self._w_inner + 1, 0, n)

        sum_hollow = (cumsum[r_out] - cumsum[l_out]) - (cumsum[r_in] - cumsum[l_in])
        cnt_hollow = (r_out - l_out) - (r_in - l_in)

        var_est = sum_hollow / np.maximum(cnt_hollow, 1)
        var_eff = np.where(np.isfinite(var_est), np.maximum(var_est, self._floor), self._floor)
        renorm_factor = np.where(np.isfinite(var_eff), 1.0 / np.sqrt(var_eff), 1.0)
        if self._max_boost is not None:
            renorm_factor = np.clip(renorm_factor, 0.0, float(self._max_boost))
        return renorm_factor.astype(self._dtype)


def dynamic_snr_renormalize(series, dt=None, window_duration=8.0, hollow_duration=0.5, variance_floor=1.0, scale=None, max_boost_factor=None):
    """Dynamically renormalize SNR series using a rolling window with hollow exclusion.

    Downweights non-stationary noise excursions and glitch tails while preserving
    short gravitational-wave signal peaks by excluding a central hollow window.
    Supports two-sided renormalization when variance_floor < 1.0 and max_boost_factor is set.

    Parameters
    ----------
    series : ndarray or Array or TimeSeries
        Complex SNR time series (or scaled reference series if scale is provided).
    dt : float, optional
        Sampling interval (seconds) of the series. If None, inferred from series.delta_t if available.
    window_duration : float, optional
        Full width of the outer rolling window (seconds, default: 8.0s).
    hollow_duration : float, optional
        Half-width of the inner hollow exclusion window (seconds, default: 0.5s,
        giving an excluded region of +/- 0.5s around each sample).
    variance_floor : float, optional
        Minimum variance threshold to renormalize (default: 1.0).
    scale : float, optional
        Amplitude scaling factor if series is an internally scaled buffer.
    max_boost_factor : float, optional
        Safety ceiling for maximum allowable SNR boost factor (default: None).

    Returns
    -------
    renormalized_series : same type as series
        Series scaled by renorm factor.
    """
    if os.environ.get('PYCBC_DISABLE_DYNAMIC_SNR_RENORM', '0') == '1':
        return series

    n = len(series)
    if n == 0:
        return series

    factor = get_dynamic_snr_renorm_factor(
        series, dt=dt,
        window_duration=window_duration,
        hollow_duration=hollow_duration,
        variance_floor=variance_floor,
        scale=scale,
        max_boost_factor=max_boost_factor,
    )
    return series * factor


