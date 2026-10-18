import numpy as np


# matchedfilter back end is required for ratio/FIR filtering.
try:
    import matchedfilter as _mf
except ImportError:
    _mf = None


def _require_matchedfilter():
    if _mf is None:
        raise ImportError("matchedfilter is required for ratio/FIR filtering but is not installed")


def decimation_factor(tap_sample_rate, data_sample_rate):
    """The integer ratio of the bank's tap sample rate to the data's; raises if it is not one."""
    exact = float(tap_sample_rate) / float(data_sample_rate)
    factor = int(np.round(exact))
    if abs(exact - factor) > 1e-5 or factor < 1:
        raise ValueError(
            f"The bank sample rate ({tap_sample_rate} Hz) must be an exact integer multiple "
            f"of the data sample rate ({data_sample_rate} Hz); the ratio is {exact:.4f}.")
    return factor


def leaf_bank_options(engine, false_dismissal, coarse_band_hz, first_stage_snr, device,
                      peak_window, tap_sample_rate, data_sample_rate):
    """Keyword arguments for a fine bank's matchedfilter TimeDomainFilterBank (the search sets
    the threshold). engine is pycbc's --ratio-filter-engine name; peak_window (s) is the
    granularity of the peak search, one peak a template per span this long (pycbc's cluster
    window), so how matchedfilter blocks the series is its own choice."""
    _require_matchedfilter()
    decimation_factor(tap_sample_rate, data_sample_rate)
    mode = {'matchedfilter': 'flat', 'matchedfilter-hierarchical': 'hier'}.get(engine, engine)
    return dict(
        tap_sample_rate=int(tap_sample_rate),
        data_sample_rate=int(data_sample_rate),
        engine=mode.lower(),
        false_dismissal=float(false_dismissal),
        first_stage_snr=float(first_stage_snr),
        coarse_band_hz=float(coarse_band_hz),
        device=device or None,
        binsize=int(round(float(peak_window) * int(data_sample_rate))) if peak_window else None,
    )


class ChisqFilters(object):
    """A fine group's chisq filters: each template's spectrum on the chisq block (length
    block_length, conjugated, centre tap at lag 0), made from its taps when first read --
    most templates never have a trigger loud enough for a chisq.

    The chisq is computed on a block of its own, centred on each trigger, so its resolution
    is pycbc's choice and independent of how matchedfilter blocks the series for filtering.
    The block must hold every filter with room for the lags evaluated around the trigger
    (reach either side, for the auto-chisq); otherwise this raises rather than lengthen the
    block, since its length is part of the statistic. valid: (T, 2) [start, stop) of each
    filter's valid lags in its block."""

    def __init__(self, taps, counts, block_length, reach=0, tap_sample_rate=2048,
                 data_sample_rate=2048):
        _require_matchedfilter()
        n = self.block_length = int(block_length)
        self._taps = taps
        self._counts = np.asarray(counts, dtype=np.int64)
        self._rates = dict(tap_sample_rate=tap_sample_rate, data_sample_rate=data_sample_rate)
        self.valid = _mf.template_valid(self._counts, n, **self._rates)
        centre = n // 2
        short = np.flatnonzero((centre - reach < self.valid[:, 0])
                               | (centre + reach >= self.valid[:, 1]))
        if len(short):
            lo, hi = self.valid[short[0]]
            raise ValueError(
                "--chisq-block-length %d is too short for this template group: a "
                "ratio filter has valid lags [%d, %d), not the block centre "
                "+-%d" % (n, lo, hi, reach))
        self._rows = {}

    def __len__(self):
        return len(self._counts)

    def __getitem__(self, i):
        i = int(i)
        row = self._rows.get(i)
        if row is None:
            spectra, _ = _mf.template_spectra(self._taps[i:i + 1], self._counts[i:i + 1],
                                              self.block_length, **self._rates)
            row = self._rows[i] = spectra[0]
        return row


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
    series of ~10^6, so the envelope is never built everywhere. The O(n) part is
    the local power (float32, as the series) and its float64 sum over every block
    of _BLOCK samples; a window sum at a requested index is then the block sums
    it spans plus the partial blocks at its edges, so every sample still enters
    at full resolution. Supports len() and indexing by integer arrays.

    A series with non-finite samples takes the power of those samples as zero.
    """

    _BLOCK = 16

    def __init__(self, series, dt=None, window_duration=8.0, hollow_duration=0.5,
                 variance_floor=1.0, scale=None, max_boost_factor=None):
        self._n = n = len(series)
        self._ones = n == 0
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
        # Local power P(t) = 0.5 * |z(t)|^2. For a REAL series x = Re z (the real top filter's
        # output, real-real-complex hierarchy) the same expectation is x^2: E x^2 = E|z|^2 / 2.
        # The constant factor (and the 1/scale^2 of a scaled series) multiplies the window
        # sums, not every sample.
        if np.iscomplexobj(arr):
            power = arr.real * arr.real
            power += arr.imag * arr.imag
            half = 0.5
        else:
            power = arr * arr
            half = 1.0
        if scale is not None and float(scale) > 0:
            self._factor = half / (float(scale) * float(scale))
        else:
            self._factor = half

        blocks = self._block_prefix(power)
        if not np.isfinite(blocks[-1]):
            np.nan_to_num(power, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
            blocks = self._block_prefix(power)
        self._power, self._blocks = power, blocks

        self._w_outer = int(round((0.5 * window_duration) / dt))
        self._w_inner = int(round(hollow_duration / dt))
        self._floor = max(float(variance_floor), 1e-12)
        self._max_boost = max_boost_factor
        self._dtype = np.float64 if arr.dtype == np.complex128 else np.float32

    _PARTIAL = np.arange(_BLOCK - 1)

    @classmethod
    def _block_prefix(cls, power):
        """blocks[j]: the float64 sum of power over samples [0, j * _BLOCK)."""
        B, n = cls._BLOCK, len(power)
        whole = n // B
        blocks = np.zeros(-(-n // B) + 1, dtype=np.float64)
        blocks[1:whole + 1] = power[:whole * B].reshape(whole, B).sum(axis=1, dtype=np.float64)
        if len(blocks) > whole + 1:
            blocks[-1] = power[whole * B:].sum(dtype=np.float64)
        return np.cumsum(blocks, out=blocks)

    def _prefix(self, m):
        """Sum of the power over samples [0, m), for a 1D integer array m in [0, n]."""
        B = self._BLOCK
        base = m - m % B
        v = self._power[np.minimum(base[:, None] + self._PARTIAL, self._n - 1)]
        v[self._PARTIAL >= (m - base)[:, None]] = 0     # the partial block [base, m) only
        return self._blocks[m // B] + v.sum(axis=1, dtype=np.float64)

    def __len__(self):
        return self._n

    def __getitem__(self, idx):
        idx = np.asarray(idx)
        if self._ones:
            return np.ones(idx.shape, dtype=self._dtype)
        n = self._n
        l_out = np.clip(idx - self._w_outer, 0, n)
        r_out = np.clip(idx + self._w_outer + 1, 0, n)
        l_in = np.clip(idx - self._w_inner, 0, n)
        r_in = np.clip(idx + self._w_inner + 1, 0, n)

        edges = self._prefix(np.concatenate([r_out.ravel(), l_out.ravel(), r_in.ravel(), l_in.ravel()]))
        e_ro, e_lo, e_ri, e_li = edges.reshape(4, *idx.shape)
        sum_hollow = self._factor * ((e_ro - e_lo) - (e_ri - e_li))
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


