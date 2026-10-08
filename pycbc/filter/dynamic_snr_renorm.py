import os
import numpy as np


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
