# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Parametric Autoregressive (Burg Method) PSD Model."""

import numpy as np
import scipy.signal as sig
from pycbc.types import TimeSeries, FrequencySeries
import pycbc.psd
from .base import BasePSDModel


def fit_burg_ar(x: np.ndarray, order: int):
    """Fit Autoregressive model coefficients using Burg's maximum entropy method.

    Guarantees stable all-pole filter with all roots strictly inside unit circle.

    Parameters
    ----------
    x : ndarray
        1D input signal.
    order : int
        AR model order P.

    Returns
    -------
    a : ndarray
        AR polynomial coefficients [1.0, a_1, a_2, ..., a_P].
    sigma2 : float
        Innovation / driving white noise variance.
    """
    N = len(x)
    f = x.copy().astype(np.float64)
    b = x.copy().astype(np.float64)
    sigma2 = np.mean(x ** 2)

    a = np.array([1.0], dtype=np.float64)

    for m in range(1, order + 1):
        num = 2.0 * np.sum(f[m:] * b[m-1:-1])
        den = np.sum(f[m:] ** 2 + b[m-1:-1] ** 2)
        if den <= 0:
            break
        k = num / den
        k = np.clip(k, -0.9999, 0.9999)

        # Update polynomial
        a_new = np.zeros(m + 1, dtype=np.float64)
        a_new[0] = 1.0
        a_new[1:m] = a[1:m] - k * a[m-1:0:-1]
        a_new[m] = -k
        a = a_new

        sigma2 *= (1.0 - k ** 2)

        # Update forward and backward prediction errors
        f_next = f[m:] - k * b[m-1:-1]
        b_next = b[m-1:-1] - k * f[m:]
        f[m:] = f_next
        b[m-1:-1] = b_next

    return a, sigma2


class AutoregressivePSDModel(BasePSDModel):
    """Parametric all-pole AR PSD model with line cataloging."""

    def __init__(self, f_low=20.0, f_high=None, ar_order=96, line_threshold=3.0):
        super().__init__(f_low=f_low, f_high=f_high)
        self.ar_order = int(ar_order)
        self.line_threshold = float(line_threshold)
        self.a_coeffs = None
        self.sigma2 = None
        self.lines = []
        self.gain_times = None
        self.gain_scales = None
        self._cached_freqs_len = 0
        self._cached_f0 = None
        self._cached_f1 = None
        self._cached_continuum = None

    def fit(self, timeseries: TimeSeries, f_low: float = None, f_high: float = None):
        if f_low is not None:
            self.f_low = float(f_low)
        self.fs = float(timeseries.sample_rate)
        if f_high is not None:
            self.f_high = float(f_high)
        elif self.f_high is None:
            self.f_high = 0.9 * (self.fs / 2.0)

        self.t_start = float(timeseries.start_time)
        self.t_end = float(timeseries.end_time)

        x = timeseries.numpy().astype(np.float64)

        # 1. Fit global AR model using Burg's algorithm
        self.a_coeffs, self.sigma2 = fit_burg_ar(x, self.ar_order)

        # 2. Extract narrow spectral lines against the AR continuum
        # Compare 16s Welch against AR smooth spectrum
        seg_len_pilot = int(round(16.0 * self.fs))
        seg_len_pilot = min(seg_len_pilot, len(timeseries) // 2)
        welch_pilot = pycbc.psd.welch(timeseries, seg_len=seg_len_pilot,
                                      seg_stride=seg_len_pilot // 2, avg_method='median')
        freqs_w = welch_pilot.sample_frequencies.numpy()
        p_w = welch_pilot.numpy()

        p_ar = self._eval_ar_continuum(freqs_w)
        ratio = p_w / np.maximum(p_ar, 1e-60)
        df_w = float(welch_pilot.delta_f)

        mask_f = (freqs_w >= self.f_low) & (freqs_w <= self.f_high)
        peaks, _ = sig.find_peaks(ratio[mask_f], height=self.line_threshold,
                                  distance=max(1, int(round(0.4 / df_w))))

        self.lines = []
        f_sub = freqs_w[mask_f]
        p_sub = p_w[mask_f]
        ar_sub = p_ar[mask_f]

        for p_idx in peaks:
            f0 = float(f_sub[p_idx])
            excess = float(p_sub[p_idx] - ar_sub[p_idx])
            fwhm = max(2.0 * df_w, 0.15)
            self.lines.append((f0, excess, fwhm))

        # 3. Dynamic RMS gain tracking across time (4s window)
        t_step = 1.0
        t_win = 4.0
        N_win = int(round(t_win * self.fs))
        times = np.arange(self.t_start + t_win / 2, self.t_end - t_win / 2 + 1e-9, t_step)
        if len(times) == 0:
            times = np.array([self.t_start + (self.t_end - self.t_start) / 2.0])

        rms_gains = []
        base_rms = np.std(x)
        w = np.hanning(N_win)
        w_norm = np.sqrt(np.mean(w ** 2))

        for t_val in times:
            s_idx = int(round((t_val - t_win / 2 - self.t_start) * self.fs))
            e_idx = s_idx + N_win
            if e_idx > len(x):
                break
            seg_rms = np.std(x[s_idx:e_idx] * w) / (w_norm + 1e-12)
            rms_gains.append(seg_rms / (base_rms + 1e-12))

        self.gain_times = times[:len(rms_gains)]
        self.gain_scales = np.array(rms_gains)

    def _eval_ar_continuum(self, freqs: np.ndarray) -> np.ndarray:
        """Evaluate AR power spectrum: S(f) = (2 * sigma2 / fs) / |A(e^j2pi f / fs)|^2."""
        if (self._cached_freqs_len == len(freqs) and
            self._cached_continuum is not None and
            len(freqs) > 0 and self._cached_f0 == freqs[0] and self._cached_f1 == freqs[-1]):
            return self._cached_continuum.copy()

        freqs = np.asarray(freqs, dtype=np.float64)
        z = np.exp(-1j * 2.0 * np.pi * freqs / self.fs)
        # Evaluate polynomial A(z) = 1 + a_1 z + a_2 z^2 + ...
        # using Horner's method
        val = np.zeros(len(freqs), dtype=np.complex128)
        for c in self.a_coeffs:
            val = val * z + c
        denom = np.abs(val) ** 2
        psd = (2.0 * self.sigma2 / self.fs) / np.maximum(denom, 1e-60)

        self._cached_freqs_len = len(freqs)
        self._cached_f0 = freqs[0] if len(freqs) > 0 else 0
        self._cached_f1 = freqs[-1] if len(freqs) > 0 else 0
        self._cached_continuum = psd.copy()
        return psd

    def eval_psd(self, freqs: np.ndarray, t: float = None) -> np.ndarray:
        freqs = np.asarray(freqs, dtype=np.float64)
        psd = self._eval_ar_continuum(freqs)

        if t is not None and self.gain_times is not None and len(self.gain_times) > 0:
            t_idx = int(np.clip(np.argmin(np.abs(self.gain_times - t)), 0, len(self.gain_times) - 1))
            gain = self.gain_scales[t_idx]
            psd = psd * (gain ** 2)

        # Add lines using searchsorted
        if len(self.lines) > 0:
            is_sorted = len(freqs) > 1 and freqs[1] > freqs[0]
            for f0, excess, fwhm in self.lines:
                w3 = 3.0 * fwhm
                if is_sorted:
                    i_lo = int(np.searchsorted(freqs, f0 - w3))
                    i_hi = int(np.searchsorted(freqs, f0 + w3))
                    if i_hi > i_lo:
                        f_sub = freqs[i_lo:i_hi]
                        lorentz = excess * (fwhm / 2.0) ** 2 / ((f_sub - f0) ** 2 + (fwhm / 2.0) ** 2)
                        psd[i_lo:i_hi] += lorentz
                else:
                    mask = np.abs(freqs - f0) <= w3
                    if np.any(mask):
                        lorentz = excess * (fwhm / 2.0) ** 2 / ((freqs[mask] - f0) ** 2 + (fwhm / 2.0) ** 2)
                        psd[mask] += lorentz

        # Low-frequency brick wall below f_low
        psd[freqs < self.f_low] = 1e30
        return psd
