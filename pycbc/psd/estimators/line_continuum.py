# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Two-component Line-Continuum Spectral Model.

Models the noise power spectral density as a smooth, shape-preserving PCHIP spline continuum
plus a sum of Cauchy-Lorentzian spectral lines, following the physical two-component
spectral decomposition of gravitational-wave detector noise.
"""

import numpy as np
import scipy.signal as sig
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d
from pycbc.types import TimeSeries, FrequencySeries
import pycbc.psd
from .base import BasePSDModel


class LineContinuumPSDModel(BasePSDModel):
    """Two-component Line-Continuum spectral model.
    
    Decomposes the PSD into:
        S_n(f) = S_spline(f) + sum_i [ A_i * gamma_i^2 / ((f - f_i)^2 + gamma_i^2) ]
    
    where S_spline(f) is a monotonic shape-preserving Hermite spline (PCHIP) across
    log-spaced control knots, and each spectral line is parameterized by central
    frequency f_i, half-width at half-maximum (HWHM) gamma_i, and peak amplitude A_i.
    Guarantees strict lower-envelope continuum preservation and eliminates inverted drops.
    """

    def __init__(self, f_low=20.0, f_high=None, n_knots=64, pilot_seg_len_s=16.0,
                 line_snr_threshold=3.5, max_lines=400):
        super().__init__(f_low=f_low, f_high=f_high)
        self.n_knots = int(n_knots)
        self.pilot_seg_len_s = float(pilot_seg_len_s)
        self.line_snr_threshold = float(line_snr_threshold)
        self.max_lines = int(max_lines)

        self.knots_f = None
        self.knots_val = None
        self.spline = None
        self.lines = []  # list of tuples: (f0, A, gamma)
        self.pilot_psd = None

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

        # 1. Pilot high-resolution Welch spectrum (16s segments, 50% overlap)
        seg_len_samples = int(round(self.pilot_seg_len_s * self.fs))
        seg_len_samples = min(seg_len_samples, len(timeseries) // 2)
        if seg_len_samples <= 0:
            seg_len_samples = len(timeseries)

        pilot_fs = pycbc.psd.welch(
            timeseries, seg_len=seg_len_samples, seg_stride=max(1, seg_len_samples // 2),
            avg_method='median'
        )
        self.pilot_psd = pilot_fs

        freqs = pilot_fs.sample_frequencies.numpy().astype(np.float64)
        p_vals = pilot_fs.numpy().astype(np.float64)
        df = float(pilot_fs.delta_f)

        mask = (freqs >= max(10.0, self.f_low - 5.0)) & (freqs <= self.f_high)
        f_sub = freqs[mask]
        p_sub = p_vals[mask]
        p_sub = np.maximum(np.nan_to_num(p_sub, nan=1e-46, posinf=1e-30, neginf=1e-60), 1e-60)

        # 2. Line extraction via smooth baseline ratio
        smooth_p = 10.0 ** gaussian_filter1d(np.log10(p_sub), sigma=max(1, int(round(2.5 / df))))
        excess_ratio = p_sub / np.maximum(smooth_p, 1e-60)

        peaks, _ = sig.find_peaks(excess_ratio, height=self.line_snr_threshold,
                                  distance=max(1, int(round(0.3 / df))))

        if len(peaks) > 0:
            order = np.argsort(excess_ratio[peaks])[::-1]
            peaks = peaks[order[:self.max_lines]]

        self.lines = []
        is_line = np.zeros(len(f_sub), dtype=bool)

        for p_idx in peaks:
            f0 = float(f_sub[p_idx])
            amp = float(p_sub[p_idx] - smooth_p[p_idx])
            if amp <= 0:
                continue
            gamma = max(1.5 * df, 0.12)
            self.lines.append((f0, amp, gamma))
            k_lo = max(0, p_idx - int(round(3.0 * gamma / df)))
            k_hi = min(len(f_sub), p_idx + int(round(3.0 * gamma / df)) + 1)
            is_line[k_lo:k_hi] = True

        # 3. Fit Shape-Preserving Spline Continuum on Line-Cleaned Data
        clean_f = f_sub[~is_line]
        clean_p = p_sub[~is_line]

        if len(clean_f) < self.n_knots:
            clean_f = f_sub
            clean_p = p_sub

        log_clean_f = np.log10(clean_f)
        log_clean_p = np.log10(np.maximum(clean_p, 1e-60))

        # Logarithmically spaced knots
        knots_log_f = np.linspace(log_clean_f[0], log_clean_f[-1], self.n_knots)
        knots_f = 10.0 ** knots_log_f
        half_w = 0.5 * (knots_log_f[1] - knots_log_f[0])

        knots_log_p = []
        for kf in knots_log_f:
            win = (log_clean_f >= kf - half_w) & (log_clean_f <= kf + half_w)
            if np.sum(win) > 2:
                knots_log_p.append(np.median(log_clean_p[win]))
            else:
                knots_log_p.append(np.interp(kf, log_clean_f, log_clean_p))

        self.knots_f = knots_f
        self.knots_val = np.nan_to_num(np.array(knots_log_p, dtype=np.float64),
                                       nan=-46.0, posinf=-30.0, neginf=-60.0)
        # Use PchipInterpolator (Piecewise Cubic Hermite Interpolating Polynomial)
        # to strictly preserve monotonicity and prevent Runge oscillations / spline sagging
        self.spline = PchipInterpolator(knots_log_f, self.knots_val)

    def eval_psd(self, freqs: np.ndarray, t: float = None, include_lines: bool = True, **kwargs) -> np.ndarray:
        freqs = np.asarray(freqs, dtype=np.float64)
        valid_f = np.clip(freqs, self.knots_f[0], self.knots_f[-1])
        cont = 10.0 ** self.spline(np.log10(valid_f))

        # Add Lorentzian lines with smooth taper on flanks to eliminate truncation step discontinuities
        if include_lines and len(self.lines) > 0:
            is_sorted = len(freqs) > 1 and freqs[1] > freqs[0]
            for f0, amp, gamma in self.lines:
                w_taper = max(10.0 * gamma, 1.0)
                w_cut = max(20.0 * gamma, 2.0)
                if is_sorted:
                    i_lo = int(np.searchsorted(freqs, f0 - w_cut))
                    i_hi = int(np.searchsorted(freqs, f0 + w_cut))
                    if i_hi > i_lo:
                        f_sub = freqs[i_lo:i_hi]
                        df_line = np.abs(f_sub - f0)
                        lorentz = amp * (gamma ** 2) / (df_line ** 2 + gamma ** 2)
                        mid_wing = df_line > w_taper
                        if np.any(mid_wing):
                            taper = 0.5 * (1.0 + np.cos(np.pi * (df_line[mid_wing] - w_taper) / (w_cut - w_taper)))
                            lorentz[mid_wing] *= taper
                        cont[i_lo:i_hi] += lorentz
                else:
                    df_line = np.abs(freqs - f0)
                    mask = df_line <= w_cut
                    if np.any(mask):
                        df_sub = df_line[mask]
                        lorentz = amp * (gamma ** 2) / (df_sub ** 2 + gamma ** 2)
                        mid_wing = df_sub > w_taper
                        if np.any(mid_wing):
                            taper = 0.5 * (1.0 + np.cos(np.pi * (df_sub[mid_wing] - w_taper) / (w_cut - w_taper)))
                            lorentz[mid_wing] *= taper
                        cont[mask] += lorentz

        # Low-frequency and high-frequency cutoff
        cont[freqs < self.f_low] = 1e30
        if self.f_high is not None:
            cont[freqs > self.f_high] = 1e30
        return cont
