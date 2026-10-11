# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Two-component Line-Continuum Spline model with unbiased median knots and smooth non-stationarity tracking."""

import numpy as np
import scipy.signal as sig
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.ndimage import gaussian_filter1d, median_filter
from pycbc.types import TimeSeries, FrequencySeries
import pycbc.psd
from .base import BasePSDModel


class SplineLinePSDModel(BasePSDModel):
    """Two-component model separating high-Q spectral lines from an unbiased
    spline continuum, with smooth frequency-blended non-stationarity tracking.
    """

    def __init__(self, f_low=20.0, f_high=None, n_knots=48, pilot_seg_len_s=16.0):
        super().__init__(f_low=f_low, f_high=f_high)
        self.n_knots = int(n_knots)
        self.pilot_seg_len_s = float(pilot_seg_len_s)
        self.knots_f = None
        self.knots_log_f = None
        self.spline_continuum = None
        self.lines = []
        self.bands = []
        self.scale_times = None
        self.band_scales = None

    def _init_bands(self, f_max):
        raw_bands = [
            (self.f_low, 45.0, "Seismic & Low-Frequency"),
            (45.0, 150.0, "Suspension & Mains"),
            (150.0, 500.0, "Intermediate Bucket"),
            (500.0, f_max, "High-Frequency Shot Noise"),
        ]
        self.bands = [
            (max(self.f_low, b_lo), min(f_max, b_hi), name)
            for b_lo, b_hi, name in raw_bands
            if b_lo < f_max and b_hi > self.f_low
        ]

    def _compute_adaptive_knots(self, log_f, log_p, curvature_alpha=1.5):
        grid_log_f = np.linspace(log_f[0], log_f[-1], 1000)
        grid_log_p = np.interp(grid_log_f, log_f, log_p)
        log_p_smooth = gaussian_filter1d(grid_log_p, sigma=15)
        d_gf = float(grid_log_f[1] - grid_log_f[0])
        d2y = np.gradient(np.gradient(log_p_smooth, d_gf), d_gf)
        curv = gaussian_filter1d(np.abs(d2y), sigma=10)

        rho = 1.0 + curvature_alpha * (curv / (np.mean(curv) + 1e-12))
        cdf = np.cumsum(rho)
        cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])

        u = np.linspace(0.0, 1.0, self.n_knots)
        knots_log_f = np.interp(u, cdf, grid_log_f)
        return 10.0 ** knots_log_f, knots_log_f

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
        self._init_bands(self.f_high)

        # 1. High-resolution pilot PSD (16s Welch median)
        seg_len_samples = int(round(self.pilot_seg_len_s * self.fs))
        seg_len_samples = min(seg_len_samples, len(timeseries) // 2)
        pilot_psd = pycbc.psd.welch(
            timeseries, seg_len=seg_len_samples, seg_stride=seg_len_samples // 2,
            avg_method='median'
        )

        freqs = pilot_psd.sample_frequencies.numpy()
        p_vals = pilot_psd.numpy()
        df = float(pilot_psd.delta_f)

        mask = (freqs >= max(10.0, self.f_low - 5.0)) & (freqs <= self.f_high)
        f_sub = freqs[mask]
        p_sub = p_vals[mask]
        log_f = np.log10(f_sub)
        log_p = np.log10(p_sub)

        # 2. Curvature adaptive knots
        self.knots_f, self.knots_log_f = self._compute_adaptive_knots(log_f, log_p)

        # 3. Iterative Line Detection & Unbiased Median Continuum Fitting
        half_w = 0.035
        # Pass 1: Median filter baseline
        p_smooth_initial = 10.0 ** gaussian_filter1d(log_p, sigma=int(round(3.0 / df)))
        line_mask = p_sub > (2.5 * p_smooth_initial)

        # Extract lines
        peaks, props = sig.find_peaks(p_sub / p_smooth_initial, height=2.5, distance=int(round(0.4 / df)))
        self.lines = []
        for p_idx in peaks:
            f0 = float(f_sub[p_idx])
            excess = float(p_sub[p_idx] - p_smooth_initial[p_idx])
            fwhm = max(2.0 * df, 0.15)
            self.lines.append((f0, excess, fwhm))

        # Pass 2: Unbiased Median Knot Fitting (Masking line peaks)
        clean_log_f = log_f[~line_mask]
        clean_log_p = log_p[~line_mask]

        knot_vals = []
        for kf in self.knots_log_f:
            win = (clean_log_f >= kf - half_w) & (clean_log_f <= kf + half_w)
            if np.sum(win) > 5:
                # 50th percentile (median) on line-cleaned data: strictly unbiased!
                knot_vals.append(np.median(clean_log_p[win]))
            else:
                win_all = (log_f >= kf - half_w) & (log_f <= kf + half_w)
                knot_vals.append(np.median(log_p[win_all]) if np.sum(win_all) > 0 else log_p[0])

        self.knots_log_psd = np.array(knot_vals, dtype=np.float64)
        self.spline_continuum = CubicSpline(self.knots_log_f, self.knots_log_psd)

        # 4. Smooth Dynamic Non-Stationarity Tracking
        self._track_nonstationarity(timeseries)

    def _track_nonstationarity(self, timeseries: TimeSeries, t_step=1.0, t_window=4.0):
        t_start = float(timeseries.start_time)
        t_end = float(timeseries.end_time)
        total_dur = t_end - t_start
        times = np.arange(t_start + t_window / 2, t_end - t_window / 2 + 1e-9, t_step)
        if len(times) == 0:
            times = np.array([t_start + total_dur / 2.0])

        n_times = len(times)
        n_bands = len(self.bands)
        fs = float(timeseries.sample_rate)

        # Baseline reference power per band from the fitted continuum
        ref_powers = np.zeros(n_bands, dtype=np.float64)
        f_eval = np.linspace(self.f_low, self.f_high, 2000)
        cont_eval = 10.0 ** self.spline_continuum(np.log10(f_eval))
        for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
            b_m = (f_eval >= b_lo) & (f_eval < b_hi)
            ref_powers[b_idx] = np.mean(cont_eval[b_m]) if np.sum(b_m) > 0 else 1.0

        raw_scales = np.ones((n_bands, n_times), dtype=np.float64)
        w = np.hanning(int(round(t_window * fs)))
        w_sum_sq = np.sum(w ** 2)

        for i, t in enumerate(times):
            s_idx = int(round((t - t_window / 2 - t_start) * fs))
            e_idx = s_idx + len(w)
            if e_idx > len(timeseries):
                break
            seg = timeseries[s_idx:e_idx].numpy()
            spec = np.abs(np.fft.rfft(w * seg)) ** 2 * (2.0 / (fs * w_sum_sq))
            f_seg = np.fft.rfftfreq(len(seg), 1.0 / fs)

            for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
                b_m = (f_seg >= b_lo) & (f_seg < b_hi)
                band_p = np.mean(spec[b_m]) if np.sum(b_m) > 0 else ref_powers[b_idx]
                raw_scales[b_idx, i] = np.sqrt(max(1e-12, band_p / (ref_powers[b_idx] + 1e-60)))

        # Clean and smooth scale factors across time (25s median filter + 4s Gaussian)
        k_med = max(3, int(round(25.0 / t_step)))
        if k_med % 2 == 0:
            k_med += 1

        smoothed_scales = np.ones_like(raw_scales)
        for b_idx in range(n_bands):
            med_s = median_filter(raw_scales[b_idx], size=k_med)
            # Clip outlier spikes
            mad = np.median(np.abs(raw_scales[b_idx] - med_s))
            clipped = np.clip(raw_scales[b_idx], med_s - 3.0 * mad, med_s + 3.0 * mad)
            smoothed_scales[b_idx] = gaussian_filter1d(clipped, sigma=max(1.0, 4.0 / t_step))

        self.scale_times = times
        self.band_scales = smoothed_scales

    def _eval_smooth_band_weights(self, freqs: np.ndarray, t: float = None) -> np.ndarray:
        """Evaluate continuous frequency-blended band scaling factors with monotonic smooth transitions."""
        if t is None or self.scale_times is None or len(self.scale_times) == 0:
            return np.ones_like(freqs, dtype=np.float64)

        t_idx = int(np.clip(np.argmin(np.abs(self.scale_times - t)), 0, len(self.scale_times) - 1))
        scales = self.band_scales[:, t_idx]

        # Band center frequencies in log10 space
        band_centers_log = np.array([
            0.5 * (np.log10(b_lo) + np.log10(b_hi))
            for (b_lo, b_hi, _) in self.bands
        ], dtype=np.float64)
        log_scales_sq = np.log10(np.maximum(scales ** 2, 1e-6))

        # Clamp freqs in log space for flat extrapolation at band endpoints
        log_f = np.log10(np.maximum(freqs, 1.0))
        clamped_log_f = np.clip(log_f, band_centers_log[0], band_centers_log[-1])

        pchip = PchipInterpolator(band_centers_log, log_scales_sq)
        smooth_w = 10.0 ** pchip(clamped_log_f)
        return smooth_w

    def eval_psd(self, freqs: np.ndarray, t: float = None, include_lines: bool = True) -> np.ndarray:
        freqs = np.asarray(freqs, dtype=np.float64)
        valid_f = np.clip(freqs, self.knots_f[0], self.knots_f[-1])
        cont = 10.0 ** self.spline_continuum(np.log10(valid_f))

        # Smooth band weighting
        band_w = self._eval_smooth_band_weights(freqs, t=t)
        total_psd = cont * band_w

        if include_lines and len(self.lines) > 0:
            is_sorted = len(freqs) > 1 and freqs[1] > freqs[0]
            for f0, excess, fwhm in self.lines:
                w3 = 3.0 * fwhm
                if is_sorted:
                    i_lo = int(np.searchsorted(freqs, f0 - w3))
                    i_hi = int(np.searchsorted(freqs, f0 + w3))
                    if i_hi > i_lo:
                        f_sub = freqs[i_lo:i_hi]
                        lorentz = excess * (fwhm / 2.0) ** 2 / ((f_sub - f0) ** 2 + (fwhm / 2.0) ** 2)
                        total_psd[i_lo:i_hi] += lorentz
                else:
                    mask = np.abs(freqs - f0) <= w3
                    if np.any(mask):
                        lorentz = excess * (fwhm / 2.0) ** 2 / ((freqs[mask] - f0) ** 2 + (fwhm / 2.0) ** 2)
                        total_psd[mask] += lorentz

        # Low-frequency wall below f_low
        total_psd[freqs < self.f_low] = 1e30
        return total_psd
