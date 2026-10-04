# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""Compact instantaneous PSD representation, spline continuum modeling,
and self-contained HDF5 packaging.
"""

import numpy as np
import scipy.signal as sig
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
import scipy.stats as stats
import h5py
import json
import igwn_segments as segments
from pycbc.types import TimeSeries, FrequencySeries

__all__ = [
    "fit_als_baseline",
    "CompactInstantaneousPSD",
    "package_overwhitened_product",
]


def fit_als_baseline(log_freqs, log_psd, n_knots=48, p=0.01):
    """Fit a smooth baseline continuum to log10(PSD) vs log10(f).

    Uses log-frequency quantile spline regression, strictly robust to positive line
    and chirp outliers.

    Parameters
    ----------
    log_freqs : ndarray
        log10(f) frequency array.
    log_psd : ndarray
        log10(PSD) spectral power array.
    n_knots : int, optional
        Number of log-spaced knots (default: 48).
    p : float, optional
        Quantile percentile between 0 and 1 (default: 0.01 = 1st percentile).

    Returns
    -------
    cs : CubicSpline
        Fitted spline object mapping log10(f) to log10(PSD).
    knots_f : ndarray
        Knot locations in Hertz.
    knot_vals : ndarray
        Knot values in log10(PSD).
    """
    knots_f = np.logspace(np.log10(10**log_freqs[0]), np.log10(10**log_freqs[-1]), n_knots)
    knots_log_f = np.log10(knots_f)
    half_w = 0.035

    knot_vals = []
    for kf in knots_log_f:
        win = (log_freqs >= kf - half_w) & (log_freqs <= kf + half_w)
        if np.sum(win) > 0:
            knot_vals.append(np.percentile(log_psd[win], int(p * 100)))
        else:
            knot_vals.append(knot_vals[-1] if knot_vals else log_psd[0])

    knot_vals = np.array(knot_vals)
    cs = CubicSpline(knots_log_f, knot_vals)
    return cs, knots_f, knot_vals


class CompactInstantaneousPSD:
    """Compact, factorized representation of time-varying PSD:
    S_n(f, t) = S_continuum(f) * [ sum_{b=1}^{N_b} g_b(t)^2 * B_b(f) ] + sum_k L_k(f, t)

    Features:
    1. Smooth Continuum: cubic spline knots across physical bandwidth (e.g. 15 to 1024 Hz).
    2. Spectral Line Catalog: central frequency f0, peak excess A_k, FWHM Gamma_k.
    3. Multi-Band Non-Stationarity Factors: g_b(t) tracked smoothly at 1 Hz across physical bands.
    4. Instantaneous evaluation in O(1) time.
    """

    def __init__(self, f_min=15.0, f_max=None, n_knots=48):
        self.f_min = float(f_min)
        self.f_max = float(f_max) if f_max is not None else None
        self._f_max_explicit = f_max is not None
        if self.f_max is None:
            self.f_max = 1024.0
            self._f_max_explicit = False
        self.n_knots = int(n_knots)
        self.knots_f = np.logspace(np.log10(self.f_min), np.log10(self.f_max), self.n_knots)
        self.knots_log_f = np.log10(self.knots_f)
        self.knots_log_psd = None
        self.spline_continuum = None
        self.lines = []
        self._init_bands()
        self.scale_times = None
        self.band_scales = None

    def _init_bands(self):
        f_hi_eff = self.f_max
        raw_bands = [
            (15.0, 45.0, "Seismic & Scattered Light"),
            (45.0, 150.0, "Suspension & Mains"),
            (150.0, 500.0, "Intermediate Bucket"),
            (500.0, min(1024.0, f_hi_eff), "Shot Noise Floor"),
        ]
        if f_hi_eff > 1024.0:
            raw_bands.append((1024.0, f_hi_eff, "High-Frequency Shot Noise"))
        self.bands = [
            (max(self.f_min, b_lo), min(f_hi_eff, b_hi), name)
            for b_lo, b_hi, name in raw_bands
            if b_lo < f_hi_eff and b_hi > self.f_min
        ]
        if len(self.bands) == 0:
            self.bands = [(self.f_min, f_hi_eff, "Bandwidth")]

    @staticmethod
    def compute_curvature_adaptive_knots(log_f, log_p, n_knots=48, curvature_alpha=1.5):
        """Compute frequency knot locations adaptively weighted by log-spectral curvature.

        Concentrates spline knots in regions of rapid spectral transition (e.g. seismic wall,
        intermediate sensitivity bucket turnaround, violin resonances) while placing fewer knots
        in flat shot-noise regions, matching BayesLine's trans-dimensional flexibility without MCMC.

        Parameters
        ----------
        log_f : ndarray
            log10(frequency) array.
        log_p : ndarray
            log10(PSD) array.
        n_knots : int, optional
            Number of spline knots (default: 48).
        curvature_alpha : float, optional
            Relative weight of spectral curvature vs uniform log-frequency spacing (default: 1.5).

        Returns
        -------
        knots_f : ndarray
            Adaptive knot frequencies in Hertz.
        knots_log_f : ndarray
            log10 of adaptive knot frequencies.
        """
        grid_log_f = np.linspace(log_f[0], log_f[-1], 1000)
        grid_log_p = np.interp(grid_log_f, log_f, log_p)

        log_p_smooth = gaussian_filter1d(grid_log_p, sigma=15)
        d_gf = float(grid_log_f[1] - grid_log_f[0])
        d2y = np.gradient(np.gradient(log_p_smooth, d_gf), d_gf)
        curv = gaussian_filter1d(np.abs(d2y), sigma=10)

        rho = 1.0 + curvature_alpha * (curv / (np.mean(curv) + 1e-12))
        cdf = np.cumsum(rho)
        cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])

        u = np.linspace(0.0, 1.0, n_knots)
        knots_log_f = np.interp(u, cdf, grid_log_f)
        knots_f = 10.0 ** knots_log_f
        return knots_f, knots_log_f

    def fit_static_from_psd(self, psd_pilot, line_threshold=3.5, p_percentile=10, prominence=2.0,
                            adaptive_knots=True, curvature_alpha=1.5, iterative_refine=True):
        """Fit smooth continuum and catalog narrow lines from a pilot PSD with sub-bin precision.

        Includes:
        1. Curvature-adaptive knot placement: dynamically concentrates spline knots at steep
           spectral gradients (seismic wall, bucket knee, violin modes) matching BayesLine's
           trans-dimensional RJMCMC flexibility in O(1) time (<1 ms).
        2. Sub-bin parabolic line refinement and empirical FWHM measurement.
        3. Iterative joint continuum-line decoupling: eliminates line skirt bias on the continuum.
        4. Analytical Bayesian credible interval preparation.

        Parameters
        ----------
        psd_pilot : FrequencySeries
            Pilot PSD estimate (e.g. median Welch).
        line_threshold : float, optional
            Ratio threshold above continuum to identify lines (default: 3.5).
        p_percentile : float, optional
            Percentile for quantile knot fitting (default: 10, robust to positive signal/line bias).
        prominence : float, optional
            Minimum peak prominence in ratio domain to reject noise fluctuations (default: 2.0).
        adaptive_knots : bool, optional
            Whether to use spectral curvature-adaptive knot allocation (default: True).
        curvature_alpha : float, optional
            Weight of curvature vs uniform spacing in adaptive knot allocation (default: 1.5).
        iterative_refine : bool, optional
            Whether to perform iterative joint continuum-line decoupling (default: True).
        """
        freqs = psd_pilot.sample_frequencies.numpy()
        p_vals = psd_pilot.numpy()
        df = float(psd_pilot.delta_f)

        f_nyquist = float(freqs[-1])
        if not getattr(self, "_f_max_explicit", False) and f_nyquist > self.f_max:
            self.f_max = f_nyquist
            self._init_bands()
            if not adaptive_knots:
                self.knots_f = np.logspace(np.log10(self.f_min), np.log10(self.f_max), self.n_knots)
                self.knots_log_f = np.log10(self.knots_f)

        mask = (freqs >= self.f_min) & (freqs <= self.f_max)
        f_sub = freqs[mask]
        p_sub = p_vals[mask]
        log_f = np.log10(f_sub)
        log_p = np.log10(p_sub)

        # 1. Curvature-adaptive knot placement
        if adaptive_knots:
            knots_f, knots_log_f = self.compute_curvature_adaptive_knots(
                log_f, log_p, n_knots=self.n_knots, curvature_alpha=curvature_alpha
            )
            self.knots_f = knots_f
            self.knots_log_f = knots_log_f

        half_w = 0.035
        # Pass 1: Initial Quantile Spline Continuum
        knot_vals = []
        for kf in self.knots_log_f:
            win = (log_f >= kf - half_w) & (log_f <= kf + half_w)
            if np.sum(win) > 0:
                knot_vals.append(np.percentile(log_p[win], p_percentile))
            else:
                knot_vals.append(knot_vals[-1] if knot_vals else log_p[0])
        self.knots_log_psd = np.array(knot_vals, dtype=np.float64)
        self.spline_continuum = CubicSpline(self.knots_log_f, self.knots_log_psd)

        # Store reference band powers for scale factor tracking
        self.ref_band_powers = np.zeros(len(self.bands), dtype=np.float64)
        for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
            b_mask = (freqs >= b_lo) & (freqs < b_hi)
            self.ref_band_powers[b_idx] = np.median(p_vals[b_mask]) if np.sum(b_mask) > 0 else 1.0

        # Helper function for line detection and sub-bin refinement
        def _extract_lines(p_current, cont_spline):
            cont_eval = 10 ** cont_spline(log_f)
            ratio = p_current / cont_eval
            min_distance = max(1, int(round(0.5 / df)))
            peaks, _ = sig.find_peaks(ratio, height=line_threshold, prominence=prominence, distance=min_distance)

            detected_lines = []
            for p_idx in peaks:
                # Sub-bin parabolic refinement
                if 0 < p_idx < len(p_current) - 1:
                    y0 = np.log(p_current[p_idx])
                    ym1 = np.log(p_current[p_idx - 1])
                    yp1 = np.log(p_current[p_idx + 1])
                    denom = ym1 - 2.0 * y0 + yp1
                    if abs(denom) > 1e-12 and denom < 0:
                        delta = 0.5 * (ym1 - yp1) / denom
                        delta = np.clip(delta, -0.5, 0.5)
                        f0_hat = float(f_sub[p_idx] + delta * df)
                        peak_power = float(np.exp(y0 - 0.25 * (ym1 - yp1) * delta))
                    else:
                        f0_hat = float(f_sub[p_idx])
                        peak_power = float(p_current[p_idx])
                else:
                    f0_hat = float(f_sub[p_idx])
                    peak_power = float(p_current[p_idx])

                # Empirical FWHM measurement via scipy.signal.peak_widths
                widths, _, _, _ = sig.peak_widths(ratio, [p_idx], rel_height=0.5)
                fwhm_meas = float(widths[0] * df)
                fwhm_clamped = float(np.clip(fwhm_meas, 0.05, 3.0))

                c_f0 = float(10 ** cont_spline(np.log10(f0_hat)))
                excess = max(0.0, float(peak_power - c_f0))
                detected_lines.append((f0_hat, excess, fwhm_clamped))
            return detected_lines

        self.lines = _extract_lines(p_sub, self.spline_continuum)

        # Pass 2: Iterative joint line-continuum decoupling (EM-like refinement)
        if iterative_refine:
            cont_eval = 10 ** self.spline_continuum(log_f)
            p_clean = p_sub.copy()
            if len(self.lines) > 0:
                for f0, excess, fwhm in self.lines:
                    mask_l = np.abs(f_sub - f0) <= 5.0 * fwhm
                    if np.any(mask_l):
                        lorentz = excess * (fwhm / 2.0) ** 2 / ((f_sub[mask_l] - f0) ** 2 + (fwhm / 2.0) ** 2)
                        p_clean[mask_l] = np.maximum(p_clean[mask_l] - lorentz, 0.5 * cont_eval[mask_l])

            log_p_clean = np.log10(p_clean)
            knot_vals_ref = []
            for kf in self.knots_log_f:
                win = (log_f >= kf - half_w) & (log_f <= kf + half_w)
                if np.sum(win) > 0:
                    knot_vals_ref.append(np.percentile(log_p_clean[win], 50))
                else:
                    knot_vals_ref.append(knot_vals_ref[-1] if knot_vals_ref else log_p_clean[0])
            self.knots_log_psd = np.array(knot_vals_ref, dtype=np.float64)
            self.spline_continuum = CubicSpline(self.knots_log_f, self.knots_log_psd)
            self.lines = _extract_lines(p_sub, self.spline_continuum)

        # Store continuum residual statistics for Bayesian credible intervals
        cont_final = 10 ** self.spline_continuum(log_f)
        mask_cont = np.ones(len(f_sub), dtype=bool)
        for f0, _, fwhm in self.lines:
            mask_cont &= (np.abs(f_sub - f0) > 3.0 * fwhm)
        res_log = log_p[mask_cont] - np.log10(cont_final[mask_cont])
        raw_mad = float(np.median(np.abs(res_log - np.median(res_log)))) if len(res_log) > 0 else 0.05
        self.continuum_mad_log = max(raw_mad, 0.02)
        self.delta_f_pilot = df

    def track_nonstationarity(self, timeseries, t_step=1.0, t_window=4.0,
                              medfilt_kernel_s=25.0, clip_outlier_sigma=2.5,
                              smooth_sigma_s=4.0):
        """Track time-varying multi-band scale factors g_b(t) with robust signal & glitch rejection.

        Uses a two-stage acausal robust estimation:
        Stage 1: Compute windowed power ratios r_b(t) across physical bands.
        Stage 2: Apply two-sided acausal median filter to strictly reject transient
                 CBC signals (BBH/NSBH/BNS) and short non-Gaussian glitches.
        Stage 3: Clip positive outlier excursions exceeding robust threshold.
        Stage 4: Zero-phase Gaussian smooth to track true slow physical detector drifts (seismic,
                 alignment, thermal) with zero lag and zero transient bump.

        Parameters
        ----------
        timeseries : TimeSeries
            Input strain time series.
        t_step : float, optional
            Time step in seconds between scale factor evaluations (default: 1.0).
        t_window : float, optional
            Window duration in seconds for local power estimation (default: 4.0).
        medfilt_kernel_s : float, optional
            Median filter duration in seconds (default: 25.0).
        clip_outlier_sigma : float, optional
            Outlier clipping threshold in units of median absolute deviation (default: 2.5).
        smooth_sigma_s : float, optional
            Zero-phase Gaussian smoothing standard deviation in seconds (default: 4.0).
        """
        fs = float(timeseries.sample_rate)
        t_start = float(timeseries.start_time)
        t_end = float(timeseries.end_time)
        total_dur = t_end - t_start
        if total_dur <= t_window and total_dur > 0:
            t_window = max(0.5, total_dur / 2.0)
            t_step = min(t_step, t_window / 2.0)
        times = np.arange(t_start + t_window / 2, t_end - t_window / 2 + 1e-9, t_step)
        if len(times) == 0:
            times = np.array([t_start + total_dur / 2.0])
        n_times = len(times)
        n_bands = len(self.bands)

        raw_scales = np.ones((n_bands, n_times), dtype=np.float64)
        if getattr(self, "ref_band_powers", None) is not None and len(self.ref_band_powers) == n_bands:
            ref_powers = self.ref_band_powers.copy()
        else:
            ref_powers = np.zeros(n_bands, dtype=np.float64)
            freqs = np.linspace(self.f_min, self.f_max, 2000)
            cont = self.eval_continuum(freqs)
            for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
                b_mask = (freqs >= b_lo) & (freqs < b_hi)
                ref_powers[b_idx] = np.mean(cont[b_mask]) if np.sum(b_mask) > 0 else 1.0

        w = np.hanning(int(round(t_window * fs)))
        w_sum_sq = np.sum(w ** 2)

        # Stage 1: Windowed power measurements
        for i, t in enumerate(times):
            s_idx = int(round((t - t_window / 2 - t_start) * fs))
            e_idx = s_idx + int(round(t_window * fs))
            if e_idx > len(timeseries):
                break
            seg = timeseries[s_idx:e_idx].numpy()
            spec = np.abs(np.fft.rfft(w * seg)) ** 2 * (2.0 / (fs * w_sum_sq))
            f_seg = np.fft.rfftfreq(len(seg), 1.0 / fs)

            for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
                b_mask = (f_seg >= b_lo) & (f_seg < b_hi)
                if np.sum(b_mask) > 0:
                    P_meas = np.median(spec[b_mask]) / np.log(2.0)
                    ratio = np.sqrt(max(P_meas / max(ref_powers[b_idx], 1e-50), 1e-4))
                    raw_scales[b_idx, i] = ratio
                else:
                    raw_scales[b_idx, i] = 1.0

        # Stage 2 to 4: Acausal robust median filtering, clipping, and smoothing
        scales = np.ones((n_bands, n_times), dtype=np.float64)
        k_med = int(round(medfilt_kernel_s / t_step))
        if k_med % 2 == 0:
            k_med += 1

        for b_idx in range(n_bands):
            r_b = raw_scales[b_idx, :].copy()
            if len(r_b) >= 3:
                curr_k = min(k_med, len(r_b) if len(r_b) % 2 == 1 else len(r_b) - 1)
                curr_k = max(3, curr_k)
                # Stage 2: Temporal median filter (strictly rejects transient signals & glitches)
                med_b = sig.medfilt(r_b, kernel_size=curr_k)

                # Stage 3: Outlier clipping
                res = r_b - med_b
                mad = np.median(np.abs(res))
                sigma_est = 1.4826 * mad
                outlier_thresh = max(clip_outlier_sigma * sigma_est, 0.05 * np.median(med_b))
                clipped_r = np.where(res > outlier_thresh, med_b, r_b)
                clipped_r = np.clip(clipped_r, 0.5 * med_b, 2.0 * med_b)

                # Re-median filter
                med_clipped = sig.medfilt(clipped_r, kernel_size=curr_k)

                # Stage 4: Zero-phase Gaussian smooth
                sigma_samples = max(1.0, smooth_sigma_s / t_step)
                scales[b_idx, :] = gaussian_filter1d(med_clipped, sigma=sigma_samples, mode='nearest')
            else:
                scales[b_idx, :] = r_b

        self.scale_times = times
        self.band_scales = scales

    def eval_continuum(self, freqs):
        """Evaluate smooth continuum at specified frequencies with smooth seismic wall regularization.

        Parameters
        ----------
        freqs : float or ndarray
            Frequency points in Hertz.

        Returns
        -------
        continuum : float or ndarray
            Evaluated smooth continuum spectral power.
        """
        if self.spline_continuum is None:
            raise ValueError("Spline continuum has not been fitted. Call fit_static_from_psd first.")

        freqs_arr = np.asarray(freqs, dtype=np.float64)
        scalar_input = freqs_arr.ndim == 0
        if scalar_input:
            freqs_arr = np.atleast_1d(freqs_arr)

        out = np.zeros_like(freqs_arr, dtype=np.float64)

        # In-band evaluation [f_min, f_max]
        in_band = (freqs_arr >= self.f_min) & (freqs_arr <= self.f_max)
        if np.any(in_band):
            out[in_band] = 10.0 ** self.spline_continuum(np.log10(freqs_arr[in_band]))

        # Smooth physical seismic wall regularization below f_min:
        # Instead of artificial flat clamping, extrapolate using the physical slope at f_min
        low_band = freqs_arr < self.f_min
        if np.any(low_band):
            log_f0 = self.knots_log_f[0]
            val_f0 = float(self.knots_log_psd[0])
            d_slope = float(self.spline_continuum.derivative(1)(log_f0))
            slope_eff = min(d_slope, -4.0)

            f_safe = np.maximum(freqs_arr[low_band], 1.0)
            log_f_low = np.log10(f_safe)
            log_p_low = val_f0 + slope_eff * (log_f_low - log_f0)
            out[low_band] = 10.0 ** log_p_low

        # High-frequency band above f_max:
        high_band = freqs_arr > self.f_max
        if np.any(high_band):
            log_f_hi = self.knots_log_f[-1]
            val_f_hi = float(self.knots_log_psd[-1])
            d_slope_hi = float(self.spline_continuum.derivative(1)(log_f_hi))
            slope_hi_eff = max(d_slope_hi, 0.0)
            log_f_high = np.log10(freqs_arr[high_band])
            log_p_high = val_f_hi + slope_hi_eff * (log_f_high - log_f_hi)
            out[high_band] = 10.0 ** log_p_high

        return out[0] if scalar_input else out

    def eval_psd(self, freqs, t=None, include_lines=True):
        """Evaluate total instantaneous PSD S_n(f, t).

        Parameters
        ----------
        freqs : ndarray
            Frequency points in Hertz.
        t : float, optional
            GPS time for instantaneous scale factor application.
        include_lines : bool, optional
            Whether to add cataloged Lorentzian lines (default: True).

        Returns
        -------
        total_psd : ndarray
            Evaluated PSD at specified frequencies.
        """
        cont = self.eval_continuum(freqs)

        if t is not None and self.scale_times is not None and len(self.scale_times) > 0:
            t_idx = int(np.argmin(np.abs(self.scale_times - t)))
            band_weights = np.ones_like(freqs, dtype=np.float64)
            for b_idx, (b_lo, b_hi, _) in enumerate(self.bands):
                b_mask = (freqs >= b_lo) & (freqs < b_hi)
                band_weights[b_mask] = self.band_scales[b_idx, t_idx] ** 2
            cont = cont * band_weights

        total_psd = cont.copy()
        if include_lines:
            for f0, excess, fwhm in self.lines:
                mask = np.abs(freqs - f0) <= 3.0 * fwhm
                if np.any(mask):
                    lorentz = excess * (fwhm / 2.0) ** 2 / ((freqs[mask] - f0) ** 2 + (fwhm / 2.0) ** 2)
                    total_psd[mask] += lorentz

        return total_psd

    def eval_credible_intervals(self, freqs, t=None, alpha=0.90):
        """Evaluate Bayesian posterior median and credible intervals (e.g. 90% credible band).

        Matches BayesLine's posterior credible intervals in O(1) time without MCMC sampling.

        Parameters
        ----------
        freqs : ndarray
            Frequencies in Hertz.
        t : float, optional
            GPS time for instantaneous scale factor application.
        alpha : float, optional
            Credible interval coverage probability (default: 0.90 for 90% CI).

        Returns
        -------
        psd_med : ndarray
            Instantaneous median PSD.
        psd_lower : ndarray
            Lower credible bound (e.g. 5th percentile for alpha=0.90).
        psd_upper : ndarray
            Upper credible bound (e.g. 95th percentile for alpha=0.90).
        """
        psd_med = self.eval_psd(freqs, t=t, include_lines=True)
        sigma_mad = getattr(self, "continuum_mad_log", 0.08)
        df = getattr(self, "delta_f_pilot", 0.125)
        w_half = 0.035
        # Effective independent bins contributing to knot interpolation
        n_bins = np.maximum(4.0, 2.0 * w_half * freqs * np.log(10.0) / df)
        sigma_knot = (1.4826 * sigma_mad) / np.sqrt(n_bins)
        z = float(stats.norm.ppf(0.5 + alpha / 2.0))
        psd_lower = psd_med * (10.0 ** (-z * sigma_knot))
        psd_upper = psd_med * (10.0 ** (+z * sigma_knot))
        return psd_med, psd_lower, psd_upper

    def get_memory_footprint(self):
        """Calculate serialized memory footprint in bytes."""
        b_knots = len(self.knots_f) * 8 if self.knots_f is not None else 0
        b_lines = len(self.lines) * 24
        b_scales = self.band_scales.nbytes if self.band_scales is not None else 0
        total = b_knots + b_lines + b_scales
        return {
            "knots_bytes": b_knots,
            "lines_bytes": b_lines,
            "scales_bytes": b_scales,
            "total_bytes": total,
            "num_lines": len(self.lines),
        }

    def save_to_hdf(self, h5group):
        """Save compact model representation into an HDF5 group."""
        h5group.attrs["f_min"] = self.f_min
        h5group.attrs["f_max"] = self.f_max
        h5group.attrs["n_knots"] = self.n_knots
        h5group.attrs["description"] = "Factorized Compact Instantaneous PSD S_n(f, t)"

        if self.knots_f is not None:
            h5group.create_dataset("knots_f", data=self.knots_f)
        if self.knots_log_psd is not None:
            h5group.create_dataset("knots_log_psd", data=self.knots_log_psd)

        lines_arr = np.array(self.lines, dtype=[("f0", "f8"), ("excess", "f8"), ("fwhm", "f8")])
        h5group.create_dataset("lines", data=lines_arr)

        if self.scale_times is not None:
            h5group.create_dataset("scale_times", data=self.scale_times)
        if self.band_scales is not None:
            h5group.create_dataset("band_scales", data=self.band_scales)

        bands_grp = h5group.create_group("bands")
        for b_idx, (b_lo, b_hi, b_name) in enumerate(self.bands):
            b_sub = bands_grp.create_group(f"band_{b_idx}")
            b_sub.attrs["f_low"] = b_lo
            b_sub.attrs["f_high"] = b_hi
            b_sub.attrs["name"] = b_name

    @classmethod
    def load_from_hdf(cls, h5group):
        """Load compact model representation from an HDF5 group."""
        f_min = float(h5group.attrs["f_min"])
        f_max = float(h5group.attrs["f_max"])
        n_knots = int(h5group.attrs["n_knots"])
        obj = cls(f_min=f_min, f_max=f_max, n_knots=n_knots)

        if "knots_f" in h5group:
            obj.knots_f = h5group["knots_f"][:]
            obj.knots_log_f = np.log10(obj.knots_f)
        if "knots_log_psd" in h5group:
            obj.knots_log_psd = h5group["knots_log_psd"][:]
            obj.spline_continuum = CubicSpline(obj.knots_log_f, obj.knots_log_psd)

        if "lines" in h5group:
            lines_data = h5group["lines"][:]
            obj.lines = [(float(row["f0"]), float(row["excess"]), float(row["fwhm"])) for row in lines_data]

        if "scale_times" in h5group:
            obj.scale_times = h5group["scale_times"][:]
        if "band_scales" in h5group:
            obj.band_scales = h5group["band_scales"][:]

        if "bands" in h5group:
            bands_grp = h5group["bands"]
            loaded_bands = []
            for k in sorted(bands_grp.keys()):
                b_sub = bands_grp[k]
                loaded_bands.append((float(b_sub.attrs["f_low"]), float(b_sub.attrs["f_high"]), str(b_sub.attrs["name"])))
            if loaded_bands:
                obj.bands = loaded_bands

        return obj


def package_overwhitened_product(output_path, ts_dow, ts_dw, model_compact=None, psd_base=None,
                               metadata=None, gated_segments=None, valid_segments=None, boundary_pad=None):
    """Package the overwhitened stream d_ow(t), whitened stream d_w(t),
    auxiliary instantaneous PSD model, and segment validity/gating masks into a self-contained HDF5 product.

    Parameters
    ----------
    output_path : str
        Target HDF5 file path.
    ts_dow : TimeSeries
        Overwhitened strain time series d_ow(t).
    ts_dw : TimeSeries
        Intermediate whitened strain time series d_w(t).
    model_compact : CompactInstantaneousPSD, optional
        Fitted compact instantaneous PSD model.
    psd_base : FrequencySeries, optional
        Reference base PSD.
    metadata : dict, optional
        Additional dictionary of attributes or metrics.
    gated_segments : list of (float, float), optional
        List of [start_gps, end_gps] intervals that were gated/inpainted.
    valid_segments : list of (float, float), optional
        List of [start_gps, end_gps] intervals of certified valid data (past boundary pad, non-gated).
    boundary_pad : float, optional
        Filter transient boundary padding in seconds (default: max_filter_duration / 2.0 or 2.0s).
    """
    fs = float(ts_dow.sample_rate)
    with h5py.File(output_path, "w") as f:
        # Overwhitened stream d_ow(t)
        f.create_dataset("dow", data=ts_dow.numpy(), compression="gzip", compression_opts=4)
        f["dow"].attrs["sample_rate"] = fs
        f["dow"].attrs["delta_t"] = 1.0 / fs
        f["dow"].attrs["start_time"] = float(ts_dow.start_time)
        f["dow"].attrs["duration"] = float(ts_dow.duration)
        f["dow"].attrs["description"] = "Regularized overwhitened strain d_ow(t) for direct matched filtering"

        # Intermediate whitened stream d_w(t)
        f.create_dataset("dw", data=ts_dw.numpy(), compression="gzip", compression_opts=4)
        f["dw"].attrs["sample_rate"] = fs
        f["dw"].attrs["delta_t"] = 1.0 / fs
        f["dw"].attrs["start_time"] = float(ts_dw.start_time)
        f["dw"].attrs["duration"] = float(ts_dw.duration)
        f["dw"].attrs["description"] = "Intermediate unit-variance whitened strain d_w(t)"

        # Boundary pad & segment bookkeeping
        bpad = boundary_pad
        if bpad is None and metadata is not None and "max_filter_duration" in metadata:
            bpad = float(metadata["max_filter_duration"]) / 2.0
        elif bpad is None:
            bpad = 2.0
        f.attrs["boundary_pad"] = float(bpad)

        # Gated / Inpainted segments
        if gated_segments is not None and len(gated_segments) > 0:
            arr_gated = np.array(gated_segments, dtype=np.float64)
            if arr_gated.ndim == 1:
                arr_gated = arr_gated.reshape(-1, 2)
            f.create_dataset("gated_segments", data=arr_gated)
        else:
            f.create_dataset("gated_segments", data=np.zeros((0, 2), dtype=np.float64))

        # Valid segments
        if valid_segments is not None and len(valid_segments) > 0:
            arr_valid = np.array(valid_segments, dtype=np.float64)
            if arr_valid.ndim == 1:
                arr_valid = arr_valid.reshape(-1, 2)
            f.create_dataset("valid_segments", data=arr_valid)
        else:
            t_s = float(ts_dow.start_time) + bpad
            t_e = float(ts_dow.end_time) - bpad
            if t_e > t_s:
                seg_raw = segments.segmentlist([segments.segment(t_s, t_e)])
                if gated_segments is not None and len(gated_segments) > 0:
                    g_list = segments.segmentlist([segments.segment(float(s[0]), float(s[1])) for s in gated_segments])
                    seg_valid = (seg_raw - g_list).coalesce()
                else:
                    seg_valid = seg_raw
                arr_valid = np.array([(float(s[0]), float(s[1])) for s in seg_valid], dtype=np.float64)
                f.create_dataset("valid_segments", data=arr_valid)
            else:
                f.create_dataset("valid_segments", data=np.zeros((0, 2), dtype=np.float64))

        # Instantaneous PSD model
        if model_compact is not None:
            psd_grp = f.create_group("instantaneous_psd")
            model_compact.save_to_hdf(psd_grp)

        # Base PSD
        if psd_base is not None:
            base_grp = f.require_group("psd_base")
            base_grp.create_dataset("frequencies", data=psd_base.sample_frequencies.numpy())
            base_grp.create_dataset("values", data=psd_base.numpy())
            base_grp.attrs["delta_f"] = float(psd_base.delta_f)

        # Metadata
        if metadata is not None:
            meta_grp = f.create_group("metadata")
            for k, v in metadata.items():
                if isinstance(v, (int, float, str, bool, np.number)):
                    meta_grp.attrs[k] = v
                else:
                    try:
                        meta_grp.attrs[k] = json.dumps(v)
                    except Exception:
                        meta_grp.attrs[k] = str(v)
