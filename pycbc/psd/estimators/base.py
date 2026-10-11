# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Base class for modern, non-stationary PSD estimators and continuous overwhitening."""

from abc import ABC, abstractmethod
import numpy as np
from pycbc.types import TimeSeries, FrequencySeries, zeros, complex64
import pycbc.psd


class BasePSDModel(ABC):
    """Abstract base class for full-job PSD estimation, dynamic tracking,
    and continuous overwhitening.
    """

    def __init__(self, f_low=20.0, f_high=None):
        self.f_low = float(f_low)
        self.f_high = float(f_high) if f_high is not None else None
        self.fs = None
        self.t_start = None
        self.t_end = None

    @abstractmethod
    def fit(self, timeseries: TimeSeries, f_low: float = None, f_high: float = None):
        """Fit the PSD model across the continuous timeseries of the job."""
        pass

    @abstractmethod
    def eval_psd(self, freqs: np.ndarray, t: float = None) -> np.ndarray:
        """Evaluate the instantaneous PSD S_n(f, t) in strain^2/Hz at given frequencies.
        If t is None, returns the baseline/time-average spectrum.
        """
        pass

    def average_psd(self, t_start: float, t_end: float, delta_f: float, f_len: int,
                    t_step: float = 2.0) -> FrequencySeries:
        """Compute the robust segment-average operating point over [t_start, t_end].

        Uses a time-slice median across evaluation timestamps in [t_start, t_end]
        to strictly reject transient glitches and non-Gaussian bursts.
        Applies low-frequency zero-inversion cutoff (S_n -> inf below f_low).

        Parameters
        ----------
        t_start : float
            GPS start time of the segment.
        t_end : float
            GPS end time of the segment.
        delta_f : float
            Frequency spacing in Hertz.
        f_len : int
            Number of frequency bins.
        t_step : float, optional
            Step size in seconds between time slices (default: 2.0s).

        Returns
        -------
        psd_fs : FrequencySeries
            The segment-median PSD as a PyCBC FrequencySeries.
        """
        freq_grid = np.arange(f_len, dtype=np.float64) * delta_f
        t_start_c = max(t_start, self.t_start) if self.t_start is not None else t_start
        t_end_c = min(t_end, self.t_end) if self.t_end is not None else t_end
        
        t_samples = np.arange(t_start_c, t_end_c + 1e-9, max(0.5, t_step))
        if len(t_samples) == 0:
            t_samples = np.array([(t_start + t_end) / 2.0])

        psd_matrix = np.zeros((len(t_samples), f_len), dtype=np.float64)
        for idx, t_val in enumerate(t_samples):
            psd_matrix[idx, :] = self.eval_psd(freq_grid, t=t_val)

        # Time-slice median across the segment
        seg_median_psd = np.median(psd_matrix, axis=0)

        # Apply low-frequency cutoff (S_n -> inf below f_low)
        k_low = int(np.floor(self.f_low / delta_f))
        if k_low > 0:
            seg_median_psd[:k_low] = 1e30

        psd_fs = FrequencySeries(seg_median_psd.astype(np.float64), delta_f=delta_f, epoch=t_start)
        return psd_fs

    def construct_fir_kernel(self, t: float = None, max_filter_duration: float = 4.0,
                             delta_f: float = 0.0625) -> np.ndarray:
        """Construct zero-phase symmetric FIR overwhitening kernel taps.

        Parameters
        ----------
        t : float, optional
            GPS time for instantaneous PSD evaluation.
        max_filter_duration : float, optional
            Maximum filter duration in seconds (default: 4.0).
        delta_f : float, optional
            Frequency resolution for kernel construction (default: 0.0625 Hz = 16s FFT).

        Returns
        -------
        fir_dow : ndarray
            1D symmetric FIR impulse response with odd length.
        """
        fs = self.fs if self.fs is not None else 2048.0
        f_nyq = fs / 2.0
        f_len = int(round(f_nyq / delta_f)) + 1
        freq_grid = np.arange(f_len, dtype=np.float64) * delta_f

        psd_vals = self.eval_psd(freq_grid, t=t)
        # Ensure positive and non-zero
        psd_vals = np.maximum(psd_vals, 1e-60)
        psd_fs = FrequencySeries(psd_vals.astype(np.float64), delta_f=delta_f)

        # Truncate inverse spectrum with Hann window
        N_filt = int(round(max_filter_duration * fs))
        if N_filt % 2 == 0:
            N_filt += 1

        psd_ist = pycbc.psd.inverse_spectrum_truncation(
            psd_fs, max_filter_len=N_filt, which_spectrum='invpsd',
            low_frequency_cutoff=self.f_low, trunc_method='hann'
        )

        # Build FIR impulse response
        # Create zero-phase symmetric kernel via IFFT
        N_fft = (f_len - 1) * 2
        inv_psd = np.zeros(f_len, dtype=np.float64)
        k_min = int(np.floor(self.f_low / delta_f))
        inv_psd[k_min:] = 1.0 / np.maximum(psd_ist.numpy()[k_min:], 1e-60)

        # Smooth high-frequency taper if f_high is set
        if self.f_high is not None and self.f_high < f_nyq:
            k_high = int(np.floor(self.f_high / delta_f))
            taper_w = max(1, int(round(0.05 * self.f_high / delta_f)))
            if k_high < f_len:
                taper = 0.5 * (1.0 + np.cos(np.linspace(0, np.pi, taper_w)))
                k_taper_end = min(f_len, k_high + taper_w)
                actual_w = k_taper_end - k_high
                inv_psd[k_high:k_taper_end] *= taper[:actual_w]
                inv_psd[k_taper_end:] = 0.0

        # Symmetrize IFFT to center peak
        h_full = np.roll(np.fft.irfft(inv_psd, n=N_fft), N_fft // 2)
        k_mid = N_fft // 2
        half_filt = N_filt // 2
        fir_dow = h_full[k_mid - half_filt : k_mid + half_filt + 1]
        return fir_dow

    def overwhiten(self, timeseries: TimeSeries, max_filter_duration: float = 4.0,
                   chunk_len_s: float = 16.0, crossfade_s: float = 2.0) -> TimeSeries:
        """Streamingly overwhiten the continuous timeseries using overlap-save
        with smoothly updated FIR filters and Hann cross-fading to eliminate seam artifacts.

        Parameters
        ----------
        timeseries : TimeSeries
            Input continuous strain time series.
        max_filter_duration : float, optional
            Duration of FIR filter in seconds (default: 4.0s).
        chunk_len_s : float, optional
            Chunk size in seconds for updating time-varying FIR filters (default: 16.0s).
        crossfade_s : float, optional
            Duration in seconds of Hann cross-fade between adjacent chunk filters (default: 2.0s).

        Returns
        -------
        out_ts : TimeSeries
            Continuously overwhitened strain time series.
        """
        fs = float(timeseries.sample_rate)
        N = len(timeseries)
        dt = 1.0 / fs
        t_start = float(timeseries.start_time)

        M = int(round(max_filter_duration * fs))
        if M % 2 == 0:
            M += 1
        k_center = (M - 1) // 2

        L = int(round(chunk_len_s * fs))
        fade_len = int(round(crossfade_s * fs))
        fade_len = min(fade_len, L // 4)

        out = np.zeros(N, dtype=np.float64)
        weight = np.zeros(N, dtype=np.float64)
        x = timeseries.numpy().astype(np.float64)

        from scipy.fft import next_fast_len

        # Process in chunks with step L - fade_len
        step = L - fade_len if fade_len > 0 else L

        for s_idx in range(0, N, step):
            e_idx = min(s_idx + L, N)
            chunk_dur = e_idx - s_idx
            t_chunk_center = t_start + (s_idx + e_idx) * 0.5 * dt

            # Compute FIR filter for this chunk timestamp
            fir = self.construct_fir_kernel(t=t_chunk_center, max_filter_duration=max_filter_duration)

            # Overlap-save block convolution
            cs = s_idx - k_center
            ce = e_idx + (M - 1 - k_center)
            pad_left = max(0, -cs)
            pad_right = max(0, ce - N)

            sub_x = x[max(0, cs) : min(N, ce)]
            if pad_left > 0 or pad_right > 0:
                pad_mode = 'reflect' if len(sub_x) > max(pad_left, pad_right) else 'edge'
                sub_x = np.pad(sub_x, (pad_left, pad_right), mode=pad_mode)

            n_fft = next_fast_len(len(sub_x) + M)
            sub_fd = np.fft.rfft(sub_x, n=n_fft)
            fir_fd = np.fft.rfft(fir, n=n_fft)
            conv_res = np.fft.irfft(sub_fd * fir_fd, n=n_fft)[M - 1 : M - 1 + chunk_dur]

            # Hann window for smooth cross-fading
            w_chunk = np.ones(chunk_dur, dtype=np.float64)
            if s_idx > 0 and fade_len > 0:
                # Fade in at left edge
                w_chunk[:fade_len] = 0.5 * (1.0 - np.cos(np.linspace(0, np.pi, fade_len)))
            if e_idx < N and fade_len > 0:
                # Fade out at right edge
                w_chunk[-fade_len:] = 0.5 * (1.0 + np.cos(np.linspace(0, np.pi, fade_len)))

            out[s_idx:e_idx] += conv_res * w_chunk
            weight[s_idx:e_idx] += w_chunk

        # Normalize overlapping weights
        valid = weight > 0
        out[valid] /= weight[valid]

        out_ts = TimeSeries(out.astype(np.float64), delta_t=dt, epoch=timeseries.start_time)
        out_ts._is_overwhitened = True
        return out_ts
