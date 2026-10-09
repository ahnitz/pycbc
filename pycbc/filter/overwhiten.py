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

"""Core routines for regularized continuous overwhitening and whitening
of gravitational-wave strain data for matched filtering.
"""

import os
import glob
from pathlib import Path
import h5py
import json
import numpy as np
import scipy.signal as sig
import igwn_segments as segments
from pycbc.types import TimeSeries, FrequencySeries, zeros
import pycbc.psd
import pycbc.fft


__all__ = [
    "construct_regularized_kernels",
    "overlap_save_filter",
    "produce_conditioned_streams",
    "overwhiten_strain",
    "overwhiten_block_strain",
    "subtract_coherent_lines",
    "matched_filter_overwhitened",
    "verify_time_varying_streaming_seams",
    "OverwhitenedData",
    "OverwhitenedDataset",
    "OverwhitenedSlice",
]


def construct_regularized_kernels(psd_base, fs, duration, f_low=18.0, f_taper=4.0,
                                  f_high=None, max_filter_duration=4.0):
    """Construct regularized overwhitening and whitening kernels.

    Constructs:
    - C^2 smooth Hann highpass taper to eliminate low-frequency wall sinc-ringing.
    - Inverse Spectrum Truncation (IST) ensuring compact time-domain impulse response <= max_filter_duration.
    - Precomputed centered zero-phase symmetric FIR filters for overlap-save continuous streaming.

    Parameters
    ----------
    psd_base : FrequencySeries
        Base PSD estimate.
    fs : float
        Sampling frequency in Hertz.
    duration : float
        Total duration in seconds of the stream (defines delta_f = 1/duration).
    f_low : float, optional
        Low frequency cutoff in Hertz (default: 18.0).
    f_taper : float, optional
        Width in Hertz of the smooth highpass taper below f_low (default: 4.0).
    f_high : float, optional
        High frequency cutoff in Hertz. Frequencies at or above f_high are zeroed out (default: None).
    max_filter_duration : float, optional
        Maximum time-domain filter impulse response duration in seconds (default: 4.0).

    Returns
    -------
    invpsd_full : ndarray
        Frequency-domain regularized overwhitening kernel values.
    w_kernel_full : ndarray
        Frequency-domain dimensionless whitening kernel values (unit discrete variance).
    fir_dow : ndarray
        Centered zero-phase FIR filter kernel for d_ow(t).
    fir_w : ndarray
        Centered zero-phase FIR filter kernel for d_w(t).
    psd_ist_ow : FrequencySeries
        Truncated inverse PSD frequency series.
    delta_f_full : float
        Frequency resolution of the interpolated full spectrum.
    """
    N_full = int(round(duration * fs))
    delta_f_full = 1.0 / duration

    # Support CompactInstantaneousPSD directly as psd_base or interpolate FrequencySeries
    if hasattr(psd_base, "eval_psd"):
        freqs_eval = np.fft.rfftfreq(N_full, 1.0 / fs)
        psd_vals = psd_base.eval_psd(freqs_eval, include_lines=True)
        psd_full = FrequencySeries(psd_vals, delta_f=delta_f_full)
    else:
        psd_full = pycbc.psd.interpolate(psd_base, delta_f_full)
    freqs_full = psd_full.sample_frequencies.numpy()

    # C^2 smooth highpass taper
    taper_full = np.ones(len(freqs_full), dtype=np.float64)
    f_start = max(0.0, f_low - f_taper)
    taper_full[freqs_full < f_start] = 0.0
    mid = (freqs_full >= f_start) & (freqs_full < f_low)
    if np.any(mid) and f_taper > 0:
        taper_full[mid] = 0.5 * (1.0 - np.cos(np.pi * (freqs_full[mid] - f_start) / f_taper))
    if f_high is None:
        f_high = 0.9 * (fs / 2.0)
    taper_full[freqs_full >= f_high] = 0.0

    # Safe filter truncation bounds: max_filter_duration cannot exceed half stream duration
    max_filt_dur_bound = max(1.0 / fs, min(float(max_filter_duration), (duration / 2.0) - 1.0 / fs))
    N_filt = min(int(round(max_filt_dur_bound * fs)), N_full - 2)
    if N_filt % 2 != 0:
        N_filt += 1

    # Truncate inverse spectrum
    psd_ist_ow = pycbc.psd.inverse_spectrum_truncation(
        psd_full, max_filter_len=N_filt, which_spectrum='invpsd',
        low_frequency_cutoff=f_low, high_frequency_cutoff=f_high,
        trunc_method='hann'
    )
    psd_ist_w = pycbc.psd.inverse_spectrum_truncation(
        psd_full, max_filter_len=N_filt, which_spectrum='invasd',
        low_frequency_cutoff=f_low, high_frequency_cutoff=f_high,
        trunc_method='hann'
    )

    # Overwhitening kernel with zero-division protection
    psd_ow_safe = np.maximum(psd_ist_ow.numpy(), 1e-60)
    invpsd_full = np.where(taper_full > 0, taper_full / psd_ow_safe, 0.0)

    # Dimensionless Whitening kernel for unit discrete variance:
    # Bandwidth delta_F = exact discrete noise-equivalent bandwidth of the regularized taper
    delta_F = float(np.sum(taper_full[1:-1] ** 2) + 0.5 * taper_full[0]**2 + 0.5 * taper_full[-1]**2) * delta_f_full
    if delta_F <= 0:
        raise ValueError(f"Sampling rate {fs} Hz is too low for f_low={f_low} Hz.")
    psd_w_safe = np.maximum(psd_ist_w.numpy(), 1e-60)
    w_kernel_full = np.where(taper_full > 0, taper_full / np.sqrt(delta_F * psd_w_safe), 0.0)

    # Precompute zero-phase symmetric FIR kernels for overlap-save
    h_dow_full = np.roll(np.fft.irfft(invpsd_full, n=N_full), N_full // 2)
    h_w_full = np.roll(np.fft.irfft(w_kernel_full, n=N_full), N_full // 2)

    k_mid = N_full // 2
    half_filt = min(N_filt // 2, k_mid)
    fir_dow = h_dow_full[k_mid - half_filt : k_mid + half_filt + 1]
    fir_w = h_w_full[k_mid - half_filt : k_mid + half_filt + 1]

    return invpsd_full, w_kernel_full, fir_dow, fir_w, psd_ist_ow, delta_f_full


def overlap_save_filter(timeseries, fir_kernel, chunk_len_s=16.0):
    """Continuous overlap-save streaming convolution for arbitrary FIR filter.

    Applies a zero-phase FIR kernel in streaming chunks using the overlap-save
    algorithm, guaranteeing exact linear convolution and eliminating seam/boundary
    artifacts at chunk transitions.

    Parameters
    ----------
    timeseries : TimeSeries
        Input PyCBC TimeSeries.
    fir_kernel : ndarray
        1D FIR filter impulse response (typically symmetric with odd length).
    chunk_len_s : float, optional
        Streaming chunk length in seconds (default: 16.0).

    Returns
    -------
    filtered : TimeSeries
        Output filtered PyCBC TimeSeries with identical delta_t and start_time.
    """
    fs = float(timeseries.sample_rate)
    N = len(timeseries)
    dt = 1.0 / fs
    x = timeseries.numpy()

    from scipy.fft import next_fast_len
    M = len(fir_kernel)
    k_center = (M - 1) // 2
    L = int(round(chunk_len_s * fs))
    out = np.zeros(N, dtype=np.float64)

    n_fft = next_fast_len(L + 2 * M)
    fir_fd = np.fft.rfft(fir_kernel, n=n_fft)

    for s in range(0, N, L):
        e = min(s + L, N)
        block_len = e - s
        cs = s - k_center
        ce = e + (M - 1 - k_center)

        pad_left = max(0, -cs)
        pad_right = max(0, ce - N)
        sub = x[max(0, cs) : min(N, ce)]
        if pad_left > 0 or pad_right > 0:
            pad_mode = 'reflect' if len(sub) > max(pad_left, pad_right) else 'edge'
            sub = np.pad(sub, (pad_left, pad_right), mode=pad_mode)

        sub_fd = np.fft.rfft(sub, n=n_fft)
        c = np.fft.irfft(sub_fd * fir_fd, n=n_fft)[M - 1 : M - 1 + block_len]
        out[s:e] = c

    return TimeSeries(out, delta_t=dt, epoch=timeseries.start_time)


def produce_conditioned_streams(timeseries, invpsd_full, w_kernel_full, fir_dow, fir_w,
                                delta_f_full=None, use_overlap_save=True, chunk_len_s=16.0):
    """Produce both overwhitened stream d_ow(t) and intermediate whitened stream d_w(t).

    Parameters
    ----------
    timeseries : TimeSeries
        Input conditioned strain time series.
    invpsd_full : ndarray
        Full frequency-domain overwhitening kernel.
    w_kernel_full : ndarray
        Full frequency-domain whitening kernel.
    fir_dow : ndarray
        Time-domain FIR filter for overwhitening.
    fir_w : ndarray
        Time-domain FIR filter for whitening.
    delta_f_full : float, optional
        Frequency spacing if use_overlap_save is False.
    use_overlap_save : bool, optional
        If True, use overlap-save FIR streaming convolution (default: True).
        If False, apply full frequency-domain circular multiplication.
    chunk_len_s : float, optional
        Chunk length in seconds for overlap-save (default: 16.0).

    Returns
    -------
    ts_dow : TimeSeries
        Overwhitened strain stream d_ow(t).
    ts_dw : TimeSeries
        Intermediate whitened strain stream d_w(t) certified to standard Gaussian white noise.
    """
    fs = float(timeseries.sample_rate)
    N = len(timeseries)
    dt = 1.0 / fs

    if not use_overlap_save:
        if delta_f_full is None:
            delta_f_full = 1.0 / (N * dt)
        fd = timeseries.to_frequencyseries()
        dow_fd = FrequencySeries(fd.numpy() * invpsd_full, delta_f=delta_f_full, epoch=timeseries.start_time)
        dw_fd = FrequencySeries(fd.numpy() * w_kernel_full, delta_f=delta_f_full, epoch=timeseries.start_time)
        return dow_fd.to_timeseries(), dw_fd.to_timeseries()

    # Overlap-save streaming with shared forward FFT
    from scipy.fft import next_fast_len
    M = len(fir_dow)
    k_center = (M - 1) // 2
    L = int(round(chunk_len_s * fs))
    x = timeseries.numpy()
    d_ow = np.zeros(N, dtype=np.float64)
    d_w = np.zeros(N, dtype=np.float64)

    n_fft = next_fast_len(L + 2 * M)
    fir_dow_fd = np.fft.rfft(fir_dow, n=n_fft)
    fir_w_fd = np.fft.rfft(fir_w, n=n_fft)

    for s in range(0, N, L):
        e = min(s + L, N)
        block_len = e - s
        cs = s - k_center
        ce = e + (M - 1 - k_center)

        pad_left = max(0, -cs)
        pad_right = max(0, ce - N)
        sub = x[max(0, cs) : min(N, ce)]
        if pad_left > 0 or pad_right > 0:
            pad_mode = 'reflect' if len(sub) > max(pad_left, pad_right) else 'edge'
            sub = np.pad(sub, (pad_left, pad_right), mode=pad_mode)

        sub_fd = np.fft.rfft(sub, n=n_fft)
        c_dow = np.fft.irfft(sub_fd * fir_dow_fd, n=n_fft)[M - 1 : M - 1 + block_len]
        c_w = np.fft.irfft(sub_fd * fir_w_fd, n=n_fft)[M - 1 : M - 1 + block_len]

        d_ow[s:e] = c_dow
        d_w[s:e] = c_w

    ts_dow = TimeSeries(d_ow, delta_t=dt, epoch=timeseries.start_time)
    ts_dw = TimeSeries(d_w, delta_t=dt, epoch=timeseries.start_time)
    return ts_dow, ts_dw


def subtract_coherent_lines(timeseries, lines=None, psd=None, f_low=18.0,
                            line_threshold=3.5, max_lines=64, bandwidth=0.2):
    """Two-stage pre-whitening time-domain coherent line subtraction.

    Identifies or accepts prominent narrow spectral lines (such as 60 Hz mains harmonics
    and 500 Hz violin resonances) and subtracts them coherently in the time domain using
    narrowband heterodyne complex demodulation, eliminating line ringing and residual
    ACF correlations without distorting broadband gravitational-wave chirps.

    Parameters
    ----------
    timeseries : TimeSeries
        Input strain time series.
    lines : list of (f0, excess, fwhm) or list of float, optional
        Precomputed line frequencies or tuples. If None, detected from psd.
    psd : FrequencySeries, optional
        PSD estimate used for line detection if lines is None.
    f_low : float, optional
        Minimum frequency for line removal in Hz (default: 18.0).
    line_threshold : float, optional
        Line prominence threshold in ratio units (default: 3.5).
    max_lines : int, optional
        Maximum number of narrow lines to subtract (default: 64).
    bandwidth : float, optional
        Demodulation filter bandwidth in Hz (default: 0.2).

    Returns
    -------
    ts_clean : TimeSeries
        Strain time series with coherent spectral lines subtracted.
    subtracted_lines : list of tuple
        Catalog of subtracted line parameters (f0, peak_excess, fwhm).
    """
    from scipy.ndimage import median_filter
    fs = float(timeseries.sample_rate)
    N = len(timeseries)
    dt = float(timeseries.delta_t)
    df = 1.0 / (N * dt)
    freqs = np.fft.rfftfreq(N, dt)

    if lines is None:
        if psd is None:
            seg_len = min(int(round(8.0 * fs)), max(int(round(2.0 * fs)), N // 4))
            seg_stride = seg_len // 2
            psd = pycbc.psd.welch(timeseries, seg_len=seg_len, seg_stride=seg_stride, avg_method='median')

        p_vals = psd.numpy()
        f_pilot = psd.sample_frequencies.numpy()
        df_pilot = float(psd.delta_f)

        k_med = max(3, int(round(4.0 / df_pilot)))
        if k_med % 2 == 0:
            k_med += 1
        p_med = median_filter(p_vals, size=k_med)
        ratio = p_vals / np.maximum(p_med, 1e-60)

        f_search_low = max(25.0, f_low)
        in_range = np.where((f_pilot >= f_search_low) & (f_pilot <= (fs / 2.0) - 5.0) & (ratio >= line_threshold))[0]
        detected = []
        if len(in_range) > 0:
            clusters = np.split(in_range, np.where(np.diff(in_range) > 1)[0] + 1)
            for c in clusters:
                cluster_width = float(len(c) * df_pilot)
                # True instrumental lines are narrow. Broad features are slopes/glitches.
                if cluster_width > 2.0:
                    continue
                peak_idx = c[np.argmax(ratio[c])]
                f0 = float(f_pilot[peak_idx])
                excess = float(ratio[peak_idx])
                fwhm = max(0.05, min(1.0, cluster_width))
                detected.append((f0, excess, fwhm))
        detected.sort(key=lambda item: item[1], reverse=True)
        lines = detected[:max_lines]

    line_list = []
    for item in lines:
        if isinstance(item, (int, float, np.number)):
            line_list.append((float(item), 4.0, bandwidth))
        elif isinstance(item, (list, tuple)) and len(item) >= 3:
            line_list.append((float(item[0]), float(item[1]), float(item[2])))
        elif isinstance(item, (list, tuple)) and len(item) == 1:
            line_list.append((float(item[0]), 4.0, bandwidth))

    X = np.fft.rfft(timeseries.numpy())
    T = np.ones(len(freqs), dtype=np.float64)

    min_sub_f = max(25.0, f_low)
    for f0, excess, fwhm in line_list:
        if f0 < min_sub_f or f0 >= fs / 2.0:
            continue
        sigma_f = max(0.01, min(0.20, fwhm / 2.355 if fwhm > 0 else bandwidth / 2.355))
        atten = 1.0 - 1.0 / np.sqrt(max(1.01, excess))
        idx_min = max(0, int((f0 - 4.0 * sigma_f) / df))
        idx_max = min(len(freqs), int((f0 + 4.0 * sigma_f) / df) + 1)
        if idx_min < idx_max:
            f_sub = freqs[idx_min:idx_max]
            gauss = np.exp(-0.5 * ((f_sub - f0) / sigma_f)**2)
            T[idx_min:idx_max] *= (1.0 - atten * gauss)

    x_clean = np.fft.irfft(X * T, n=N)
    ts_clean = TimeSeries(x_clean, delta_t=dt, epoch=timeseries.start_time)
    return ts_clean, line_list


def overwhiten_strain(timeseries, psd=None, f_low=18.0, f_taper=4.0, f_high=None, max_filter_duration=4.0,
                      use_overlap_save=True, chunk_len_s=16.0, seg_len=None, seg_stride=None,
                      avg_method='median', model_compact=None, track_nonstationarity=False,
                      subtract_lines=False, line_threshold=3.5, max_lines=64,
                      inpaint_edges=False, edge_pad_duration=None):
    """Convenience driver to overwhiten strain data in a single call.

    Parameters
    ----------
    timeseries : TimeSeries
        Input PyCBC TimeSeries.
    psd : FrequencySeries, optional
        Precomputed base PSD. If None, estimated using Welch's method.
    f_low : float, optional
        Low frequency cutoff in Hz (default: 18.0).
    f_taper : float, optional
        Taper width below f_low in Hz (default: 4.0).
    f_high : float, optional
        High frequency cutoff in Hz (default: None).
    max_filter_duration : float, optional
        FIR filter length in seconds (default: 4.0).
    use_overlap_save : bool, optional
        Whether to use continuous overlap-save convolution (default: True).
    chunk_len_s : float, optional
        Chunk size in seconds (default: 16.0).
    seg_len : int, optional
        Welch segment length in samples if PSD is estimated (default: 4 * fs).
    seg_stride : int, optional
        Welch segment stride in samples if PSD is estimated (default: 2 * fs).
    avg_method : str, optional
        Welch average method: 'median', 'mean', or 'median-mean' (default: 'median').
    model_compact : CompactInstantaneousPSD, optional
        Precomputed compact instantaneous PSD model with tracked scale factors.
    track_nonstationarity : bool, optional
        Whether to fit compact instantaneous PSD and track 1 Hz multi-band scale factors (default: False).
    subtract_lines : bool, optional
        Whether to subtract coherent sinusoidal lines (default: False).
    line_threshold : float, optional
        Threshold for line detection (default: 3.5).
    max_lines : int, optional
        Maximum number of lines to subtract (default: 64).
    inpaint_edges : bool, optional
        If True, pad the contiguous block with zero-padding on both ends and inpaint
        the outer margins using regularized hole-filling before filtering, eliminating
        boundary step transients and impulse response ringing (default: False).
    edge_pad_duration : float, optional
        Duration of edge inpainting padding in seconds. If None, defaults to min(max_filter_duration, 2.0).

    Returns
    -------
    ts_dow : TimeSeries
        Overwhitened strain.
    ts_dw : TimeSeries
        Intermediate whitened strain.
    metadata : dict
        Dictionary containing psd_base, fir_dow, fir_w, and kernel arrays.
    """
    fs = float(timeseries.sample_rate)
    duration = float(timeseries.duration)

    subtracted_lines = []
    if subtract_lines:
        timeseries, subtracted_lines = subtract_coherent_lines(
            timeseries, psd=psd, f_low=f_low, line_threshold=line_threshold, max_lines=max_lines
        )

    if psd is None:
        if seg_len is None:
            max_seg = int(round(max_filter_duration * fs))
            seg_len = max(max_seg, int(round(4.0 * fs)))
            if seg_len > len(timeseries) // 2:
                seg_len = max(int(round(4.0 * fs)), len(timeseries) // 4)
        if seg_stride is None:
            seg_stride = seg_len // 2
        psd = pycbc.psd.welch(timeseries, seg_len=seg_len, seg_stride=seg_stride, avg_method=avg_method)

    if track_nonstationarity and model_compact is None:
        from pycbc.psd.model import CompactInstantaneousPSD
        model_compact = CompactInstantaneousPSD(f_min=max(15.0, f_low - f_taper), f_max=min(1024.0, fs / 2.0), n_knots=48)
        model_compact.fit_static_from_psd(psd)
        model_compact.track_nonstationarity(timeseries, t_step=1.0, t_window=4.0)

    invpsd_full, w_kernel_full, fir_dow, fir_w, psd_ist_ow, delta_f_full = construct_regularized_kernels(
        psd, fs, duration, f_low=f_low, f_taper=f_taper, f_high=f_high, max_filter_duration=max_filter_duration
    )

    if model_compact is not None and model_compact.band_scales is not None:
        bands = model_compact.bands
        n_bands = len(bands)
        N_full = len(timeseries)
        freqs_full = np.fft.rfftfreq(N_full, 1.0 / fs)
        n_freqs = len(freqs_full)

        # Multi-band partition of unity phi_b(f) spanning full Nyquist bandwidth
        phi = np.zeros((n_bands, n_freqs), dtype=np.float64)
        for b, (flo, fhi, _) in enumerate(bands):
            f_hi_eff = fhi if b < n_bands - 1 else max(fhi, freqs_full[-1] + 1.0)
            f_lo_eff = flo if b > 0 else min(flo, 0.0)
            mask = (freqs_full >= f_lo_eff) & (freqs_full < f_hi_eff)
            phi[b, mask] = 1.0
        tot_phi = np.sum(phi, axis=0)
        tot_phi[tot_phi == 0] = 1.0
        for b in range(n_bands):
            phi[b] /= tot_phi

        # Sub-band FIR filters for whitening and overwhitening
        N_filt = int(round(max_filter_duration * fs))
        k_mid = N_full // 2
        fir_w_bands = []
        fir_dow_bands = []
        for b in range(n_bands):
            w_band = w_kernel_full * phi[b]
            dow_band = invpsd_full * phi[b]
            h_w_b = np.roll(np.fft.irfft(w_band, n=N_full), N_full // 2)
            h_dow_b = np.roll(np.fft.irfft(dow_band, n=N_full), N_full // 2)
            fir_w_bands.append(h_w_b[k_mid - N_filt // 2 : k_mid + N_filt // 2 + 1])
            fir_dow_bands.append(h_dow_b[k_mid - N_filt // 2 : k_mid + N_filt // 2 + 1])

        # Sample-wise continuous scale factors g_b(t)
        t_samples = timeseries.sample_times.numpy()
        g_samples = np.zeros((n_bands, N_full), dtype=np.float64)
        for b in range(n_bands):
            g_samples[b] = np.interp(t_samples, model_compact.scale_times, model_compact.band_scales[b])
            g_samples[b] = np.maximum(g_samples[b], 1e-4)

        x = timeseries.numpy()
        if use_overlap_save:
            from scipy.fft import next_fast_len
            L = int(round(chunk_len_s * fs))
            M = len(fir_w_bands[0])
            k_center = (M - 1) // 2
            n_fft = next_fast_len(L + 2 * M)

            fir_w_bands_fd = [np.fft.rfft(fir, n=n_fft) for fir in fir_w_bands]
            fir_dow_bands_fd = [np.fft.rfft(fir, n=n_fft) for fir in fir_dow_bands]

            y_w = np.zeros(N_full, dtype=np.float64)
            y_dow = np.zeros(N_full, dtype=np.float64)

            for s in range(0, N_full, L):
                e = min(s + L, N_full)
                block_len = e - s
                cs = s - k_center
                ce = e + (M - 1 - k_center)

                pad_left = max(0, -cs)
                pad_right = max(0, ce - N_full)
                sub = x[max(0, cs) : min(N_full, ce)]
                if pad_left > 0 or pad_right > 0:
                    pad_mode = 'reflect' if len(sub) > max(pad_left, pad_right) else 'edge'
                    sub = np.pad(sub, (pad_left, pad_right), mode=pad_mode)

                sub_fd = np.fft.rfft(sub, n=n_fft)
                chunk_w = np.zeros(block_len, dtype=np.float64)
                chunk_dow = np.zeros(block_len, dtype=np.float64)
                for b in range(n_bands):
                    c_w_b = np.fft.irfft(sub_fd * fir_w_bands_fd[b], n=n_fft)[M - 1 : M - 1 + block_len]
                    c_dow_b = np.fft.irfft(sub_fd * fir_dow_bands_fd[b], n=n_fft)[M - 1 : M - 1 + block_len]
                    chunk_w += c_w_b / g_samples[b, s:e]
                    chunk_dow += c_dow_b / g_samples[b, s:e]

                y_w[s:e] = chunk_w
                y_dow[s:e] = chunk_dow

            ts_dw = TimeSeries(y_w, delta_t=1.0 / fs, epoch=timeseries.start_time)
            ts_dow = TimeSeries(y_dow, delta_t=1.0 / fs, epoch=timeseries.start_time)
        else:
            y_w = np.zeros(N_full, dtype=np.float64)
            y_dow = np.zeros(N_full, dtype=np.float64)
            for b in range(n_bands):
                c_w_b = sig.fftconvolve(x, fir_w_bands[b], mode='same')
                c_dow_b = sig.fftconvolve(x, fir_dow_bands[b], mode='same')
                y_w += c_w_b / g_samples[b]
                y_dow += c_dow_b / g_samples[b]
            ts_dw = TimeSeries(y_w, delta_t=1.0 / fs, epoch=timeseries.start_time)
            ts_dow = TimeSeries(y_dow, delta_t=1.0 / fs, epoch=timeseries.start_time)
    else:
        if inpaint_edges:
            from pycbc.strain.gate import gate_and_paint
            pad_dur = min(max_filter_duration, 2.0) if edge_pad_duration is None else float(edge_pad_duration)
            pad_n = int(round(pad_dur * fs))
            N_orig = len(timeseries)

            ts_ext = timeseries.copy()
            ts_ext.prepend_zeros(pad_n)
            ts_ext.append_zeros(pad_n)
            N_tot = len(ts_ext)
            dur_ext = float(ts_ext.duration)

            psd_ext = pycbc.psd.interpolate(psd, ts_ext.delta_f)
            invpsd_ext, w_ext, fir_dow_ext, fir_w_ext, psd_ist_ow, df_ext = construct_regularized_kernels(
                psd_ext, fs, dur_ext, f_low=f_low, f_taper=f_taper, max_filter_duration=max_filter_duration
            )
            invpsd_fs = FrequencySeries(invpsd_ext, delta_f=ts_ext.delta_f)
            if hasattr(ts_ext, 'precision') and ts_ext.precision == 'single':
                invpsd_fs = invpsd_fs.astype(np.float32)

            painted = gate_and_paint(ts_ext, 0, pad_n, invpsd_fs, method='toeplitz')
            painted = gate_and_paint(painted, pad_n + N_orig, N_tot, invpsd_fs, method='toeplitz')

            dow_ext, dw_ext = produce_conditioned_streams(
                painted, invpsd_ext, w_ext, fir_dow_ext, fir_w_ext,
                delta_f_full=df_ext, use_overlap_save=False
            )

            ts_dow = dow_ext[pad_n : pad_n + N_orig]
            ts_dw = dw_ext[pad_n : pad_n + N_orig]
            ts_dow.start_time = timeseries.start_time
            ts_dw.start_time = timeseries.start_time

            invpsd_full = invpsd_ext
            w_kernel_full = w_ext
            fir_dow = fir_dow_ext
            fir_w = fir_w_ext
            delta_f_full = df_ext
        else:
            ts_dow, ts_dw = produce_conditioned_streams(
                timeseries, invpsd_full, w_kernel_full, fir_dow, fir_w,
                delta_f_full=delta_f_full, use_overlap_save=use_overlap_save, chunk_len_s=chunk_len_s
            )

    metadata = {
        "psd_base": psd,
        "psd_ist_ow": psd_ist_ow,
        "invpsd_full": invpsd_full,
        "w_kernel_full": w_kernel_full,
        "fir_dow": fir_dow,
        "fir_w": fir_w,
        "delta_f_full": delta_f_full,
        "f_low": f_low,
        "f_taper": f_taper,
        "max_filter_duration": max_filter_duration,
        "chunk_len_s": chunk_len_s,
        "model_compact": model_compact,
        "subtracted_lines": subtracted_lines,
        "inpainted_edges": inpaint_edges,
        "edge_pad_duration": pad_dur if inpaint_edges else 0.0,
    }
    ts_dow._is_overwhitened = True
    return ts_dow, ts_dw, metadata


def overwhiten_block_strain(strain_block, injections=None, psd=None, f_low=18.0, f_taper=4.0,
                           max_filter_duration=4.0, use_overlap_save=True, chunk_len_s=16.0,
                           seg_len=None, seg_stride=None, avg_method='median',
                           model_compact=None, track_nonstationarity=False,
                           subtract_lines=False, line_threshold=3.5, max_lines=64,
                           inpaint_edges=False, edge_pad_duration=None):
    """Condition a raw strain block in-memory for search pipelines (e.g. pycbc_inspiral).

    Takes a continuous raw strain block (typically 256s - 2048s), optionally adds
    software injections into the raw strain data, and performs high-throughput
    regularized overwhitening and whitening in memory (>300x real-time).

    Parameters
    ----------
    strain_block : TimeSeries
        Raw detector strain block.
    injections : TimeSeries, iterable of TimeSeries, or callable, optional
        Software injection(s) to add to the raw strain block prior to conditioning.
        If callable, called with strain_block as `strain_block = injections(strain_block)`.
    psd : FrequencySeries, optional
        Base PSD estimate. If None, estimated using Welch median over the block.
    f_low : float, optional
        Low frequency cutoff in Hz (default: 18.0).
    f_taper : float, optional
        Taper width below f_low in Hz (default: 4.0).
    max_filter_duration : float, optional
        FIR filter length in seconds (default: 4.0).
    use_overlap_save : bool, optional
        Whether to use continuous overlap-save convolution (default: True).
    chunk_len_s : float, optional
        Chunk size in seconds (default: 16.0).
    seg_len : int, optional
        Welch segment length in samples if PSD is estimated.
    seg_stride : int, optional
        Welch segment stride in samples if PSD is estimated.
    avg_method : str, optional
        Welch averaging method (default: 'median').
    model_compact : CompactInstantaneousPSD, optional
        Precomputed compact instantaneous PSD model.
    track_nonstationarity : bool, optional
        Whether to track multi-band non-stationarity (default: False).
    subtract_lines : bool, optional
        Whether to subtract coherent sinusoidal lines (default: False).
    line_threshold : float, optional
        Threshold for line detection (default: 3.5).
    max_lines : int, optional
        Maximum number of lines to subtract (default: 64).

    Returns
    -------
    ts_dow : TimeSeries
        Overwhitened strain block ready for matched filtering.
    ts_dw : TimeSeries
        Whitened strain block (calibrated discrete unit variance).
    metadata : dict
        Whitening metadata dictionary, including block duration and sample rate.
    """
    if not isinstance(strain_block, TimeSeries):
        raise TypeError("strain_block must be a PyCBC TimeSeries")

    raw_block = strain_block.copy()
    if injections is not None:
        if callable(injections):
            raw_block = injections(raw_block)
        elif isinstance(injections, TimeSeries):
            if len(injections) == len(raw_block) and abs(injections.sample_rate - raw_block.sample_rate) < 1e-6:
                raw_block = raw_block + injections
            else:
                t0_inj = float(injections.start_time)
                t0_raw = float(raw_block.start_time)
                dt = float(raw_block.delta_t)
                idx_start = int(round((t0_inj - t0_raw) / dt))
                idx_end = idx_start + len(injections)
                overlap_start = max(0, idx_start)
                overlap_end = min(len(raw_block), idx_end)
                if overlap_start < overlap_end:
                    inj_start = overlap_start - idx_start
                    inj_end = inj_start + (overlap_end - overlap_start)
                    raw_block.data[overlap_start:overlap_end] += injections.data[inj_start:inj_end]
        elif hasattr(injections, '__iter__'):
            for inj in injections:
                if isinstance(inj, TimeSeries):
                    if len(inj) == len(raw_block) and abs(inj.sample_rate - raw_block.sample_rate) < 1e-6:
                        raw_block = raw_block + inj
                    else:
                        t0_inj = float(inj.start_time)
                        t0_raw = float(raw_block.start_time)
                        dt = float(raw_block.delta_t)
                        idx_start = int(round((t0_inj - t0_raw) / dt))
                        idx_end = idx_start + len(inj)
                        overlap_start = max(0, idx_start)
                        overlap_end = min(len(raw_block), idx_end)
                        if overlap_start < overlap_end:
                            inj_start = overlap_start - idx_start
                            inj_end = inj_start + (overlap_end - overlap_start)
                            raw_block.data[overlap_start:overlap_end] += inj.data[inj_start:inj_end]

    ts_dow, ts_dw, metadata = overwhiten_strain(
        raw_block, psd=psd, f_low=f_low, f_taper=f_taper,
        max_filter_duration=max_filter_duration,
        use_overlap_save=use_overlap_save, chunk_len_s=chunk_len_s,
        seg_len=seg_len, seg_stride=seg_stride, avg_method=avg_method,
        model_compact=model_compact, track_nonstationarity=track_nonstationarity,
        subtract_lines=subtract_lines, line_threshold=line_threshold,
        max_lines=max_lines,
        inpaint_edges=inpaint_edges, edge_pad_duration=edge_pad_duration
    )
    metadata["block_duration"] = float(strain_block.duration)
    metadata["sample_rate"] = float(strain_block.sample_rate)
    return ts_dow, ts_dw, metadata


def matched_filter_overwhitened(template, ts_dow, snr_optimal=None, psd=None, f_low=18.0, f_high=None):
    """Direct matched filter against pre-overwhitened strain d_ow(t).

    Calculates:
        rho(t) = (4 * delta_f / snr_optimal) * IFFT[ h_tilde*(f) * d_ow_tilde(f) ]

    Since d_ow has already been filtered by K_reg(f) = 1 / S_n(f), the matched filter
    is a pure static template correlation without needing per-segment PSD re-weighting.

    Parameters
    ----------
    template : FrequencySeries or TimeSeries
        Gravitational wave template waveform.
    ts_dow : TimeSeries or FrequencySeries
        Overwhitened continuous strain time series or frequency series d_ow.
    snr_optimal : float, optional
        Optimal template SNR (sqrt(sigmasq)) for normalization.
        If None, calculated from template and psd.
    psd : FrequencySeries, optional
        Base PSD used to calculate snr_optimal if snr_optimal is None.
    f_low : float, optional
        Low frequency cutoff in Hz (default: 18.0).
    f_high : float, optional
        High frequency cutoff in Hz (default: Nyquist).

    Returns
    -------
    snr : TimeSeries
        Complex matched filter SNR time series rho(t).
    """
    from pycbc.filter.matchedfilter import make_frequency_series, get_cutoff_indices

    if isinstance(ts_dow, FrequencySeries):
        N = (len(ts_dow) - 1) * 2
        delta_f = ts_dow.delta_f
        fs = float(N * delta_f)
        dow_fd = ts_dow
        epoch = ts_dow.epoch
    elif isinstance(ts_dow, TimeSeries):
        N = len(ts_dow)
        fs = float(ts_dow.sample_rate)
        delta_f = ts_dow.delta_f
        dow_fd = ts_dow.to_frequencyseries()
        epoch = ts_dow.start_time
    else:
        raise TypeError("ts_dow must be a PyCBC TimeSeries or FrequencySeries")

    htilde = make_frequency_series(template)
    if len(htilde) != N // 2 + 1:
        htilde = htilde.copy()
        htilde.resize(N // 2 + 1)

    if snr_optimal is None:
        if psd is not None:
            if abs(psd.delta_f - htilde.delta_f) > 1e-9:
                psd = pycbc.psd.interpolate(psd, htilde.delta_f)
            sigmasq_val = pycbc.filter.sigmasq(htilde, psd=psd, low_frequency_cutoff=f_low, high_frequency_cutoff=f_high)
            snr_optimal = float(np.sqrt(sigmasq_val))
        else:
            raise ValueError("Either snr_optimal or psd must be provided to normalize matched_filter_overwhitened.")

    kmin, kmax = get_cutoff_indices(f_low, f_high, delta_f, N)

    qtilde = zeros(N, dtype=np.complex128)
    qtilde[kmin:kmax] = htilde[kmin:kmax].conj() * dow_fd[kmin:kmax]

    _q = zeros(N, dtype=np.complex128)
    pycbc.fft.ifft(qtilde, _q)

    norm = (4.0 * delta_f) / snr_optimal
    snr = TimeSeries(_q, delta_t=1.0 / fs, epoch=epoch) * norm
    return snr


def verify_time_varying_streaming_seams(ts, psd_base, model_compact, chunk_len_s=16.0, filter_duration=4.0):
    """Test and verify continuous overlap-save streaming with time-varying multi-band scale factors g_b(t).

    Parameters
    ----------
    ts : TimeSeries
        Input strain time series.
    psd_base : FrequencySeries
        Base static PSD.
    model_compact : CompactInstantaneousPSD
        Fitted compact instantaneous PSD model with tracked band scales.
    chunk_len_s : float, optional
        Chunk length in seconds (default: 16.0).
    filter_duration : float, optional
        Filter duration in seconds (default: 4.0).

    Returns
    -------
    result : dict
        Dictionary of seam errors and pass/fail status.
    """
    fs = float(ts.sample_rate)
    N_full = len(ts)
    dt = 1.0 / fs
    duration = N_full / fs

    _, w_kernel_full, _, _, _, _ = construct_regularized_kernels(
        psd_base, fs, duration, f_low=18.0, f_taper=4.0, max_filter_duration=filter_duration
    )

    bands = model_compact.bands
    n_bands = len(bands)
    freqs = np.fft.rfftfreq(N_full, dt)
    n_freqs = len(freqs)

    # Multi-band partition of unity phi_b(f)
    phi = np.zeros((n_bands, n_freqs), dtype=np.float64)
    for b, (flo, fhi, _) in enumerate(bands):
        mask = (freqs >= flo) & (freqs < fhi)
        phi[b, mask] = 1.0
    tot_phi = np.sum(phi, axis=0)
    tot_phi[tot_phi == 0] = 1.0
    for b in range(n_bands):
        phi[b] /= tot_phi

    # Sub-band FIR filters
    N_filt = int(round(filter_duration * fs))
    k_mid = N_full // 2
    fir_w_bands = []
    for b in range(n_bands):
        w_band = w_kernel_full * phi[b]
        h_w_b = np.roll(np.fft.irfft(w_band, n=N_full), N_full // 2)
        fir_w_bands.append(h_w_b[k_mid - N_filt // 2 : k_mid + N_filt // 2 + 1])

    # Sample-wise continuous scale factors g_b(t)
    t_samples = ts.sample_times.numpy()
    g_samples = np.zeros((n_bands, N_full), dtype=np.float64)
    for b in range(n_bands):
        g_samples[b] = np.interp(t_samples, model_compact.scale_times, model_compact.band_scales[b])

    # Continuous reference (full unchunked convolution with sample-wise scaling)
    x = ts.numpy()
    y_ref = np.zeros(N_full, dtype=np.float64)
    for b in range(n_bands):
        y_b = sig.fftconvolve(x, fir_w_bands[b], mode='same')
        y_ref += y_b / g_samples[b]

    # Chunked overlap-save streaming
    L = int(round(chunk_len_s * fs))
    M = len(fir_w_bands[0])
    k_center = (M - 1) // 2

    y_stream = np.zeros(N_full, dtype=np.float64)
    for s in range(0, N_full, L):
        e = min(s + L, N_full)
        block_len = e - s
        cs = s - k_center
        ce = e + (M - 1 - k_center)

        pad_left = max(0, -cs)
        pad_right = max(0, ce - N_full)
        sub = x[max(0, cs) : min(N_full, ce)]
        if pad_left > 0 or pad_right > 0:
            sub = np.pad(sub, (pad_left, pad_right), mode='constant')

        y_chunk = np.zeros(block_len, dtype=np.float64)
        for b in range(n_bands):
            c_b = sig.fftconvolve(sub, fir_w_bands[b], mode='valid')[:block_len]
            y_chunk += c_b / g_samples[b, s:e]
        y_stream[s:e] = y_chunk

    # Evaluate seam error across interior chunk boundaries
    edge_skip = int(round((filter_duration + 2.0) * fs))
    valid_slice = slice(edge_skip, N_full - edge_skip)
    std_ref = float(np.std(y_ref[valid_slice]))
    diff = np.abs(y_stream[valid_slice] - y_ref[valid_slice])
    max_rel_err = float(np.max(diff) / std_ref)

    seam_errs = []
    chunk_times_s = range(int(chunk_len_s), int(duration - chunk_len_s), int(chunk_len_s))
    for c_time in chunk_times_s:
        idx = int(round(c_time * fs))
        if idx - 5 >= 0 and idx + 5 <= N_full:
            seam_diff = np.max(np.abs(y_stream[idx-5:idx+5] - y_ref[idx-5:idx+5]))
            seam_errs.append(float(seam_diff / std_ref))

    max_seam_err = float(np.max(seam_errs)) if seam_errs else 0.0
    mean_seam_err = float(np.mean(seam_errs)) if seam_errs else 0.0

    return {
        "max_interior_relative_error": max_rel_err,
        "max_seam_relative_error": max_seam_err,
        "mean_seam_relative_error": mean_seam_err,
        "chunk_len_s": chunk_len_s,
        "passed": bool(max_seam_err < 1e-4)
    }


class OverwhitenedSlice:
    """Lightweight queried slice of conditioned overwhitened data.

    Attributes
    ----------
    dow : TimeSeries
        Overwhitened strain stream d_ow(t) for direct matched filtering.
    dw : TimeSeries
        Companion whitened strain stream d_w(t) for statistical verification.
    times : ndarray
        Sample timestamps in GPS seconds.
    is_valid : ndarray of bool
        Boolean mask where True denotes certified valid data (past boundary padding and not gated).
    is_gated : ndarray of bool
        Boolean mask where True denotes gated / inpainted data.
    is_boundary : ndarray of bool
        Boolean mask where True denotes boundary transient padding.
    """
    def __init__(self, dow, dw, times, is_valid, is_gated, is_boundary, psd_func=None):
        self.dow = dow
        self.dw = dw
        self.times = times
        self.is_valid = is_valid
        self.is_gated = is_gated
        self.is_boundary = is_boundary
        self._psd_func = psd_func

    @property
    def start_time(self):
        return float(self.dow.start_time)

    @property
    def end_time(self):
        return float(self.dow.end_time)

    @property
    def duration(self):
        return float(self.dow.duration)

    @property
    def sample_rate(self):
        return float(self.dow.sample_rate)

    @property
    def delta_t(self):
        return float(self.dow.delta_t)

    def get_psd(self, delta_f=None, f_low=None, f_high=None, length=None):
        """Evaluate the instantaneous PSD at the midpoint of this slice."""
        if self._psd_func is not None:
            t_mid = (self.start_time + self.end_time) / 2.0
            return self._psd_func(t_mid, delta_f=delta_f, f_low=f_low, f_high=f_high, length=length)
        raise ValueError("Instantaneous PSD function not available for this slice.")

    def __repr__(self):
        valid_pct = np.mean(self.is_valid) * 100 if len(self.is_valid) > 0 else 0.0
        gated_pct = np.mean(self.is_gated) * 100 if len(self.is_gated) > 0 else 0.0
        return (f"<OverwhitenedSlice GPS [{self.start_time:.2f}, {self.end_time:.2f}], "
                f"dur={self.duration:.2f}s, valid={valid_pct:.1f}%, gated={gated_pct:.1f}%>")


class _ChunkInfo:
    """Internal metadata descriptor for a single HDF5 chunk file."""
    def __init__(self, file_path, start_time, duration, sample_rate, delta_t,
                 boundary_pad, gated_segments, valid_segments, has_psd_model, has_psd_base):
        self.file_path = file_path
        self.start_time = float(start_time)
        self.duration = float(duration)
        self.end_time = self.start_time + self.duration
        self.sample_rate = float(sample_rate)
        self.delta_t = float(delta_t)
        self.boundary_pad = float(boundary_pad)
        self.gated_segments = gated_segments
        self.valid_segments = valid_segments
        self.has_psd_model = bool(has_psd_model)
        self.has_psd_base = bool(has_psd_base)
        self._cached_psd_model = None
        self._cached_psd_base = None

    def get_psd_model(self):
        if not self.has_psd_model:
            return None
        if self._cached_psd_model is None:
            from pycbc.psd.model import CompactInstantaneousPSD
            with h5py.File(self.file_path, "r") as f:
                if "instantaneous_psd" in f:
                    self._cached_psd_model = CompactInstantaneousPSD.load_from_hdf(f["instantaneous_psd"])
        return self._cached_psd_model

    def get_psd_base(self):
        if not self.has_psd_base:
            return None
        if self._cached_psd_base is None:
            with h5py.File(self.file_path, "r") as f:
                if "psd_base" in f:
                    vals = f["psd_base/values"][:]
                    df = float(f["psd_base"].attrs.get("delta_f", 1.0))
                    self._cached_psd_base = FrequencySeries(vals, delta_f=df)
        return self._cached_psd_base


class OverwhitenedData:
    """High-level, user-friendly interface to query continuous overwhitened data,
    companion whitened data, instantaneous PSD models, and data validity/gating status.

    Parameters
    ----------
    source : str, Path, list of str/Path, or h5py.File
        File path to an overwhitened HDF5 product, directory containing HDF5 chunk files,
        or list of file paths.
    boundary_pad : float, optional
        Filter transient boundary padding in seconds (overrides file metadata if provided).

    Examples
    --------
    >>> data = OverwhitenedData.open("conditioned_strain.hdf5")
    >>> print(data.start_time, data.end_time, data.valid_segments)
    >>> dow = data.get_overwhitened(start=1369375100, end=1369375200)
    >>> psd = data.get_psd(t=1369375150.0, delta_f=0.125)
    >>> is_clean = data.is_valid(1369375150.0)
    """
    def __init__(self, source, boundary_pad=None):
        self.chunks = []
        self._open_handles = {}

        if isinstance(source, (str, Path)):
            src_path = Path(source)
            if src_path.is_dir():
                file_list = sorted(list(src_path.glob("*.hdf5")) + list(src_path.glob("*.h5")) + list(src_path.glob("*.hdf")))
                if not file_list:
                    raise FileNotFoundError(f"No HDF5 files found in directory {source}")
                self._index_files(file_list, boundary_pad)
            elif src_path.is_file():
                self._index_files([src_path], boundary_pad)
            else:
                glob_files = sorted(glob.glob(str(source)))
                if glob_files:
                    self._index_files(glob_files, boundary_pad)
                else:
                    raise FileNotFoundError(f"Source file or directory not found: {source}")
        elif isinstance(source, (list, tuple)):
            self._index_files(source, boundary_pad)
        else:
            raise TypeError(f"Unsupported source type: {type(source)}")

        if not self.chunks:
            raise ValueError("No valid overwhitened data chunks could be indexed.")

        self.chunks.sort(key=lambda c: c.start_time)
        self.start_time = self.chunks[0].start_time
        self.end_time = self.chunks[-1].end_time
        self.duration = self.end_time - self.start_time
        self.sample_rate = self.chunks[0].sample_rate
        self.delta_t = self.chunks[0].delta_t
        self.boundary_pad = boundary_pad if boundary_pad is not None else self.chunks[0].boundary_pad

        # Aggregate segment lists
        self.all_segments = segments.segmentlist([
            segments.segment(c.start_time, c.end_time) for c in self.chunks
        ]).coalesce()

        all_gated = []
        for c in self.chunks:
            all_gated.extend(c.gated_segments)
        self.gated_segments = segments.segmentlist(all_gated).coalesce()

        all_boundary = []
        for c in self.chunks:
            bpad = c.boundary_pad
            all_boundary.append(segments.segment(c.start_time, min(c.end_time, c.start_time + bpad)))
            all_boundary.append(segments.segment(max(c.start_time, c.end_time - bpad), c.end_time))
        self.boundary_segments = segments.segmentlist(all_boundary).coalesce()

        all_valid = []
        for c in self.chunks:
            all_valid.extend(c.valid_segments)
        self.valid_segments = segments.segmentlist(all_valid).coalesce()

    def _index_files(self, file_paths, user_pad=None):
        """Index metadata from an array of HDF5 file paths."""
        for fp in file_paths:
            fp_str = str(fp)
            with h5py.File(fp_str, "r") as f:
                if "dow" not in f:
                    continue
                dset = f["dow"]
                st = float(dset.attrs["start_time"])
                dur = float(dset.attrs["duration"])
                fs = float(dset.attrs.get("sample_rate", 2048.0))
                dt = float(dset.attrs.get("delta_t", 1.0 / fs))

                if user_pad is not None:
                    bpad = float(user_pad)
                elif "boundary_pad" in f.attrs:
                    bpad = float(f.attrs["boundary_pad"])
                elif "metadata" in f and "max_filter_duration" in f["metadata"].attrs:
                    bpad = float(f["metadata"].attrs["max_filter_duration"]) / 2.0
                else:
                    bpad = 2.0

                gated_list = []
                if "gated_segments" in f:
                    arr_g = f["gated_segments"][:]
                    for row in arr_g:
                        if len(row) == 2 and row[1] > row[0]:
                            gated_list.append(segments.segment(float(row[0]), float(row[1])))
                seglist_gated = segments.segmentlist(gated_list).coalesce()

                valid_list = []
                if "valid_segments" in f:
                    arr_v = f["valid_segments"][:]
                    for row in arr_v:
                        if len(row) == 2 and row[1] > row[0]:
                            valid_list.append(segments.segment(float(row[0]), float(row[1])))
                    seglist_valid = (segments.segmentlist(valid_list) - seglist_gated).coalesce()
                else:
                    raw_valid = segments.segmentlist([
                        segments.segment(st + bpad, max(st + bpad, st + dur - bpad))
                    ])
                    seglist_valid = (raw_valid - seglist_gated).coalesce()

                has_psd_model = "instantaneous_psd" in f
                has_psd_base = "psd_base" in f

                chunk = _ChunkInfo(
                    file_path=fp_str, start_time=st, duration=dur,
                    sample_rate=fs, delta_t=dt, boundary_pad=bpad,
                    gated_segments=seglist_gated, valid_segments=seglist_valid,
                    has_psd_model=has_psd_model, has_psd_base=has_psd_base
                )
                self.chunks.append(chunk)

    @classmethod
    def open(cls, source, boundary_pad=None):
        """Classmethod factory opening an OverwhitenedData instance."""
        return cls(source, boundary_pad=boundary_pad)

    def close(self):
        """Close any cached open file handles."""
        for h in self._open_handles.values():
            try:
                h.close()
            except Exception:
                pass
        self._open_handles.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def valid_intervals(self):
        """List of (start_gps, end_gps) tuples for valid data."""
        return [(float(s[0]), float(s[1])) for s in self.valid_segments]

    @property
    def gated_intervals(self):
        """List of (start_gps, end_gps) tuples for gated/inpainted data."""
        return [(float(s[0]), float(s[1])) for s in self.gated_segments]

    def _get_h5_file(self, file_path):
        if file_path not in self._open_handles:
            self._open_handles[file_path] = h5py.File(file_path, "r")
        return self._open_handles[file_path]

    def is_valid(self, t):
        """Check if GPS time `t` falls in certified valid data (past boundary, non-gated)."""
        if np.isscalar(t):
            t_flt = float(t)
            return any(s[0] <= t_flt < s[1] for s in self.valid_segments)
        times = np.asarray(t, dtype=np.float64)
        mask = np.zeros(len(times), dtype=bool)
        for s in self.valid_segments:
            mask |= (times >= s[0]) & (times < s[1])
        return mask

    def is_gated(self, t):
        """Check if GPS time `t` falls within a gated / inpainted segment."""
        if np.isscalar(t):
            t_flt = float(t)
            return any(s[0] <= t_flt < s[1] for s in self.gated_segments)
        times = np.asarray(t, dtype=np.float64)
        mask = np.zeros(len(times), dtype=bool)
        for s in self.gated_segments:
            mask |= (times >= s[0]) & (times < s[1])
        return mask

    def is_boundary(self, t):
        """Check if GPS time `t` falls within filter transient boundary padding."""
        if np.isscalar(t):
            t_flt = float(t)
            return any(s[0] <= t_flt < s[1] for s in self.boundary_segments)
        times = np.asarray(t, dtype=np.float64)
        mask = np.zeros(len(times), dtype=bool)
        for s in self.boundary_segments:
            mask |= (times >= s[0]) & (times < s[1])
        return mask

    def get_valid_mask(self, start_time=None, end_time=None):
        """Compute sample-wise boolean validity mask over a given time range."""
        ts_ref = self.get_overwhitened(start_time, end_time)
        return self.is_valid(ts_ref.sample_times.numpy())

    def get_gated_mask(self, start_time=None, end_time=None):
        """Compute sample-wise boolean gating mask over a given time range."""
        ts_ref = self.get_overwhitened(start_time, end_time)
        return self.is_gated(ts_ref.sample_times.numpy())

    def _read_stream_chunk(self, chunk, group_name, start_time, end_time):
        """Read time slice from a single chunk."""
        f = self._get_h5_file(chunk.file_path)
        dset = f[group_name]
        fs = chunk.sample_rate
        dt = chunk.delta_t

        s_idx = max(0, int(round((start_time - chunk.start_time) * fs)))
        e_idx = min(len(dset), int(round((end_time - chunk.start_time) * fs)))
        if s_idx >= e_idx:
            return None

        arr = dset[s_idx:e_idx].astype(np.float64)
        t_epoch = chunk.start_time + s_idx * dt
        return TimeSeries(arr, delta_t=dt, epoch=t_epoch)

    def _query_stream(self, group_name, start_time=None, end_time=None):
        """Query either dow or dw over specified interval."""
        if start_time is None:
            start_time = self.start_time
        if end_time is None:
            end_time = self.end_time

        start_time = float(start_time)
        end_time = float(end_time)
        if end_time <= start_time:
            raise ValueError(f"end_time ({end_time}) must be greater than start_time ({start_time})")

        matching_chunks = [
            c for c in self.chunks
            if not (c.end_time <= start_time or c.start_time >= end_time)
        ]
        if not matching_chunks:
            raise ValueError(f"Requested time range [{start_time}, {end_time}] is outside covered range [{self.start_time}, {self.end_time}]")

        slices = []
        for c in matching_chunks:
            s_ts = self._read_stream_chunk(c, group_name, start_time, end_time)
            if s_ts is not None and len(s_ts) > 0:
                slices.append(s_ts)

        if not slices:
            raise ValueError(f"No samples found in range [{start_time}, {end_time}]")

        if len(slices) == 1:
            return slices[0]

        total_len = sum(len(s) for s in slices)
        combined = np.empty(total_len, dtype=np.float64)
        cur = 0
        for s in slices:
            combined[cur : cur + len(s)] = s.numpy()
            cur += len(s)
        res_ts = TimeSeries(combined, delta_t=slices[0].delta_t, epoch=slices[0].start_time)
        res_ts.gating_info = [(float(s[0]), float(s[1])) for s in self.gated_segments]
        res_ts.valid_segments = self.valid_segments
        res_ts.boundary_pad = self.boundary_pad
        return res_ts

    def get_overwhitened(self, start_time=None, end_time=None):
        """Query overwhitened strain stream d_ow(t) over specified GPS range."""
        return self._query_stream("dow", start_time=start_time, end_time=end_time)

    def get_whitened(self, start_time=None, end_time=None):
        """Query companion whitened strain stream d_w(t) over specified GPS range."""
        return self._query_stream("dw", start_time=start_time, end_time=end_time)

    def get_series(self, which="dow", start_time=None, end_time=None):
        """Retrieve TimeSeries for specified stream ('dow' or 'dw') over range."""
        if which == "dow":
            return self.get_overwhitened(start_time, end_time)
        elif which == "dw":
            return self.get_whitened(start_time, end_time)
        else:
            raise ValueError(f"Unknown stream '{which}'. Choose 'dow' or 'dw'.")

    def time_slice(self, start_time, end_time, which="dow"):
        """PyCBC-idiomatic time_slice interface.

        Parameters
        ----------
        start_time : float
            Start GPS time.
        end_time : float
            End GPS time.
        which : {'dow', 'dw', 'both', 'slice'}, optional
            Stream to return: 'dow' (default), 'dw', or 'both'/'slice' for OverwhitenedSlice.
        """
        if which == "dow":
            return self.get_overwhitened(start_time, end_time)
        elif which == "dw":
            return self.get_whitened(start_time, end_time)
        else:
            return self.query(start_time, end_time)

    def __getitem__(self, key):
        """Python slice indexing syntax: data[start:end] -> OverwhitenedSlice."""
        if isinstance(key, slice):
            return self.query(key.start, key.stop)
        raise TypeError(f"Invalid key type {type(key)}. Expected slice: data[start:end].")

    @property
    def segmentlistdict(self):
        """Return igwn_segments.segmentlistdict of VALID, GATED, and BOUNDARY segments."""
        d = segments.segmentlistdict()
        d["VALID"] = self.valid_segments
        d["GATED"] = self.gated_segments
        d["BOUNDARY"] = self.boundary_segments
        return d

    def get_segment_dict(self, ifo=None):
        """Return segmentlistdict keyed by detector IFO for PyCBC workflow/veto compatibility."""
        det = ifo or getattr(self, "ifo", "STRAIN")
        d = segments.segmentlistdict()
        d[det] = self.valid_segments
        return d

    def get_psd(self, t, length=None, delta_f=None, low_frequency_cutoff=None,
                low_freq_cutoff=None, f_low=None, f_high=None):
        """Evaluate the instantaneous PSD S_n(f, t) at GPS time `t`.

        Follows standard PyCBC PSD function argument conventions (length, delta_f, low_frequency_cutoff).

        Parameters
        ----------
        t : float
            GPS time.
        length : int, optional
            Exact number of frequency points (e.g. for matched filtering).
        delta_f : float, optional
            Frequency resolution in Hz.
        low_frequency_cutoff, low_freq_cutoff, f_low : float, optional
            Frequencies below this cutoff in Hz are set to zero.
        f_high : float, optional
            High frequency cutoff in Hz.

        Returns
        -------
        psd : FrequencySeries
            Instantaneous PSD at time `t`.
        """
        if low_frequency_cutoff is not None:
            f_low = low_frequency_cutoff
        elif low_freq_cutoff is not None:
            f_low = low_freq_cutoff
        t = float(t)
        target_chunk = None
        for c in self.chunks:
            if c.start_time <= t <= c.end_time:
                target_chunk = c
                break
        if target_chunk is None:
            dists = [min(abs(c.start_time - t), abs(c.end_time - t)) for c in self.chunks]
            target_chunk = self.chunks[int(np.argmin(dists))]

        model = target_chunk.get_psd_model()
        fs = target_chunk.sample_rate

        if model is not None:
            if delta_f is None:
                psd_b = target_chunk.get_psd_base()
                delta_f = float(psd_b.delta_f) if psd_b is not None else 0.125
            else:
                delta_f = float(delta_f)

            if length is not None:
                freqs = np.arange(length, dtype=np.float64) * delta_f
            else:
                f_max_eff = float(f_high) if f_high is not None else fs / 2.0
                freqs = np.arange(0.0, f_max_eff + delta_f / 2.0, delta_f, dtype=np.float64)

            vals = model.eval_psd(freqs, t=t, include_lines=True)
            if f_low is not None:
                vals[freqs < float(f_low)] = 0.0

            return FrequencySeries(vals, delta_f=delta_f, epoch=t)

        psd_b = target_chunk.get_psd_base()
        if psd_b is not None:
            if delta_f is not None and abs(float(psd_b.delta_f) - delta_f) > 1e-9:
                return pycbc.psd.interpolate(psd_b, delta_f)
            return psd_b

        raise ValueError(f"No PSD model or base PSD available in chunk for time {t}")

    def query(self, start_time, end_time):
        """Query a continuous time window, returning a unified OverwhitenedSlice."""
        dow = self.get_overwhitened(start_time, end_time)
        dw = self.get_whitened(start_time, end_time)
        times = dow.sample_times.numpy()

        is_valid = self.is_valid(times)
        is_gated = self.is_gated(times)
        is_boundary = self.is_boundary(times)

        return OverwhitenedSlice(
            dow=dow, dw=dw, times=times,
            is_valid=is_valid, is_gated=is_gated, is_boundary=is_boundary,
            psd_func=self.get_psd
        )

    def __repr__(self):
        return (f"<OverwhitenedData GPS [{self.start_time:.2f}, {self.end_time:.2f}], "
                f"dur={self.duration:.2f}s, chunks={len(self.chunks)}, "
                f"valid_segs={len(self.valid_segments)}, gated_segs={len(self.gated_segments)}>")


# Standard aliases
OverwhitenedDataset = OverwhitenedData

