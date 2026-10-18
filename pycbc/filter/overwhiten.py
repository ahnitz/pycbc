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

"""Regularized continuous overwhitening, overlap-save streaming convolution,
and pre-whitening coherent line subtraction.
"""

import numpy as np

import pycbc.psd
from pycbc.types import FrequencySeries, TimeSeries


def construct_regularized_kernels(
    psd_base,
    fs,
    duration,
    f_low=18.0,
    f_taper=4.0,
    f_high=None,
    max_filter_duration=4.0,
    synthesize_full_spectrum=True,
):
    """Construct regularized overwhitening and whitening kernels.

    Constructs:
    - C^2 smooth Hann highpass taper to eliminate low-frequency wall sinc-ringing.
    - Inverse Spectrum Truncation (IST) ensuring compact time-domain impulse
      response <= max_filter_duration.
    - Precomputed centered zero-phase symmetric FIR filters for overlap-save
      continuous streaming.

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
        High frequency cutoff in Hertz. Defaults to 0.9 * (fs / 2).
    max_filter_duration : float, optional
        Maximum duration in seconds of the time-domain filter impulse response
        (default: 4.0).
    synthesize_full_spectrum : bool, optional
        Whether to synthesize full-length frequency-domain arrays when compact
        grid optimization is active (default: True).

    Returns
    -------
    invpsd_full : ndarray or None
        Frequency-domain regularized overwhitening kernel values (None if
        synthesize_full_spectrum=False).
    w_kernel_full : ndarray or None
        Frequency-domain dimensionless whitening kernel values (None if
        synthesize_full_spectrum=False).
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

    # For large duration streams, compute inverse spectrum truncation and FIR taps
    # on a compact frequency grid (e.g. 64s). The FIR filter impulse response is bounded
    # by max_filter_duration (typically <= 4.0s), so computing IST across millions of
    # frequency bins is mathematically redundant.
    dur_compact = min(float(duration), max(64.0, 4.0 * float(max_filter_duration)))
    if dur_compact < duration:
        invpsd_c, w_c, fir_dow, fir_w, psd_ist_ow_c, _ = construct_regularized_kernels(
            psd_base,
            fs,
            dur_compact,
            f_low=f_low,
            f_taper=f_taper,
            f_high=f_high,
            max_filter_duration=max_filter_duration,
            synthesize_full_spectrum=True,
        )
        if not synthesize_full_spectrum:
            return None, None, fir_dow, fir_w, psd_ist_ow_c, delta_f_full

        half_filt = (len(fir_dow) - 1) // 2

        # Synthesize full-grid frequency kernels directly from the compact zero-phase FIR taps
        h_pad = np.zeros(N_full, dtype=np.float64)
        h_pad[: half_filt + 1] = fir_dow[half_filt:]
        h_pad[-half_filt:] = fir_dow[:half_filt]
        invpsd_full = np.fft.rfft(h_pad).real

        h_w_pad = np.zeros(N_full, dtype=np.float64)
        h_w_pad[: half_filt + 1] = fir_w[half_filt:]
        h_w_pad[-half_filt:] = fir_w[:half_filt]
        w_kernel_full = np.fft.rfft(h_w_pad).real

        psd_ist_ow = FrequencySeries(
            1.0 / np.maximum(invpsd_full, 1e-60), delta_f=delta_f_full
        )
        return invpsd_full, w_kernel_full, fir_dow, fir_w, psd_ist_ow, delta_f_full

    psd_full = pycbc.psd.interpolate(psd_base, delta_f_full)
    freqs_full = psd_full.sample_frequencies.numpy()

    # C^2 smooth highpass taper
    taper_full = np.ones(len(freqs_full), dtype=np.float64)
    f_start = max(0.0, f_low - f_taper)
    taper_full[freqs_full < f_start] = 0.0
    mid = (freqs_full >= f_start) & (freqs_full < f_low)
    if np.any(mid) and f_taper > 0:
        taper_full[mid] = 0.5 * (
            1.0 - np.cos(np.pi * (freqs_full[mid] - f_start) / f_taper)
        )
    if f_high is not None:
        taper_full[freqs_full >= f_high] = 0.0

    # Safe filter truncation bounds: max_filter_duration cannot exceed half stream duration
    max_filt_dur_bound = max(
        1.0 / fs, min(float(max_filter_duration), (duration / 2.0) - 1.0 / fs)
    )
    N_filt = min(int(round(max_filt_dur_bound * fs)), N_full - 2)
    if N_filt % 2 != 0:
        N_filt += 1

    # Truncate inverse spectrum
    psd_ist_ow = pycbc.psd.inverse_spectrum_truncation(
        psd_full,
        max_filter_len=N_filt,
        which_spectrum="invpsd",
        low_frequency_cutoff=f_low,
        high_frequency_cutoff=f_high,
        trunc_method="hann",
    )
    psd_ist_w = pycbc.psd.inverse_spectrum_truncation(
        psd_full,
        max_filter_len=N_filt,
        which_spectrum="invasd",
        low_frequency_cutoff=f_low,
        high_frequency_cutoff=f_high,
        trunc_method="hann",
    )

    # Overwhitening kernel with zero-division protection
    psd_ow_safe = np.maximum(psd_ist_ow.numpy(), 1e-60)
    invpsd_full = np.where(taper_full > 0, taper_full / psd_ow_safe, 0.0)

    # Dimensionless Whitening kernel for unit discrete variance:
    # Bandwidth delta_F = exact discrete noise-equivalent bandwidth of the regularized taper
    delta_F = (
        float(
            np.sum(taper_full[1:-1] ** 2)
            + 0.5 * taper_full[0] ** 2
            + 0.5 * taper_full[-1] ** 2
        )
        * delta_f_full
    )
    if delta_F <= 0:
        raise ValueError(f"Sampling rate {fs} Hz is too low for f_low={f_low} Hz.")
    psd_w_safe = np.maximum(psd_ist_w.numpy(), 1e-60)
    w_kernel_full = np.where(
        taper_full > 0, taper_full / np.sqrt(delta_F * psd_w_safe), 0.0
    )

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
    from scipy.fft import next_fast_len

    fs = float(timeseries.sample_rate)
    N = len(timeseries)
    dt = 1.0 / fs
    x = timeseries.numpy()

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
            pad_mode = "reflect" if len(sub) > max(pad_left, pad_right) else "edge"
            sub = np.pad(sub, (pad_left, pad_right), mode=pad_mode)

        sub_fd = np.fft.rfft(sub, n=n_fft)
        c = np.fft.irfft(sub_fd * fir_fd, n=n_fft)[M - 1 : M - 1 + block_len]
        out[s:e] = c

    return TimeSeries(out, delta_t=dt, epoch=timeseries.start_time)


def produce_conditioned_streams(
    timeseries,
    invpsd_full,
    w_kernel_full,
    fir_dow,
    fir_w,
    delta_f_full=None,
    use_overlap_save=True,
    chunk_len_s=16.0,
):
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
        dow_fd = FrequencySeries(
            fd.numpy() * invpsd_full,
            delta_f=delta_f_full,
            epoch=timeseries.start_time,
        )
        dw_fd = FrequencySeries(
            fd.numpy() * w_kernel_full,
            delta_f=delta_f_full,
            epoch=timeseries.start_time,
        )
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
            pad_mode = "reflect" if len(sub) > max(pad_left, pad_right) else "edge"
            sub = np.pad(sub, (pad_left, pad_right), mode=pad_mode)

        sub_fd = np.fft.rfft(sub, n=n_fft)
        c_dow = np.fft.irfft(sub_fd * fir_dow_fd, n=n_fft)[M - 1 : M - 1 + block_len]
        c_w = np.fft.irfft(sub_fd * fir_w_fd, n=n_fft)[M - 1 : M - 1 + block_len]

        d_ow[s:e] = c_dow
        d_w[s:e] = c_w

    ts_dow = TimeSeries(d_ow, delta_t=dt, epoch=timeseries.start_time)
    ts_dw = TimeSeries(d_w, delta_t=dt, epoch=timeseries.start_time)
    return ts_dow, ts_dw


def subtract_coherent_lines(
    timeseries,
    lines=None,
    psd=None,
    f_low=18.0,
    line_threshold=3.5,
    max_lines=64,
    bandwidth=0.2,
):
    """Two-stage pre-whitening time-domain coherent line subtraction.

    Identifies or accepts prominent narrow spectral lines (such as 60 Hz mains harmonics
    and 500 Hz violin resonances) and subtracts them coherently in the time domain,
    eliminating line ringing and residual ACF correlations without distorting
    broadband gravitational-wave chirps.

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
            psd = pycbc.psd.welch(
                timeseries,
                seg_len=seg_len,
                seg_stride=seg_stride,
                avg_method="median",
            )

        p_vals = psd.numpy()
        f_pilot = psd.sample_frequencies.numpy()
        df_pilot = float(psd.delta_f)

        k_med = max(3, int(round(4.0 / df_pilot)))
        if k_med % 2 == 0:
            k_med += 1
        p_med = median_filter(p_vals, size=k_med)
        ratio = p_vals / np.maximum(p_med, 1e-60)

        f_search_low = max(25.0, f_low)
        in_range = np.where(
            (f_pilot >= f_search_low)
            & (f_pilot <= (fs / 2.0) - 5.0)
            & (ratio >= line_threshold)
        )[0]
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
            gauss = np.exp(-0.5 * ((f_sub - f0) / sigma_f) ** 2)
            T[idx_min:idx_max] *= 1.0 - atten * gauss

    x_clean = np.fft.irfft(X * T, n=N)
    ts_clean = TimeSeries(x_clean, delta_t=dt, epoch=timeseries.start_time)
    return ts_clean, line_list


def overwhiten_strain(
    timeseries,
    psd=None,
    f_low=18.0,
    f_taper=4.0,
    f_high=None,
    max_filter_duration=4.0,
    use_overlap_save=True,
    chunk_len_s=16.0,
    seg_len=None,
    seg_stride=None,
    avg_method="median",
    subtract_lines=False,
    line_threshold=3.5,
    max_lines=64,
    inpaint_edges=False,
    edge_pad_duration=None,
    psd_model_type="welch",
):
    """Convenience driver to overwhiten strain data in a single call.

    Parameters
    ----------
    timeseries : TimeSeries
        Input PyCBC TimeSeries.
    psd : FrequencySeries, optional
        Precomputed base PSD. If None, estimated using the requested PSD method.
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
        Duration of edge inpainting padding in seconds. If None, defaults to
        min(max_filter_duration, 2.0).
    psd_model_type : str, optional
        PSD estimation method when `psd` is None: 'welch', 'trimmed-welch',
        'multitaper', or 'line-continuum' (default: 'welch').

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
    pad_dur = (
        min(max_filter_duration, 2.0)
        if edge_pad_duration is None
        else float(edge_pad_duration)
    )

    subtracted_lines = []
    if subtract_lines:
        timeseries, subtracted_lines = subtract_coherent_lines(
            timeseries,
            psd=psd,
            f_low=f_low,
            line_threshold=line_threshold,
            max_lines=max_lines,
        )

    if psd is None:
        if psd_model_type in ("line-continuum", "line_continuum"):
            psd = pycbc.psd.line_continuum(
                timeseries,
                f_low=f_low,
                f_high=f_high,
                seg_len=seg_len,
                seg_stride=seg_stride,
            )
        elif psd_model_type == "multitaper":
            psd = pycbc.psd.multitaper(
                timeseries, seg_len=seg_len, seg_stride=seg_stride
            )
        elif psd_model_type in ("trimmed-welch", "trimmed_welch"):
            psd = pycbc.psd.trimmed_welch(
                timeseries, seg_len=seg_len, seg_stride=seg_stride
            )
        else:
            if seg_len is None:
                max_seg = int(round(max_filter_duration * fs))
                seg_len = max(max_seg, int(round(4.0 * fs)))
                if seg_len > len(timeseries) // 2:
                    seg_len = max(int(round(4.0 * fs)), len(timeseries) // 4)
            if seg_stride is None:
                seg_stride = seg_len // 2
            psd = pycbc.psd.welch(
                timeseries,
                seg_len=seg_len,
                seg_stride=seg_stride,
                avg_method=avg_method,
            )

    synthesize_full = not use_overlap_save
    (
        invpsd_full,
        w_kernel_full,
        fir_dow,
        fir_w,
        psd_ist_ow,
        delta_f_full,
    ) = construct_regularized_kernels(
        psd,
        fs,
        duration,
        f_low=f_low,
        f_taper=f_taper,
        f_high=f_high,
        max_filter_duration=max_filter_duration,
        synthesize_full_spectrum=synthesize_full,
    )

    if inpaint_edges and not use_overlap_save:
        from pycbc.strain.gate import gate_and_paint

        pad_n = int(round(pad_dur * fs))
        N_orig = len(timeseries)

        ts_ext = timeseries.copy()
        ts_ext.prepend_zeros(pad_n)
        ts_ext.append_zeros(pad_n)
        N_tot = len(ts_ext)
        dur_ext = float(ts_ext.duration)

        psd_ext = pycbc.psd.interpolate(psd, ts_ext.delta_f)
        (
            invpsd_ext,
            w_ext,
            fir_dow_ext,
            fir_w_ext,
            psd_ist_ow,
            df_ext,
        ) = construct_regularized_kernels(
            psd_ext,
            fs,
            dur_ext,
            f_low=f_low,
            f_taper=f_taper,
            max_filter_duration=max_filter_duration,
        )
        invpsd_fs = FrequencySeries(invpsd_ext, delta_f=ts_ext.delta_f)
        if hasattr(ts_ext, "precision") and ts_ext.precision == "single":
            invpsd_fs = invpsd_fs.astype(np.float32)

        painted = gate_and_paint(ts_ext, 0, pad_n, invpsd_fs, method="toeplitz")
        painted = gate_and_paint(
            painted, pad_n + N_orig, N_tot, invpsd_fs, method="toeplitz"
        )

        dow_ext, dw_ext = produce_conditioned_streams(
            painted,
            invpsd_ext,
            w_ext,
            fir_dow_ext,
            fir_w_ext,
            delta_f_full=df_ext,
            use_overlap_save=False,
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
            timeseries,
            invpsd_full,
            w_kernel_full,
            fir_dow,
            fir_w,
            delta_f_full=delta_f_full,
            use_overlap_save=use_overlap_save,
            chunk_len_s=chunk_len_s,
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
        "subtracted_lines": subtracted_lines,
        "inpainted_edges": bool(inpaint_edges),
        "edge_pad_duration": pad_dur if inpaint_edges else 0.0,
    }
    return ts_dow, ts_dw, metadata
