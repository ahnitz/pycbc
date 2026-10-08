# Copyright (C) 2012 Tito Dal Canton
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
"""Utilites to estimate PSDs from data.
"""

import numpy
from pycbc.types import Array, FrequencySeries, TimeSeries, zeros
from pycbc.types import real_same_precision_as, complex_same_precision_as
from pycbc.fft import fft, ifft

# Change to True in front-end if you want this function to use caching
# This is a mostly-hidden optimization option that most users will not want
# to use. It is used in PyCBC Live
USE_CACHING_FOR_WELCH_FFTS = False
USE_CACHING_FOR_INV_SPEC_TRUNC = False
# If using caching we want output to be unique if called at different places
# (and if called from different modules/functions), these unique IDs acheive
# that. The numbers are not significant, only that they are unique.
WELCH_UNIQUE_ID = 438716587
INVSPECTRUNC_UNIQUE_ID = 100257896

def median_bias(n):
    """Calculate the bias of the median average PSD computed from `n` segments.

    Parameters
    ----------
    n : int
        Number of segments used in PSD estimation.

    Returns
    -------
    ans : float
        Calculated bias.

    Raises
    ------
    ValueError
        For non-integer or non-positive `n`.

    Notes
    -----
    See arXiv:gr-qc/0509116 appendix B for details.
    """
    if type(n) is not int or n <= 0:
        raise ValueError('n must be a positive integer')
    if n >= 1000:
        return numpy.log(2)
    ans = 1
    for i in range(1, (n - 1) // 2 + 1):
        ans += 1.0 / (2*i + 1) - 1.0 / (2*i)
    return ans

def welch(timeseries, seg_len=4096, seg_stride=2048, window='hann',
          avg_method='median', num_segments=None, require_exact_data_fit=False):
    """PSD estimator based on Welch's method.

    Parameters
    ----------
    timeseries : TimeSeries
        Time series for which the PSD is to be estimated.
    seg_len : int
        Segment length in samples.
    seg_stride : int
        Separation between consecutive segments, in samples.
    window : {'hann', numpy.ndarray}
        Function used to window segments before Fourier transforming, or
        a `numpy.ndarray` that specifies the window.
    avg_method : {'median', 'mean', 'median-mean'}
        Method used for averaging individual segment PSDs.

    Returns
    -------
    psd : FrequencySeries
        Frequency series containing the estimated PSD.

    Raises
    ------
    ValueError
        For invalid choices of `seg_len`, `seg_stride` `window` and
        `avg_method` and for inconsistent combinations of len(`timeseries`),
        `seg_len` and `seg_stride`.

    Notes
    -----
    See arXiv:gr-qc/0509116 for details.
    """
    from pycbc.strain.strain import execute_cached_fft

    window_map = {
        'hann': numpy.hanning
    }

    # sanity checks
    if isinstance(window, numpy.ndarray) and window.size != seg_len:
        raise ValueError('Invalid window: incorrect window length')
    if not isinstance(window, numpy.ndarray) and window not in window_map:
        raise ValueError('Invalid window: unknown window {!r}'.format(window))
    if avg_method not in ('mean', 'median', 'median-mean'):
        raise ValueError('Invalid averaging method')
    if type(seg_len) is not int or type(seg_stride) is not int \
        or seg_len <= 0 or seg_stride <= 0:
        raise ValueError('Segment length and stride must be positive integers')

    if timeseries.precision == 'single':
        fs_dtype = numpy.complex64
    elif timeseries.precision == 'double':
        fs_dtype = numpy.complex128

    num_samples = len(timeseries)
    if num_segments is None:
        num_segments = int(num_samples // seg_stride)
        # NOTE: Is this not always true?
        if (num_segments - 1) * seg_stride + seg_len > num_samples:
            num_segments -= 1

    if not require_exact_data_fit:
        data_len = (num_segments - 1) * seg_stride + seg_len

        # Get the correct amount of data
        if data_len < num_samples:
            diff = num_samples - data_len
            start = diff // 2
            end = num_samples - diff // 2
            # Want this to be integers so if diff is odd, catch it here.
            if diff % 2:
                start = start + 1

            timeseries = timeseries[start:end]
            num_samples = len(timeseries)
        if data_len > num_samples:
            err_msg = "I was asked to estimate a PSD on %d " %(data_len)
            err_msg += "data samples. However the data provided only contains "
            err_msg += "%d data samples." %(num_samples)

    if num_samples != (num_segments - 1) * seg_stride + seg_len:
        raise ValueError('Incorrect choice of segmentation parameters')

    if not isinstance(window, numpy.ndarray):
        window = window_map[window](seg_len)
    w = Array(window.astype(timeseries.dtype))

    # calculate psd of each segment
    delta_f = 1. / timeseries.delta_t / seg_len
    if not USE_CACHING_FOR_WELCH_FFTS:
        segment_tilde = FrequencySeries(
            numpy.zeros(int(seg_len / 2 + 1)),
            delta_f=delta_f,
            dtype=fs_dtype,
        )

    segment_psds = []
    for i in range(num_segments):
        segment_start = i * seg_stride
        segment_end = segment_start + seg_len
        segment = timeseries[segment_start:segment_end]
        assert len(segment) == seg_len
        if not USE_CACHING_FOR_WELCH_FFTS:
            fft(segment * w, segment_tilde)
        else:
            segment_tilde = execute_cached_fft(segment * w,
                                               uid=WELCH_UNIQUE_ID)
        seg_psd = abs(segment_tilde * segment_tilde.conj()).numpy()

        #halve the DC and Nyquist components to be consistent with TO10095
        seg_psd[0] /= 2
        seg_psd[-1] /= 2

        segment_psds.append(seg_psd)

    segment_psds = numpy.array(segment_psds)

    if avg_method == 'mean':
        psd = numpy.mean(segment_psds, axis=0)
    elif avg_method == 'median':
        psd = numpy.median(segment_psds, axis=0) / median_bias(num_segments)
    elif avg_method == 'median-mean':
        odd_psds = segment_psds[::2]
        even_psds = segment_psds[1::2]
        odd_median = numpy.median(odd_psds, axis=0) / \
            median_bias(len(odd_psds))
        even_median = numpy.median(even_psds, axis=0) / \
            median_bias(len(even_psds))
        psd = (odd_median + even_median) / 2

    w = w.numpy()
    psd *= 2 * delta_f * seg_len / (w*w).sum()

    return FrequencySeries(psd, delta_f=delta_f, dtype=timeseries.dtype,
                           epoch=timeseries.start_time)

def inverse_spectrum_truncation(psd, max_filter_len, which_spectrum='invasd',
                                low_frequency_cutoff=None, 
                                low_frequency_fill_value=0., trunc_method=None):
    """Modify a PSD such that the impulse response associated with its inverse
    square root is no longer than `max_filter_len` time samples. In practice
    this corresponds to a coarse graining or smoothing of the PSD.

    Parameters
    ----------
    psd : FrequencySeries
        PSD whose inverse spectrum is to be truncated.
    max_filter_len : int
        Maximum length of the time-domain filter in samples.
    which_spectrum : {'invasd', 'invpsd'}
        Which spectrum to truncate. If 'invasd' (default), apply truncation to
        the inverse ASD. If 'invpsd', apply to the inverse PSD.
    low_frequency_cutoff : {None, int}
        Frequencies below `low_frequency_cutoff` are set to value specified by 
        `low_frequency_fill_value`.
    low_frequency_fill_value : {float, 'fmin'}
        Value to set PSD to at frequencies below `low_frequency_cutoff`.
        Default 0. If 'fmin', set to the value of the PSD at the low frequency
        cutoff index.
    trunc_method : {None, 'hann'}
        Function used for truncating the time-domain filter.
        None produces a hard truncation at `max_filter_len`.
    

    Returns
    -------
    psd : FrequencySeries
        PSD whose inverse spectrum has been truncated.

    Raises
    ------
    ValueError
        For invalid types or values of `max_filter_len` and `low_frequency_cutoff`.

    Notes
    -----
    See arXiv:gr-qc/0509116 for details.
    """
    from pycbc.strain.strain import execute_cached_fft, execute_cached_ifft

    # sanity checks
    if type(max_filter_len) is not int or max_filter_len <= 0:
        raise ValueError('max_filter_len must be a positive integer')
    if low_frequency_cutoff is not None and \
            (low_frequency_cutoff < 0. or
             low_frequency_cutoff > psd.sample_frequencies[-1]):
        raise ValueError('low_frequency_cutoff must be within the bandwidth of '
                         'the PSD')

    N = (len(psd)-1)*2

    inv_spectrum = FrequencySeries(zeros(len(psd)), delta_f=psd.delta_f, \
        dtype=complex_same_precision_as(psd))

    kmin = 1
    if low_frequency_cutoff:
        kmin = int(low_frequency_cutoff / psd.delta_f)
    
    # set values below low frequency cutoff
    if low_frequency_fill_value != 0.:
        if low_frequency_fill_value == 'fmin':
            low_frequency_fill_value = 1./psd[kmin]
        inv_spectrum[:kmin] = float(low_frequency_fill_value)

    inv_spectrum[kmin:N//2] = (1.0 / psd[kmin:N//2])

    # if truncating asd, take sqrt
    if which_spectrum == 'invasd':
        inv_spectrum[:N//2] = inv_spectrum[:N//2]**0.5
    elif which_spectrum != 'invpsd':
        raise ValueError(f'Invalid which_spectrum input {which_spectrum}; '
                         f'input must be either "invpsd" or "invasd"')

    if not USE_CACHING_FOR_INV_SPEC_TRUNC:
        q = TimeSeries(
            numpy.zeros(N),
            delta_t=(N / psd.delta_f),
            dtype=real_same_precision_as(psd)
        )
        ifft(inv_spectrum, q)
    else:
        q = execute_cached_ifft(inv_spectrum, copy_output=False,
                                uid=INVSPECTRUNC_UNIQUE_ID)

    trunc_start = max_filter_len // 2
    trunc_end = N - max_filter_len // 2
    if trunc_end < trunc_start:
        raise ValueError('Invalid value in inverse_spectrum_truncation')

    if trunc_method == 'hann':
        trunc_window = Array(numpy.hanning(max_filter_len), dtype=q.dtype)
        q[0:trunc_start] *= trunc_window[-trunc_start:]
        q[trunc_end:N] *= trunc_window[0:max_filter_len//2]

    if trunc_start < trunc_end:
        q[trunc_start:trunc_end] = 0
    if not USE_CACHING_FOR_INV_SPEC_TRUNC:
        psd_trunc = FrequencySeries(
            numpy.zeros(len(psd)),
            delta_f=psd.delta_f,
            dtype=complex_same_precision_as(psd)
        )
        fft(q, psd_trunc)
    else:
        psd_trunc = execute_cached_fft(q, copy_output=False,
                                       uid=INVSPECTRUNC_UNIQUE_ID)
    if which_spectrum == 'invasd':
        psd_trunc *= psd_trunc.conj()
    psd_out = 1. / abs(psd_trunc)

    return psd_out

def interpolate(series, delta_f, length=None):
    """Return a new PSD that has been interpolated to the desired delta_f.

    Parameters
    ----------
    series : FrequencySeries
        Frequency series to be interpolated.
    delta_f : float
        The desired delta_f of the output
    length : None or int
        The desired number of frequency samples. The default is None,
        so it will be calculated from the given `series` and `delta_f`.
        But this will cause an inconsistency issue of length sometimes,
        so if `length` is given, then just use it.

    Returns
    -------
    interpolated series : FrequencySeries
        A new FrequencySeries that has been interpolated.
    """
    if length is None:
        new_n = (len(series)-1) * series.delta_f / delta_f + 1
    else:
        new_n = length
    samples = numpy.arange(0, numpy.rint(new_n)) * delta_f
    interpolated_series = numpy.interp(samples, series.sample_frequencies.numpy(), series.numpy())
    return FrequencySeries(interpolated_series, epoch=series.epoch,
                           delta_f=delta_f, dtype=series.dtype)

def bandlimited_interpolate(series, delta_f):
    """Return a new PSD that has been interpolated to the desired delta_f.

    Parameters
    ----------
    series : FrequencySeries
        Frequency series to be interpolated.
    delta_f : float
        The desired delta_f of the output

    Returns
    -------
    interpolated series : FrequencySeries
        A new FrequencySeries that has been interpolated.
    """
    series = FrequencySeries(series, dtype=complex_same_precision_as(series), delta_f=series.delta_f)

    N = (len(series) - 1) * 2
    delta_t = 1.0 / series.delta_f / N

    new_N = int(1.0 / (delta_t * delta_f))
    new_n = new_N // 2 + 1

    series_in_time = TimeSeries(zeros(N), dtype=real_same_precision_as(series), delta_t=delta_t)
    ifft(series, series_in_time)

    padded_series_in_time = TimeSeries(zeros(new_N), dtype=series_in_time.dtype, delta_t=delta_t)
    padded_series_in_time[0:N//2] = series_in_time[0:N//2]
    padded_series_in_time[new_N-N//2:new_N] = series_in_time[N//2:N]

    interpolated_series = FrequencySeries(zeros(new_n), dtype=series.dtype, delta_f=delta_f)
    fft(padded_series_in_time, interpolated_series)

    return interpolated_series


def trimmed_welch(timeseries, seg_len=4096, seg_stride=2048, alpha=0.20):
    """Robust Trimmed-Mean Welch PSD estimator.
    Discards the lowest alpha and highest alpha fraction of periodograms at each bin,
    mitigating both transient loud glitches/signals and deep nulls.
    Corrects for the finite-sample trimmed mean bias on Chi-squared (2 DOF) noise.

    Parameters
    ----------
    timeseries : TimeSeries
        Input PyCBC TimeSeries.
    seg_len : int, optional
        Segment length in samples (default: 4096).
    seg_stride : int, optional
        Segment stride in samples (default: 2048).
    alpha : float, optional
        Fraction to trim from each end (default: 0.20).

    Returns
    -------
    psd : FrequencySeries
        Bias-corrected estimated PSD.
    """
    fs = float(timeseries.sample_rate)
    delta_f = 1.0 / (seg_len / fs)
    num_segs = (len(timeseries) - seg_len) // seg_stride + 1
    if num_segs < 1:
        raise ValueError("timeseries too short for seg_len")
    w = numpy.hanning(seg_len)
    w_sum_sq = numpy.sum(w**2)

    seg_specs = []
    for i in range(num_segs):
        seg = timeseries[i * seg_stride : i * seg_stride + seg_len].numpy()
        seg_spec = numpy.abs(numpy.fft.rfft(w * seg))**2 * (2.0 / (fs * w_sum_sq))
        # Halve DC and Nyquist components to be consistent with T010095 & PyCBC convention
        seg_spec[0] /= 2.0
        seg_spec[-1] /= 2.0
        seg_specs.append(seg_spec)
    seg_specs = numpy.array(seg_specs)

    k_trim = int(num_segs * alpha)
    sorted_specs = numpy.sort(seg_specs, axis=0)
    if k_trim > 0 and 2 * k_trim < num_segs:
        trimmed = sorted_specs[k_trim : num_segs - k_trim, :]
        raw_psd = numpy.mean(trimmed, axis=0)
        # Analytical expectation of trimmed standard exponential distribution:
        bias_factor = float(1.0 + (alpha * numpy.log(alpha) - (1.0 - alpha) * numpy.log(1.0 - alpha)) / (1.0 - 2.0 * alpha))
    else:
        raw_psd = numpy.median(sorted_specs, axis=0)
        bias_factor = median_bias(num_segs)

    corrected_psd = raw_psd / bias_factor
    return FrequencySeries(corrected_psd, delta_f=delta_f, epoch=timeseries.start_time)


def multitaper(timeseries, seg_len=4096, seg_stride=2048, NW=3.0, avg_method='median'):
    """Multitaper Spectral Estimator using Discrete Prolate Spheroidal Sequences (DPSS).
    Applies K = 2*NW - 1 orthonormal Slepian tapers per segment to minimize spectral
    leakage across steep dynamic range gradients (e.g. seismic wall at 15-20 Hz).
    Combines segments via median or mean for non-Gaussian glitch robustness.

    Parameters
    ----------
    timeseries : TimeSeries
        Input PyCBC TimeSeries.
    seg_len : int, optional
        Segment length in samples (default: 4096).
    seg_stride : int, optional
        Segment stride in samples (default: 2048).
    NW : float, optional
        Time-halfbandwidth product (default: 3.0).
    avg_method : {'median', 'mean'}, optional
        Segment averaging method (default: 'median').

    Returns
    -------
    psd : FrequencySeries
        Bias-corrected estimated PSD.
    """
    from scipy.signal.windows import dpss
    import scipy.stats as stats

    fs = float(timeseries.sample_rate)
    delta_f = 1.0 / (seg_len / fs)
    K = int(2 * NW - 1)
    tapers, _ = dpss(seg_len, NW=NW, Kmax=K, return_ratios=True)

    num_segs = (len(timeseries) - seg_len) // seg_stride + 1
    if num_segs < 1:
        raise ValueError("timeseries too short for seg_len")
    seg_psds = []
    for i in range(num_segs):
        seg = timeseries[i * seg_stride : i * seg_stride + seg_len].numpy()
        tapered = tapers * seg[None, :]
        spec = numpy.abs(numpy.fft.rfft(tapered, axis=1))**2 * (2.0 / fs)
        # Halve DC and Nyquist components to be consistent with T010095 & PyCBC convention
        spec[:, 0] /= 2.0
        spec[:, -1] /= 2.0
        seg_psd = numpy.mean(spec, axis=0)
        seg_psds.append(seg_psd)
    seg_psds = numpy.array(seg_psds)

    if avg_method == 'median':
        psd_raw = numpy.median(seg_psds, axis=0)
        bias = float(stats.gamma.ppf(0.5, a=K) / K)
        psd_est = psd_raw / bias
    elif avg_method == 'mean':
        psd_est = numpy.mean(seg_psds, axis=0)
    else:
        raise ValueError(f"Unknown avg_method: {avg_method}")

    return FrequencySeries(psd_est, delta_f=delta_f, epoch=timeseries.start_time)


# Backwards compatibility aliases
estimate_psd_trimmed_welch = trimmed_welch
estimate_psd_multitaper = multitaper



