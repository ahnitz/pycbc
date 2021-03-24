""" PSD Variation """

import numpy
import logging
import scipy.signal as sig
import pycbc.psd
from numpy.fft import rfft, irfft
from pycbc.types import TimeSeries, complex_same_precision_as
from pycbc.filter import resample_to_delta_t, sigma


def mean_square(data, delta_t, srate, short_stride, stride):
    """ Calculate mean square of given time series once per stride

    First of all this function calculate the mean square of given time
    series once per short_stride. This is used to find and remove
    outliers due to short glitches. Here an outlier is defined as any
    element which is greater than two times the average of its closest
    neighbours. Every outlier is substituted with the average of the
    corresponding adjacent elements.
    Then, every second the function compute the mean square of the
    smoothed time series, within the stride.

    Parameters
    ----------
    data : numpy.ndarray
    delta_t : float
        Duration of the time series
    srate : int
        Sample rate of the data were it given as a TimeSeries
    short_stride : float
        Stride duration for outlier removal
    stride ; float
        Stride duration

    Returns
    -------
    m_s: List
        Mean square of given time series
    """

    # Calculate mean square of data once per short stride and replace
    # outliers
    short_ms = numpy.mean(data.reshape(-1, int(srate * short_stride)) ** 2,
                          axis=1)
    # Define an array of averages that is used to substitute outliers
    ave = 0.5 * (short_ms[2:] + short_ms[:-2])
    outliers = short_ms[1:-1] > (2. * ave)
    short_ms[1:-1][outliers] = ave[outliers]

    # Calculate mean square of data every step within a window equal to
    # stride seconds
    m_s = []
    inv_time = int(1. / short_stride)
    for index in range(int(delta_t - stride + 1)):
        m_s.append(numpy.mean(short_ms[inv_time * index:inv_time *
                                       int(index+stride)]))
    return m_s


def calc_filt_psd_variation(strain, segment, short_segment, psd_long_segment,
                            psd_duration, psd_stride, psd_avg_method, low_freq,
                            high_freq):
    """ Calculates time series of PSD variability

    This function first splits the segment up into 512 second chunks. It
    then calculates the PSD over this 512 second. The PSD is used to
    to create a filter that is the composition of three filters:
    1. Bandpass filter between f_low and f_high.
    2. Weighting filter which gives the rough response of a CBC template.
    3. Whitening filter.
    Next it makes the convolution of this filter with the stretch of data.
    This new time series is given to the "mean_square" function, which
    computes the mean square of the timeseries within an 8 seconds window,
    once per second.
    The result, which is the variance of the S/N in that stride for the
    Parseval theorem, is then stored in a timeseries.

    Parameters
    ----------
    strain : TimeSeries
        Input strain time series to estimate PSDs
    segment : {float, 8}
        Duration of the segments for the mean square estimation in seconds.
    short_segment : {float, 0.25}
        Duration of the short segments for the outliers removal.
    psd_long_segment : {float, 512}
        Duration of the long segments for PSD estimation in seconds.
    psd_duration : {float, 8}
        Duration of FFT segments for long term PSD estimation, in seconds.
    psd_stride : {float, 4}
        Separation between FFT segments for long term PSD estimation, in
        seconds.
    psd_avg_method : {string, 'median'}
        Method for averaging PSD estimation segments.
    low_freq : {float, 20}
        Minimum frequency to consider the comparison between PSDs.
    high_freq : {float, 480}
        Maximum frequency to consider the comparison between PSDs.

    Returns
    -------
    psd_var : TimeSeries
        Time series of the variability in the PSD estimation
    """
    # Calculate strain precision
    if strain.precision == 'single':
        fs_dtype = numpy.float32
    elif strain.precision == 'double':
        fs_dtype = numpy.float64

    # Convert start and end times immediately to floats
    start_time = numpy.float(strain.start_time)
    end_time = numpy.float(strain.end_time)

    # Resample the data
    strain = resample_to_delta_t(strain, 1.0 / 2048)
    srate = int(strain.sample_rate)

    # Fix the step for the PSD estimation and the time to remove at the
    # edge of the time series.
    step = 1.0
    strain_crop = 8.0

    # Find the times of the long segments
    times_long = numpy.arange(start_time, end_time,
                              psd_long_segment - 2 * strain_crop
                              - segment + step)

    # Create a bandpass filter between low_freq and high_freq
    filt = sig.firwin(4 * srate, [low_freq, high_freq], pass_zero=False,
                      window='hann', nyq=srate / 2)
    filt.resize(int(psd_duration * srate))
    # Fourier transform the filter and take the absolute value to get
    # rid of the phase.
    filt = abs(rfft(filt))

    psd_var_list = []
    for tlong in times_long:
        # Calculate PSD for long segment
        if tlong + psd_long_segment <= float(end_time):
            astrain = strain.time_slice(tlong, tlong + psd_long_segment)
            plong = pycbc.psd.welch(
                astrain,
                seg_len=int(psd_duration * strain.sample_rate),
                seg_stride=int(psd_stride * strain.sample_rate),
                avg_method=psd_avg_method)
        else:
            astrain = strain.time_slice(tlong, end_time)
            plong = pycbc.psd.welch(
                           strain.time_slice(end_time - psd_long_segment,
                                             end_time),
                           seg_len=int(psd_duration * strain.sample_rate),
                           seg_stride=int(psd_stride * strain.sample_rate),
                           avg_method=psd_avg_method)
        astrain = astrain.numpy()
        freqs = numpy.array(plong.sample_frequencies, dtype=fs_dtype)
        plong = plong.numpy()

        # Make the weighting filter - bandpass, which weight by f^-7/6,
        # and whiten. The normalization is chosen so that the variance
        # will be one if this filter is applied to white noise which
        # already has a variance of one.
        fweight = freqs ** (-7./6.) * filt / numpy.sqrt(plong)
        fweight[0] = 0.
        norm = (sum(abs(fweight) ** 2) / (len(fweight) - 1.)) ** -0.5
        fweight = norm * fweight
        fwhiten = numpy.sqrt(2. / srate) / numpy.sqrt(plong)
        fwhiten[0] = 0.
        full_filt = sig.hann(int(psd_duration * srate)) * numpy.roll(
            irfft(fwhiten * fweight), int(psd_duration / 2) * srate)
        # Convolve the filter with long segment of data
        wstrain = sig.fftconvolve(astrain, full_filt, mode='same')
        wstrain = wstrain[int(strain_crop * srate):-int(strain_crop * srate)]
        # compute the mean square of the chunk of data
        delta_t = len(wstrain) * strain.delta_t
        variation = mean_square(wstrain, delta_t, srate, short_segment, segment)
        psd_var_list.append(numpy.array(variation, dtype=wstrain.dtype))

    # Package up the time series to return
    psd_var = TimeSeries(numpy.concatenate(psd_var_list), delta_t=step,
                         epoch=start_time + strain_crop + segment)

    return psd_var


def find_trigger_value(psd_var, idx, start, sample_rate):
    """ Find the PSD variation value at a particular time with the filter
    method. If the time is outside the timeseries bound, 1. is given.

    Parameters
    ----------
    psd_var : TimeSeries
        Time series of the varaibility in the PSD estimation
    idx : numpy.ndarray
        Time indices of the triggers
    start : float
        GPS start time
    sample_rate : float
        Sample rate defined in ini file

    Returns
    -------
    vals : Array
        PSD variation value at a particular time
    """
    # Find gps time of the trigger
    time = start + idx / sample_rate
    # Extract the PSD variation at trigger time through linear
    # interpolation
    if not hasattr(psd_var, 'cached_psd_var_interpolant'):
        from scipy import interpolate
        psd_var.cached_psd_var_interpolant = \
            interpolate.interp1d(psd_var.sample_times.numpy(), psd_var.numpy(),
                                 fill_value=1.0, bounds_error=False)
    vals = psd_var.cached_psd_var_interpolant(time)



def sigma_var(strain, start, end, step=2.0, duration=16, subsegment=1, fmin=20,
            fmax=480, mass=10.0, time_offset=1.0):
    """ Calculate how much the noise weighted reference signal amplitude
    varies compared to the average.

    Parameters
    ----------
    strain: pycbc.type.TimeSeries
        Data containing the requrested start / end times to calculate the
        sigma variation for.
    start: float
        Earliest time allowed to be used.
    end: float
        Latest time allowed to be used
    step: {float, 2.0}
        The duration in seconds between samples of the sigma value.
    duration: {float, 16.0}
        The time in seconds used to calculate a local PSD estimate
    subsegment: {float, 1}
        The length of a psd sample when calculating the local PSD estimate
    fmin: {float, 20}
        The minimum frequency to use in the sigma calculation
    fmax: {float, 480}
        The maximum frequency to us in the sigma calculation
    mass: {float, 10.0}
        The component mass to use (assumed equal mass) for the reference
        signal.
    time_offset: {float, 1.0}
        Time offset in seconds to add to the assigned times for each sigma
        value. If zero, the time chosen is the central time used for the local
        psd estimate.

    Returns
    --------
    scale: scipy.interpolate.interp1d
        Interpolator that returns the sigma scaling for a given time. Assumed
        to be 1.0 outside the possible interval.
    """
    from pycbc.waveform import get_fd_waveform
    from scipy.interpolate import interp1d

    # segment data and calculate psd in each short period
    tstarts  = numpy.arange(start, end - duration, step)
    tends = tstarts + duration
    pshorts = [strain.time_slice(s, e).psd(subsegment)
               for s, e in zip(tstarts, tends)]

    # Calculate the scale factors
    m1 = m2 = mass
    href, _ = get_fd_waveform(approximant="IMRPhenomD", mass1=m1, mass2=m2,
                              f_lower=fmin, delta_f=pshorts[0].delta_f)
    href = href.astype(complex_same_precision_as(pshorts[0]))
    sigmas = [sigma(href, psd=p, low_frequency_cutoff=fmin,
                    high_frequency_cutoff=fmax) for p in pshorts]
    sigmas = numpy.array(sigmas)
    sigmas /= numpy.median(sigmas)

    # choose central time of each short sigma estimate offset by preset value
    sigtime = ((tends + tstarts) / 2.0 + time_offset).astype(numpy.float64)
    scale = interp1d(sigtime, sigmas, fill_value=1.0, bounds_error=False)
    return scale


class SigmaVariation(object):
    """Class that encapsulates handling of sigma variation calculation within
    Inspiral-type executables.
    """
    def __init__(self, segments, durations, subsegments, time_offsets, step,
                 fmin, fmax, mass): 
        self.sample_rate = segments.strain.sample_rate
                 
        # Set the naming and return format of the saved stats
        self._returns = {}
        self.names = {}
        for dur in durations:
            name = 'sigvar_{}'.format(int(dur))
            self.names[dur] = name
            self._returns[name] = numpy.float32
                 
        # precalculate all the sigma variation interpolating functions
        self.scales = []
        for i, seg in enumerate(segments.segment_slices):
            scale = {}
            group = zip(durations, subsegments, time_offsets)
            for dur, subs, toff in group:
                logging.info('precalc sigma var for seg %s, dur %s', i, dur)
                segd = segments.strain[seg].copy()
                segd.start_time = 0 # avoids doing time accounting when applying
                scale[dur] = sigma_var(segd, 0, segd.end_time,
                                       fmin=fmin, fmax=fmax,
                                       step=step,
                                       duration=dur,
                                       subsegment=subs,
                                       time_offset=toff)                                  
            self.scales.append(scale)
            
    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Sigma Variation")
        group.add_argument("--sigvar-durations", nargs='+', type=float, 
                          help="Durations of time used for local psd estimates")
        group.add_argument("--sigvar-subsegments", nargs='+', type=float, 
                          help=" The length of a psd sample when calculating"
                               " the local PSD estimate")
        group.add_argument("--sigvar-time-offsets", nargs='+', type=float, 
                          help="Seconds to offset each estimates time")
        group.add_argument("--sigvar-step", type=float, 
                          help="Seconds between estimates")                        
        group.add_argument("--sigvar-fmin", type=float, 
                          help="low frequency cutoff")
        group.add_argument("--sigvar-fmax", type=float, 
                          help="high frequency cutoff")
        group.add_argument("--sigvar-mass", type=float, 
                          help="low frequency cutoff")
                               
    @classmethod
    def from_cli(cls, args, segments):
        if args.sigvar_durations is None:
            return None
        return cls(segments, 
                   args.sigvar_durations,
                   args.sigvar_subsegments,
                   args.sigvar_time_offsets,
                   args.sigvar_fmin,
                   args.sigvar_step,
                   args.sigvar_fmax,
                   args.sigvar_mass)

    @property
    def returns(self):
        return self._returns

    def values(self, seg_num, idx):
        out_vals = {}
        for dur in self.scales[seg_num]:
            name = self.names[dur]
            times = idx / float(self.sample_rate)
            out_vals[name] = self.scales[seg_num][dur](times)
        return out_vals

