# Copyright (C) 2015 Alex Nitz
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This modules contains functions for calculating and manipulating
coincident triggers.
"""
import numpy, logging, h5py, pycbc.pnutils, copy, lal
from itertools import izip
from scipy.interpolate import interp1d

def background_bin_from_string(background_bins, data):
    """ Return template ids for each bin as defined by the format string

    Parameters
    ----------
    bins: list of strings
        List of strings which define how a background bin is taken from the
        list of templates.
    data: dict of numpy.ndarrays
        Dict with parameter key values and numpy.ndarray values which define
        the parameters of the template bank to bin up.

    Returns
    -------
    bins: dict
        Dictionary of location indices indexed by a bin name
    """
    used = numpy.array([], dtype=numpy.uint32)
    bins = {}
    for mbin in background_bins:
        name, bin_type, boundary = tuple(mbin.split(':'))

        if boundary[0:2] == 'lt':
            member_func = lambda vals, bd=boundary : vals < float(bd[2:])
        elif boundary[0:2] == 'gt':
            member_func = lambda vals, bd=boundary : vals > float(bd[2:])
        else:
            raise RuntimeError("Can't parse boundary condition! Must begin "
                               "with 'lt' or 'gt'")

        if bin_type == 'component' and boundary[0:2] == 'lt':
            # maximum component mass is less than boundary value
            vals = numpy.maximum(data['mass1'], data['mass2'])
        if bin_type == 'component' and boundary[0:2] == 'gt':
            # minimum component mass is greater than bdary
            vals = numpy.minimum(data['mass1'], data['mass2'])
        elif bin_type == 'total':
            vals = data['mass1'] + data['mass2']
        elif bin_type == 'chirp':
            vals = pycbc.pnutils.mass1_mass2_to_mchirp_eta(
                                               data['mass1'], data['mass2'])[0]
        elif bin_type == 'SEOBNRv2Peak':
            vals = pycbc.pnutils.get_freq('fSEOBNRv2Peak',
                  data['mass1'], data['mass2'], data['spin1z'], data['spin2z'])
        elif bin_type == 'duration':
            # pnutils vectorized function shadowing pycbc.waveform
            vals = pycbc.pnutils.get_seobnrrom_duration(data['mass1'],
                data['mass2'], data['spin1z'], data['spin2z'], data['f_lower'])
        else:
            raise ValueError('Invalid bin type %s' % bin_type)

        locs = member_func(vals)
        del vals

        # make sure we don't reuse anything from an earlier bin
        locs = numpy.where(locs)[0]
        locs = numpy.delete(locs, numpy.where(numpy.in1d(locs, used))[0])
        used = numpy.concatenate([used, locs])
        bins[name] = locs

    return bins


def calculate_n_louder(bstat, fstat, dec, skip_background=False):
    """ Calculate for each foreground event the number of background events
    that are louder than it.

    Parameters
    ----------
    bstat: numpy.ndarray
        Array of the background statistic values
    fstat: numpy.ndarray
        Array of the foreground statitsic values
    dec: numpy.ndarray
        Array of the decimation factors for the background statistics
    skip_background: optional, {boolean, False}
        Skip calculating cumulative numbers for background triggers

    Returns
    -------
    cum_back_num: numpy.ndarray
        The cumulative array of background triggers. Does not return this
        argument if skip_background == True
    fore_n_louder: numpy.ndarray
        The number of background triggers above each foreground trigger
    """
    sort = bstat.argsort()
    bstat = bstat[sort]
    dec = dec[sort]

    # calculate cumulative number of triggers louder than the trigger in
    # a given index. We need to subtract the decimation factor, as the cumsum
    # includes itself in the first sum (it is inclusive of the first value)
    n_louder = dec[::-1].cumsum()[::-1] - dec

    # Determine how many values are louder than the foreground ones
    # We need to subtract one from the index, to be consistent with the definition
    # of n_louder, as here we do want to include the background value at the
    # found index
    idx = numpy.searchsorted(bstat, fstat, side='left') - 1

    # If the foreground are *quieter* than the background or at the same value
    # then the search sorted alorithm will choose position -1, which does not exist
    # We force it back to zero.
    if isinstance(idx, numpy.ndarray): # Handle the case where our input is an array
        idx[idx < 0] = 0
    else: # Handle the case where we are simply given a scalar value
        if idx < 0:
            idx = 0

    fore_n_louder = n_louder[idx]

    if not skip_background:
        unsort = sort.argsort()
        back_cum_num = n_louder[unsort]
        return back_cum_num, fore_n_louder
    else:
        return fore_n_louder

def timeslide_durations(start1, start2, end1, end2, timeslide_offsets):
    """ Find the coincident time for each timeslide.

    Find the coincident time for each timeslide, where the first time vector
    is slid to the right by the offset in the given timeslide_offsets vector.

    Parameters
    ----------
    start1: numpy.ndarray
        Array of the start of valid analyzed times for detector 1
    start2: numpy.ndarray
        Array of the start of valid analyzed times for detector 2
    end1: numpy.ndarray
        Array of the end of valid analyzed times for detector 1
    end2: numpy.ndarray
        Array of the end of valid analyzed times for detector 2
    timseslide_offset: numpy.ndarray
        Array of offsets (in seconds) for each timeslide

    Returns
    --------
    durations: numpy.ndarray
        Array of coincident time for each timeslide in the offset array
    """
    from . import veto
    durations = []
    seg2 = veto.start_end_to_segments(start2, end2)
    for offset in timeslide_offsets:
        seg1 = veto.start_end_to_segments(start1 + offset, end1 + offset)
        durations.append(abs((seg1 & seg2).coalesce()))
    return numpy.array(durations)

def time_coincidence(t1, t2, window, slide_step=0):
    """ Find coincidences by time window

    Parameters
    ----------
    t1 : numpy.ndarray
        Array of trigger times from the first detector
    t2 : numpy.ndarray
        Array of trigger times from the second detector
    window : float
        The coincidence window in seconds
    slide_step : optional, {None, float}
        If calculating background coincidences, the interval between background
        slides in seconds.

    Returns
    -------
    idx1 : numpy.ndarray
        Array of indices into the t1 array.
    idx2 : numpy.ndarray
        Array of indices into the t2 array.
    slide : numpy.ndarray
        Array of slide ids
    """
    if slide_step:
        fold1 = t1 % slide_step
        fold2 = t2 % slide_step
    else:
        fold1 = t1
        fold2 = t2

    sort1 = fold1.argsort()
    sort2 = fold2.argsort()
    fold1 = fold1[sort1]
    fold2 = fold2[sort2]

    if slide_step:
        fold2 = numpy.concatenate([fold2 - slide_step, fold2, fold2 + slide_step])
        sort2 = numpy.concatenate([sort2, sort2, sort2])

    left = numpy.searchsorted(fold2, fold1 - window)
    right = numpy.searchsorted(fold2, fold1 + window)

    idx1 = numpy.repeat(sort1, right-left)
    idx2 = [sort2[l:r] for l,r in zip(left, right)]

    if len(idx2) > 0:
        idx2 = numpy.concatenate(idx2)
    else:
        idx2 = numpy.array([])

    if slide_step:
        diff = ((t1 / slide_step)[idx1] - (t2 / slide_step)[idx2])
        slide = numpy.rint(diff)
    else:
        slide = numpy.zeros(len(idx1))

    return idx1.astype(numpy.uint32), idx2.astype(numpy.uint32), slide.astype(numpy.int32)


def cluster_coincs(stat, time1, time2, timeslide_id, slide, window, argmax=numpy.argmax):
    """Cluster coincident events for each timeslide separately, across
    templates, based on the ranking statistic

    Parameters
    ----------
    stat: numpy.ndarray
        vector of ranking values to maximize
    time1: numpy.ndarray
        first time vector
    time2: numpy.ndarray
        second time vector
    timeslide_id: numpy.ndarray
        vector that determines the timeslide offset
    slide: float
        length of the timeslides offset interval
    window: float
        length to cluster over

    Returns
    -------
    cindex: numpy.ndarray
        The set of indices corresponding to the surviving coincidences.
    """
    logging.info('clustering coinc triggers over %ss window' % window)

    if len(time1) == 0 or len(time2) == 0:
        logging.info('No coinc triggers in one, or both, ifos.')
        return numpy.array([])

    indices = []
    if numpy.isfinite(slide):
        time = (time2 + (time1 + timeslide_id * slide)) / 2
    else:
        time = 0.5 * (time2 + time1)

    tslide = timeslide_id.astype(numpy.float128)
    time = time.astype(numpy.float128)

    span = (time.max() - time.min()) + window * 10
    time = time + span * tslide

    time_sorting = time.argsort()
    stat = stat[time_sorting]
    time = time[time_sorting]

    logging.info('sorting...')
    left = numpy.searchsorted(time, time - window)
    right = numpy.searchsorted(time, time + window)
    logging.info('done sorting')
    indices = numpy.zeros(len(left), dtype=numpy.uint32)

    # i is the index we are inspecting, j is the next one to save
    i = 0
    j = 0
    while i < len(left):
        l = left[i]
        r = right[i]

        # If there are no other points to compare it is obviously the max
        if (r - l) == 1:
            indices[j] = i
            j += 1
            i += 1
            continue

        # Find the location of the maximum within the time interval around i
        max_loc = argmax(stat[l:r]) + l

        # If this point is the max, we can skip to the right boundary
        if max_loc == i:
            indices[j] = i
            i = r
            j += 1

        # If the max is later than i, we can skip to it
        elif max_loc > i:
            i = max_loc

        elif max_loc < i:
            i += 1

    indices = indices[:j]

    logging.info('done clustering coinc triggers: %s triggers remaining' % len(indices))
    return time_sorting[indices]

class MultiRingBuffer(object):
    def __init__(self, num_rings, max_length, dtype=numpy.float32):
        self.max_length = max_length
        self.pad_count = self.max_length + 2
        self.buffer = numpy.zeros((num_rings, self.pad_count), dtype=dtype)

        self.buffer_expire = numpy.zeros((num_rings, self.pad_count), dtype=numpy.int32) 
        self.buffer_expire -= self.max_length * 2

        self.start = numpy.zeros(num_rings, dtype=numpy.uint32)
        self.index = numpy.zeros(num_rings, dtype=numpy.uint32)
        self.ladder = numpy.arange(0, num_rings, dtype=numpy.uint32)    

        self.size = 0
        self.expire = 0

    def __len__(self):
        return self.size

    def num_elements(self):
        count = self.index - self.start
        count[self.index < self.start] += self.pad_count
        total = count.sum()
        return total

    def size_increment(self):
        if self.size < self.max_length:
            self.size += 1
        self.expire += 1

        idx = self.buffer_expire[self.ladder, self.start] < self.expire - self.max_length
        self.start[numpy.logical_and(idx, self.start != self.index)] += 1
        self.start[self.start == self.pad_count] = 0

    def add(self, indices, values):
        index = self.index[indices]

        self.buffer[indices, index] = values
        self.buffer_expire[indices, index] = self.expire

        index += 1
        index[index == self.pad_count] = 0
        self.index[indices] = index
        self.size_increment()        

    def expire_vector(self, buffer_index):
        buffer_part = self.buffer_expire[buffer_index]
        start = self.start[buffer_index]
        end = self.index[buffer_index]
        
        if start <= end:
            return buffer_part[start:end]
        else:
            return numpy.concatenate([buffer_part[start:], buffer_part[:end]])        
    
    def data(self, buffer_index):
        buffer_part = self.buffer[buffer_index]
        start = self.start[buffer_index]
        end = self.index[buffer_index]
        
        if start <= end:
            return buffer_part[start:end]
        else:
            return numpy.concatenate([buffer_part[start:], buffer_part[:end]])

class CoincExpireBuffer(object):
    def __init__(self, expiration, ifos,
                       initial_size=2**20, dtype=numpy.float32):
        self.expiration = expiration
        self.buffer = numpy.zeros(initial_size, dtype=dtype)
        self.index = 0
        self.ifos = ifos

        self.time = {}
        self.timer = {}
        for ifo in self.ifos:
            self.time[ifo] = 0
            self.timer[ifo] = numpy.zeros(initial_size, dtype=numpy.int32)

    def add(self, values, times, ifos):
        # Resize the internal buffer if we need more space
        if self.index + len(values) >= len(self.buffer):
            newlen = len(self.buffer) * 2
            self.timer.resize(newlen)
            self.buffer.resize(newlen)

        self.buffer[self.index:self.index+len(values)] = values
        for ifo in self.ifos:
            self.timer[ifo][self.index:self.index+len(values)] = times[ifo]

        self.index += len(values)

        # Remove the expired old elements
        keep = None
        for ifo in ifos:
            kt = self.timer[ifo][:self.index] > self.time[ifo] - self.expiration
            keep = numpy.logical_and(keep, kt) if keep is not None else kt

        self.index = keep.sum()
        self.buffer[:self.index] = self.buffer[keep]

        for ifo in self.ifos:
            self.timer[ifo][:self.index] = self.timer[ifo][keep]

        for ifo in ifos:
            self.time[ifo] += 1

    def num_greater(self, value):
        return (self.buffer[:self.index] > value).sum()
        
    @property    
    def data(self):
        return self.buffer[:self.index]

class LiveCoincTimeslideBackgroundEstimator(object):
    def __init__(self, num_templates, analysis_block, background_statistic, stat_files, ifos, 
                 ifar_limit=100,
                 timeslide_interval=.035,
                 coinc_threshold=0.002, return_background=False):
        from pycbc import detector
        from . import stat
        self.analysis_block = analysis_block
        self.stat_calculator = stat.get_statistic(background_statistic)(stat_files)
        self.timeslide_interval = timeslide_interval
        self.return_background = return_background

        self.ifos = ifos
        if len(self.ifos) != 2:
            raise ValueError("Only a two ifo analysis is supported at this time")

        self.lookback_time = (ifar_limit * lal.YRJUL_SI * timeslide_interval) ** 0.5
        self.buffer_size = int(numpy.ceil(self.lookback_time / analysis_block))

        det0, det1 = detector.Detector(ifos[0]), detector.Detector(ifos[1])
        self.time_window = det0.light_travel_time_to_detector(det1) + coinc_threshold

        # Preallocate ring buffers to hold all the single detector triggers
        # from each template separately. The ring buffer naturally expires
        # old triggers when they wraparound.
        self.singles = {}
        self.singles_dtype = [('time', numpy.float64), ('stat', numpy.float32)]
        for ifo in self.ifos:
            self.singles[ifo] = MultiRingBuffer(num_templates,
                                                self.buffer_size,
                                                dtype=self.singles_dtype)
        self.coincs = CoincExpireBuffer(self.buffer_size, self.ifos)

    @property
    def background_time(self):
        time = 1.0 / self.timeslide_interval
        for ifo in self.singles:
            time *= len(self.singles[ifo]) * self.analysis_block
        return time

    def ifar(self, coinc_stat):
        """ Return the far that would be associated with the coincident given.
        """
        n = self.coincs.num_greater(coinc_stat)
        return self.background_time / lal.YRJUL_SI / (n + 1)

    def add_singles(self, results):
        # convert to single detector trigger values 
        # FIXME Currently configured to use pycbc live output 
        # where chisq is the reduced chisq and chisq_dof is the actual DOF
        logging.info("adding singles to the background estimate...")
        singles_data = {}
        for ifo in results:
            trigs = results[ifo]
            trigs = copy.copy(trigs)
            trigs['chisq'] = trigs['chisq'] * trigs['chisq_dof']
            trigs['chisq_dof'] = (trigs['chisq_dof'] + 2) / 2

            if len(trigs['snr'] > 0):
                single_stat = self.stat_calculator.single(trigs)
            else:
                single_stat = numpy.array([], ndmin=1, dtype=numpy.float32)

            # prepare the data that will go in the single detector buffers
            data = numpy.array(numpy.zeros(len(single_stat),
                               dtype=self.singles_dtype), ndmin=1)
            data['stat'] = single_stat
            data['time'] = trigs['end_time']
            singles_data[ifo] = (trigs['template_id'], data)

            # add each single detector trigger to the 
            # ring buffer associated with 
            # the template it was found in.
            self.singles[ifo].add(*singles_data[ifo])
            
        # for each single detector trigger find the allowed coincidences
        cstat = []
        offsets = []
        ctimes = {self.ifos[0]:[], self.ifos[1]:[]}
        single_block = {self.ifos[0]:[], self.ifos[1]:[]}
        for ifo in results:
            trigs = results[ifo]
            for i in range(len(trigs['end_time'])):
                trig_stat = singles_data[ifo][1]['stat'][i]
                trig_time = trigs['end_time'][i]
                template = trigs['template_id'][i]

                oifo = self.ifos[1] if self.ifos[0] == ifo else self.ifos[0]
                times = self.singles[oifo].data(template)['time']
                stats = self.singles[oifo].data(template)['stat']

                end_times = {ifo:numpy.array(trig_time, ndmin=1),
                            oifo:times}             

                i1, i2, slide = time_coincidence(end_times[self.ifos[0]],
                                                 end_times[self.ifos[1]], self.time_window,
                                                 self.timeslide_interval)
                idx = i1 if self.ifos[0] == oifo else i2

                c = self.stat_calculator.coinc(trig_stat, stats[idx],
                                               slide, self.timeslide_interval)
                offsets.append(slide)
                cstat.append(c)
                ctimes[oifo].append(times[idx])
                ctimes[ifo].append(numpy.zeros(len(c), dtype=numpy.float64))
                ctimes[ifo][-1].fill(trig_time)        

                single_block[oifo].append(self.singles[oifo].expire_vector(template)[idx])
                single_block[ifo].append(numpy.zeros(len(c), dtype=numpy.float64))
                single_block[ifo][-1].fill(self.singles[ifo].expire - 1)
   

        # cluster the triggers we've found
        # (both zerolag and non handled together)
        num_zerolag = 0
        if len(cstat) > 0:
            cstat = numpy.concatenate(cstat)
            logging.info('%s background and zerolag coincs', len(cstat))
            if len(cstat) > 0:
                offsets = numpy.concatenate(offsets)
                ctime0 = numpy.concatenate(ctimes[self.ifos[0]]).astype(numpy.float64)
                ctime1 = numpy.concatenate(ctimes[self.ifos[1]]).astype(numpy.float64)
                
                cidx = cluster_coincs(cstat, ctime0, ctime1, offsets, 
                                          self.timeslide_interval,
                                          self.analysis_block)

                cstat = cstat[cidx]
                offsets = offsets[cidx]
                ctime0 = ctime0[cidx]
                ctime1 = ctime1[cidx]

                zerolag_idx = (offsets == 0)
                bkg_idx = (offsets != 0)

                for ifo in self.ifos:
                    single_block[ifo] = numpy.concatenate(single_block[ifo])[cidx][bkg_idx]

                self.coincs.add(cstat[bkg_idx], single_block, results.keys())
                num_zerolag = zerolag_idx.sum()

        coinc_results = {}
        # Collect coinc results for saving
        if num_zerolag > 0:
            zerolag_results = {}
            zerolag_ctime0 = ctime0[zerolag_idx]
            zerolag_ctime1 = ctime1[zerolag_idx]

            zerolag_cstat = cstat[zerolag_idx]
            ifar = numpy.array([self.ifar(c) for c in zerolag_cstat], dtype=numpy.float32)
            tid0 = numpy.array([numpy.where(t == results[self.ifos[0]]['end_time'])[0]
                                for t in zerolag_ctime0], dtype=numpy.uint32)                    
            tid1 = numpy.array([numpy.where(t == results[self.ifos[1]]['end_time'])[0]
                                for t in zerolag_ctime1], dtype=numpy.uint32)
            zerolag_results['foreground/ifar'] = ifar
            zerolag_results['foreground/stat'] = zerolag_cstat
            zerolag_results['foreground/%s/end_time' % self.ifos[0]] = zerolag_ctime0
            zerolag_results['foreground/%s/end_time' % self.ifos[1]] = zerolag_ctime1
            zerolag_results['foreground/%s/trigger_id1' % self.ifos[0]] = tid0
            zerolag_results['foreground/%s/trigger_id2' % self.ifos[1]] = tid1
            coinc_results.update(zerolag_results)

        if self.return_background:
            coinc_results['background/stat'] = self.coincs.data
            coinc_results['background/time'] = numpy.array([self.background_time])

            for ifo in self.singles:
                coinc_results['background/%s/count' % ifo] = numpy.array(self.singles[ifo].num_elements())            
        return coinc_results

