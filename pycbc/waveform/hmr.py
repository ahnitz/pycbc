import h5py
from pycbc.types import load_timeseries, TimeSeries
from pycbc.conversions import mchirp_from_mass1_mass2, mass1_from_mchirp_q, mass2_from_mchirp_q

import numpy
from scipy.interpolate import interp1d
import os


eobfile = os.environ['HMR_FILE']
if 'HMR_FILE_SHARED' in os.environ:
    eobfile_src = eobfile
    eobfile = '/dev/shm/HMR_FILE.hdf'

    from filelock import FileLock
    import os.path
    import shutil
    with FileLock(eobfile + '.lock'):
        recopy = False
        if os.path.exists(eobfile):
            try:
                f = h5py.File(eobfile, 'r')
                test = f['0']
            except:
                recopy = True

        if not os.path.exists(eobfile) or recopy:
            shutil.copy(eobfile_src, eobfile)

f = h5py.File(eobfile, 'r')
m1s = numpy.array([f[k].attrs['m1'] for k in f])
m2s = numpy.array([f[k].attrs['m2'] for k in f])
qs = m1s / m2s
fl = numpy.array([f[k].attrs['flow'] for k in f])
keys = [k for k in f]
f.close()

def getfeob(**kwds):
    if 'duration' in kwds:
        duration = kwds['duration']
    else:
        duration = 100.0
    hp = feob(kwds['mass1'], kwds['mass2'],
                kwds['delta_f'],
                duration=duration).astype(numpy.complex128)
    kmin = int(kwds['f_lower'] / float(kwds['delta_f']))
    hp[:kmin].clear()
    hp /= kwds['distance']
    return hp, hp*1.0j

def ieob(m1, m2, delta_t, duration=100.0):
    # get closest mass ratio in set
    q = m1 / m2
    j = abs(qs - q).argmin()

    ts = load_timeseries(eobfile, group=keys[j])
    M_ref = mchirp_from_mass1_mass2(m1s[j], m2s[j])
    M_targ = mchirp_from_mass1_mass2(m1, m2)
    Mr = (M_targ / M_ref)

#    sx = ts.sample_times * Mr
    st = ts.sample_times
    inter = interp1d(st, ts, fill_value=0.0, kind='cubic', bounds_error=False)

    x = numpy.arange(-duration / Mr, st[-1]+delta_t, delta_t/Mr)
    y = inter(x)
    return TimeSeries(y, epoch=x[0], delta_t=delta_t)

def feob(m1, m2, delta_f, duration=100.0):
    # get closest mass ratio in set
    q = max(float(m1) / float(m2), float(m2) / float(m1))
    j = abs(qs - q).argmin()
    ts = load_timeseries(eobfile, group=keys[j])
    M_ref = mchirp_from_mass1_mass2(m1s[j], m2s[j])
    M_targ = mchirp_from_mass1_mass2(m1, m2)
    Mr = (M_targ / M_ref)
    st = ts.sample_times

    #print(ts.duration, duration, Mr)
    if ts.duration > duration / Mr:
        ts = ts.time_slice(st[-1] - duration/Mr, st[-1])

    sr = ts.sample_rate / Mr
    bf = 1.0 / delta_f

    tlen = sr * bf
    if len(ts) >= tlen:
        ts = ts[len(ts) - tlen:]
    else:
        ts.resize(sr * bf)

    from pycbc.waveform.utils import taper_timeseries
    ts = taper_timeseries(ts, tapermethod='start')

    fs = ts.to_frequencyseries().cyclic_time_shift(ts.start_time) * Mr**2.0
    fs._epoch = 0
    fs._delta_f = delta_f
    return fs
