import pycbc.types
import numpy

def fd_sine_gaussian(amp, quality, central_frequency, fmin, fmax, delta_f):
    f = numpy.arange(fmin, fmax, delta_f)
    kmax = int(fmax / delta_f)
    kmin = int(fmin / delta_f)
    tau = quality / 2 / numpy.pi / central_frequency
    d = amp * numpy.pi ** 0.5 / 2 * tau 
    d *= numpy.exp(-(numpy.pi  * tau  * (f - central_frequency))**2.0)
    d *= (1 + numpy.exp(-quality ** 2.0 * f / central_frequency))
    v = numpy.zeros(kmax, dtype=numpy.complex128)
    v[kmin:kmax] = d[:]
    return pycbc.types.FrequencySeries(v, delta_f=delta_f, copy=False)
