import pycbc.types
import numpy

def fd_sine_gaussian(amp, quality, central_frequency, fmin, fmax, delta_f):
    f = numpy.arange(0, fmax, delta_f)
    tau = quality / 2 / numpy.pi / central_frequency
    d = amp * numpy.pi ** 0.5 / 2 * tau 
    d *= numpy.exp(-(numpy.pi  * tau  * (f - central_frequency))**2.0)
    d *= (1 + numpy.exp(-quality ** 2.0 * f / central_frequency))
    d[0:int(fmin / delta_f)] = 0
    return pycbc.types.FrequencySeries(d.astype(numpy.complex128), delta_f = delta_f)

