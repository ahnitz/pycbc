import pycbc.types
import numpy

def fd_sine_gaussian(amp, quality, central_frequency, fmin, fmax, delta_f):
    f = numpy.arange(0, fmax, delta_f)
    tau = quality / 2 / numpy.pi / central_frequency
    d = amp * numpy.pi ** 0.5 / 2 * tau * numpy.exp(-numpy.pi ** 2.0 * tau * tau * (f - central_frequency) * (1 + numpy.exp(-quality ** 2.0 * f / central_frequency))
    return pycbc.types.FrequencySeries(d, delta_f = delta_f)
