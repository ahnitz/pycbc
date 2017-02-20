""" Utitlities for introducing tidal effects into waveform approximants
"""
import pycbc.pnutils
import numpy

def nonlinear_phase_difference(f, f0, A, n, m1, m2):
    """ Implmenents the phase difference approximation of nonlinear
    tides in Essick, et al. https://arxiv.org/pdf/1609.06362v2.pdf
    """
    x0 = f0 / 100.0
    x = f / 100.0
    mc, et = pycbc.pnutils.mass1_mass2_to_mchirp_eta(m1, m2)
    dphi = 0.4 * (mc / 1.2) ** (-10.0 / 3) * A * 10**8
    dphi *= (x0**(n - 3.0) - x**(n - 3.0)) / (n - 3.0)
    return dphi

def nonlinear_tidal_spa(**kwds):
    from . import spa_tmplt
    from pycbc.types import zeros, Array

    # We start with the standard TaylorF2 based waveform
    tmplt = spa_tmplt.spa_tmplt(**kwds)

    # Add the phasing difference from the nonlinear tides
    kmin = int((kwds['f0'] / tmplt.delta_f))
    f = numpy.arange(kmin, len(tmplt)) * tmplt.delta_f
    tmplt[kmin:] *= Array(numpy.exp(1.0j * nonlinear_phase_difference(f,
               kwds['f0'], kwds['A'], kwds['n'], kwds['mass1'], kwds['mass2'])), dtype=tmplt.dtype)
    return tmplt

def nonlinear_tidal_spa_full(**kwds):
    v = nonlinear_tidal_spa(**kwds)
    return v, v * 1.0j
