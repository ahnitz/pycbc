# Copyright (C) 2026 Alexander Nitz
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

"""Tests of the multi-signal likelihood of the relative binning model.

When several signals overlap in the same data the likelihood ratio picks up a
cross term for every pair, -Re<h_i|h_j>, which is what the relative binning
model calculates in ``multi_loglikelihood``. That term is zero-error at the
point the summary data was built at, so it can be badly wrong everywhere else
and still look right in a spot check; these tests therefore compare against a
direct, unbinned inner product at parameters displaced from the fiducial ones.
"""

import unittest
import numpy

from pycbc.detector import Detector
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import FrequencySeries
from pycbc.waveform import get_fd_waveform_sequence
from pycbc.inference.models import Relative
from utils import parse_args_cpu_only, simple_exit

parse_args_cpu_only("Relative binning multi-signal model")

IFOS = ['H1', 'L1']
SEGLEN = 8
DF = 1.0 / SEGLEN
FLOW = 25.0
FHIGH = 300.0
FNYQ = 512.0
TC = 1234567890.0

# Three overlapping signals. They differ in sky position and time as well as
# in mass, so the pieces the cross term has to line up -- the reference
# waveform of one model against the time shift of another -- are all distinct.
SIGNALS = [
    dict(mass1=40., mass2=35., spin1z=0.1, spin2z=0.05, distance=500.,
         inclination=0.4, coa_phase=1.1, polarization=0.3,
         ra=1.7, dec=-0.5, tc=TC),
    dict(mass1=22., mass2=18., spin1z=-0.2, spin2z=0.3, distance=400.,
         inclination=1.1, coa_phase=0.2, polarization=2.0,
         ra=0.6, dec=0.9, tc=TC + 0.01),
    dict(mass1=30., mass2=12., spin1z=0.4, spin2z=-0.1, distance=450.,
         inclination=0.8, coa_phase=2.4, polarization=1.2,
         ra=4.1, dec=0.2, tc=TC - 0.02),
]
FIXED = dict(f_lower=20., approximant='IMRPhenomD')
# The same signals with their higher modes. The modes beat against each
# other across the band, so the ratio of signal to reference carries
# structure the dominant mode alone does not, which is what the cross term
# has to interpolate onto a grid neither model owns.
HIGHER_MODES = dict(f_lower=20., approximant='IMRPhenomXHM')
VARIABLE = ['distance', 'coa_phase', 'tc']

# Displacements applied to every signal in turn: a shift in coalescence
# phase, a factor on the distance and a shift in coalescence time. The
# first is the fiducial point itself, where any cross term agrees. The time
# shift is what makes the ratio of signal to reference vary across the band.
OFFSETS = [(0.0, 1.0, 0.0), (0.2, 1.05, 0.001),
           (1.0, 0.8, -0.002), (2.5, 1.3, 0.003)]

# How far the binned likelihood may sit from the unbinned one, in nats. This
# is the truncation error of the binning at the epsilon used below, measured
# at 0.0057 nats for one signal and 0.0071 for three; the defects it guards
# against were 138 and 600 nats at these same displacements.
TOL = 0.02


def project(params, ifo, freqs, end_time, fixed=FIXED):
    """The detector frame waveform on the given frequencies."""
    hp, hc = get_fd_waveform_sequence(sample_points=freqs,
                                      **dict(params, **fixed))
    det = Detector(ifo)
    fp, fc = det.antenna_pattern(params['ra'], params['dec'],
                                 params['polarization'], params['tc'])
    dt = det.time_delay_from_earth_center(params['ra'], params['dec'],
                                          params['tc'])
    h = fp * hp.numpy() + fc * hc.numpy()
    return h * numpy.exp(-2.0j * numpy.pi * numpy.array(freqs)
                         * (params['tc'] + dt - end_time))


class TestRelbinMultiSignal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(FNYQ / DF) + 1
        cls.psd = aLIGOZeroDetHighPower(flen, DF, 15.0)
        cls.psds = {ifo: cls.psd for ifo in IFOS}
        cls.data = cls.make_data(SIGNALS, FIXED)

    @classmethod
    def make_data(cls, signals, fixed):
        """Zero noise data holding all of the signals at once."""
        flen = int(FNYQ / DF) + 1
        data = {}
        for ifo in IFOS:
            d = FrequencySeries(numpy.zeros(flen, dtype=numpy.complex128),
                                delta_f=DF, epoch=TC - SEGLEN + 2)
            freqs = d.sample_frequencies.numpy()
            sl = slice(int(FLOW / DF), int(FHIGH / DF) + 1)
            for params in signals:
                d.data[sl] += project(params, ifo, freqs[sl],
                                      float(d.end_time), fixed)
            data[ifo] = d
        return data

    def make_model(self, params, epsilon=0.05, flow=FLOW, fhigh=FHIGH,
                   data=None, fixed=FIXED):
        static = {k: v for k, v in params.items() if k not in VARIABLE}
        static.update(fixed)
        data = self.data if data is None else data
        return Relative(VARIABLE,
                        {ifo: data[ifo].copy() for ifo in IFOS},
                        low_frequency_cutoff={ifo: flow for ifo in IFOS},
                        high_frequency_cutoff={ifo: fhigh for ifo in IFOS},
                        psds=self.psds, static_params=static,
                        fiducial_params=dict(params, **fixed),
                        epsilon=epsilon, marginalize_phase=False)

    def direct_loglr(self, models, plist, data=None, fixed=FIXED):
        """<d|h> - 0.5<h|h> for the summed signal, without any binning.

        The sum runs over exactly the samples the binned models cover, so
        that the only difference left is the binning itself.
        """
        total = 0.0
        data = self.data if data is None else data
        for ifo in IFOS:
            d = data[ifo]
            freqs = d.sample_frequencies.numpy()
            sl = slice(max(m.edges[ifo][0] for m in models),
                       min(m.edges[ifo][-1] for m in models))
            h = sum(project(p, ifo, freqs[sl], float(d.end_time), fixed)
                    for p in plist)
            psd = self.psd.numpy()[sl]
            sh = 4 * DF * (numpy.conjugate(h) * d.numpy()[sl] / psd).sum()
            hh = 4 * DF * ((numpy.abs(h) ** 2) / psd).sum()
            total += sh.real - 0.5 * hh.real
        return total

    def displaced(self, index, offset, scale, dt):
        params = dict(SIGNALS[index])
        params['coa_phase'] = params['coa_phase'] + offset
        params['distance'] = params['distance'] * scale
        params['tc'] = params['tc'] + dt
        return params

    def check_models(self, models, tolerance, data=None, fixed=FIXED):
        """Multi-signal likelihood versus the summed waveform likelihood."""
        errors = []
        for offset, scale, dt in OFFSETS:
            plist = [self.displaced(i, offset * (i + 1), scale, dt / (i + 1))
                     for i in range(len(models))]
            for model, params in zip(models, plist):
                model.update(distance=params['distance'],
                             coa_phase=params['coa_phase'], tc=params['tc'])
            loglr = models[0].multi_loglikelihood(models[1:]) \
                - models[0].lognl
            errors.append(loglr - self.direct_loglr(models, plist,
                                                    data, fixed))
        largest = numpy.abs(errors).max()
        self.assertLess(largest, tolerance,
                        "multi-signal loglr is off by up to {} nats; "
                        "errors at each point: {}".format(largest, errors))
        return errors

    def test_single_signal_reference(self):
        """The unbinned reference has to agree with one model on its own.

        Without this the other tests could be comparing two different
        conventions for the waveform rather than testing the cross term.
        """
        model = self.make_model(SIGNALS[0])
        for offset, scale, dt in OFFSETS:
            params = self.displaced(0, offset, scale, dt)
            model.update(distance=params['distance'],
                         coa_phase=params['coa_phase'], tc=params['tc'])
            self.assertAlmostEqual(model.loglr,
                                   self.direct_loglr([model], [params]),
                                   delta=TOL)

    def test_two_signals(self):
        """Two overlapping signals sharing a bin layout."""
        models = [self.make_model(p) for p in SIGNALS[:2]]
        self.check_models(models, TOL)

    def test_three_signals(self):
        """All three pairs of cross terms have to be right at once."""
        models = [self.make_model(p) for p in SIGNALS]
        self.check_models(models, TOL)

    def band_loglr(self, models, plist):
        """The reference when the models do not analyze the same band.

        A signal contributes where its own model is looking, and a pair of
        signals interfere only where both are, so the two terms are summed
        over different ranges. With one band this reduces to direct_loglr.
        """
        total = 0.0
        for ifo in IFOS:
            d = self.data[ifo]
            freqs = d.sample_frequencies.numpy()
            psd = self.psd.numpy()
            waves, spans = [], []
            for model, params in zip(models, plist):
                lo, hi = model.edges[ifo][0], model.edges[ifo][-1]
                spans.append((lo, hi))
                waves.append(project(params, ifo, freqs[lo:hi],
                                     float(d.end_time)))
            for (lo, hi), h in zip(spans, waves):
                sh = 4 * DF * (numpy.conjugate(h) * d.numpy()[lo:hi]
                               / psd[lo:hi]).sum()
                hh = 4 * DF * ((numpy.abs(h) ** 2) / psd[lo:hi]).sum()
                total += sh.real - 0.5 * hh.real
            for i in range(len(models)):
                for j in range(i + 1, len(models)):
                    lo = max(spans[i][0], spans[j][0])
                    hi = min(spans[i][1], spans[j][1])
                    if hi <= lo:
                        continue
                    hi_i = waves[i][lo - spans[i][0]:hi - spans[i][0]]
                    hj = waves[j][lo - spans[j][0]:hi - spans[j][0]]
                    cross = 4 * DF * (numpy.conjugate(hi_i) * hj
                                      / psd[lo:hi]).sum()
                    total -= cross.real
        return total

    def test_higher_modes(self):
        """Signals carrying their higher modes, not the dominant one alone.

        The modes beat against each other across the band, so the ratio of
        signal to reference is far less linear than for the dominant mode
        alone, and the cross term has to interpolate that ratio onto a grid
        neither model owns. That makes the cross term the demanding part:
        at epsilon 0.05 a single higher-mode signal is good to 0.015 nats,
        while the pair is off by 1.7. Both reach 0.0008 at 0.01, so a
        multi-signal analysis of waveforms with higher modes wants finer
        bins than the same waveforms one at a time.
        """
        data = self.make_data(SIGNALS[:2], HIGHER_MODES)
        models = [self.make_model(p, epsilon=0.01, data=data,
                                  fixed=HIGHER_MODES)
                  for p in SIGNALS[:2]]
        self.check_models(models, TOL, data, HIGHER_MODES)

    def test_differing_bands(self):
        """Submodels do not have to analyze the same frequencies.

        Where only one model is looking, only that signal contributes; the
        two interfere only where both are. The reference waveform of a
        model is zero outside its own band, so an edge taken from the
        other model's range would be divided by zero if it were kept.
        """
        models = [self.make_model(SIGNALS[0]),
                  self.make_model(SIGNALS[1], flow=FLOW + 15.,
                                  fhigh=FHIGH - 60.)]
        for ifo in IFOS:
            self.assertGreater(models[1].edges[ifo][0],
                               models[0].edges[ifo][0])
            self.assertLess(models[1].edges[ifo][-1],
                            models[0].edges[ifo][-1])

        errors = []
        for offset, scale, dt in OFFSETS:
            plist = [self.displaced(i, offset * (i + 1), scale, dt / (i + 1))
                     for i in range(len(models))]
            for model, params in zip(models, plist):
                model.update(distance=params['distance'],
                             coa_phase=params['coa_phase'], tc=params['tc'])
            loglr = models[0].multi_loglikelihood(models[1:]) \
                - models[0].lognl
            errors.append(loglr - self.band_loglr(models, plist))
        largest = numpy.abs(errors).max()
        self.assertTrue(numpy.isfinite(largest),
                        "the cross term is not finite: {}".format(errors))
        self.assertLess(largest, TOL,
                        "multi-signal loglr is off by up to {} nats; "
                        "errors at each point: {}".format(largest, errors))

    def test_differing_bin_layouts(self):
        """Submodels do not have to agree on where the bin edges are.

        Any difference in the band, the reference waveform or the binning
        tolerance leaves the two models with unrelated sets of bin edges,
        and the cross term then has to be evaluated on a grid that is
        native to neither of them. Different epsilons is the cheapest way
        to force that here.
        """
        models = [self.make_model(p, epsilon=e)
                  for p, e in zip(SIGNALS[:2], [0.05, 0.02])]
        for ifo in IFOS:
            self.assertNotEqual(len(models[0].edges[ifo]),
                                len(models[1].edges[ifo]))
        self.check_models(models, TOL)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestRelbinMultiSignal))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
