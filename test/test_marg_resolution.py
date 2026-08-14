# Copyright (C) 2026  Alex Nitz
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
"""Tests choosing the time resolution of the marginalization automatically.

sample_rate sets the spacing of the signal to noise series the times are
drawn from, and what it costs is the scatter of the marginalized
likelihood between one evaluation of a point and the next. The rate is
chosen by asking for that scatter to be no larger than an accuracy, so
what these check is the accuracy: that the chosen rate delivers the
scatter it promised, that asking for less scatter buys a finer rate, and
that a coarser rate would not have done.

The signal here is deliberately not the one the rule was measured on --
it is a neutron star black hole binary in a different noise realisation,
where the peak of the likelihood has a different width -- so passing
means the rule carried over rather than that a number was fitted.
"""

import copy
import unittest

import numpy
from scipy.special import logsumexp
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 30., 16, 2048
HALFWIDTH = 0.004
INJ = dict(mass1=3.0, mass2=1.4, distance=130., inclination=0.9,
           ra=0.3, dec=1.1, polarization=2.2, coa_phase=0.4)


class TestMargResolution(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='IMRPhenomD', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        # the generator puts coalescence at t=0, so this puts it at tc
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        seed = 91
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 53
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='IMRPhenomD',
                          ra=INJ['ra'], dec=INJ['dec'],
                          polarization=INJ['polarization'])
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination']}

    VARIABLE = ['distance', 'inclination', 'tc']

    def prior(self):
        return JointDistribution(
            list(self.VARIABLE), SinAngle(inclination=None),
            Uniform(distance=(10, 300)),
            Uniform(tc=(TC - HALFWIDTH, TC + HALFWIDTH)))

    def model(self, sample_rate='auto', accuracy=0.005, vsamples=2000):
        return models.RelativeTimeDom(
            list(self.VARIABLE), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior(),
            fiducial_params={'mass1': INJ['mass1'], 'mass2': INJ['mass2'],
                             'tc': TC},
            epsilon=0.1, marginalize_vector_params='tc',
            sample_rate=sample_rate, marginalization_accuracy=accuracy,
            marginalize_vector_samples=vsamples)

    def repeated(self, model, ndraw=64, seed=4):
        """The same point evaluated over and over, which is what the times
        being drawn afresh each call makes scatter."""
        state = numpy.random.get_state()
        try:
            numpy.random.seed(seed)
            values = []
            for _ in range(ndraw):
                model.update(**self.point)
                values.append(float(model.loglr))
        finally:
            numpy.random.set_state(state)
        return values

    def scatter(self, model, **kwargs):
        return numpy.std(self.repeated(model, **kwargs), ddof=1)

    def integrated(self, npoint=2001):
        """The marginal computed by summing the likelihood over the prior,
        which needs no marginalization machinery to be correct."""
        model = models.Relative(
            list(self.VARIABLE), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior(),
            fiducial_params={'mass1': INJ['mass1'], 'mass2': INJ['mass2'],
                             'tc': TC}, epsilon=0.1)
        values = []
        for tc in numpy.linspace(TC - HALFWIDTH, TC + HALFWIDTH, npoint):
            model.update(tc=tc, **self.point)
            values.append(model.loglr)
        # a uniform prior makes the integral the mean over the window
        return logsumexp(values) - numpy.log(npoint)

    def test_the_chosen_rate_gets_the_right_answer(self):
        """Being quiet is not the same as being right.

        The scatter the rate is chosen from says how far apart two
        evaluations land, not where they land, so the answer is also
        checked against direct integration of the unmarginalized model.
        """
        values = self.repeated(self.model(), ndraw=16, seed=11)
        error = numpy.std(values, ddof=1) / len(values) ** 0.5
        self.assertAlmostEqual(numpy.mean(values), self.integrated(),
                               delta=max(0.02, 4 * error))

    def test_the_chosen_rate_delivers_the_accuracy(self):
        """The whole point: the rate is chosen from a promise about the
        scatter, so the scatter has to be the promised size.

        The tolerance is loose by half because the rule is a fit good to
        about 20%, and 64 draws pin a standard deviation to about 10%.
        """
        accuracy = 0.005
        model = self.model(accuracy=accuracy)
        self.assertLess(self.scatter(model), 1.5 * accuracy,
                        "chose %s Hz" % model.sample_rate)

    def test_the_knobs_point_the_right_way(self):
        """Resolution and the number of draws buy the same thing.

        The scatter goes as one over the root of their product, so asking
        for less of it must buy a finer rate, and asking for the same with
        four times the draws must cost a coarser series rather than the
        same one.
        """
        rates = [self.model(accuracy=a).sample_rate
                 for a in (0.02, 0.005, 0.002)]
        self.assertTrue(rates[0] < rates[1] < rates[2], "%s" % rates)

        few = self.model(vsamples=1000).sample_rate
        many = self.model(vsamples=4000).sample_rate
        self.assertLess(many, few, "%s against %s" % (many, few))

    def test_a_coarser_rate_would_not_have_done(self):
        """The series is rebuilt at every likelihood evaluation, so a rate
        beyond what the accuracy calls for is a running cost, and one below
        it really is noisier.

        The first check is on the prediction, which is the same statement
        one step earlier and costs no likelihood calls; the second measures
        the scatter itself, on a signal the rule was not measured on,
        rather than trusting that it tracks. Cutting the rate by eight cuts
        the samples across the peak by eight, so the rule says the scatter
        grows by getting on for three. It is compared loosely because the
        finer scatter is near the floor that no rate removes, which makes
        the ratio smaller than the rule alone says.
        """
        accuracy = 0.005
        model = self.model(accuracy=accuracy)
        coarser = self.model(sample_rate=model.sample_rate / 2)
        predicted = coarser.marginalization_error(
            coarser.resolved_samples(coarser.ref_snr))
        self.assertGreater(predicted, accuracy,
                           "%s Hz would have done" % coarser.sample_rate)

        coarse = self.model(sample_rate=model.sample_rate / 8)
        self.assertGreater(self.scatter(coarse), 1.5 * self.scatter(model))

    def test_an_explicit_rate_is_left_alone(self):
        """Asking for a rate must still get exactly that rate."""
        self.assertEqual(self.model(sample_rate=8192).sample_rate, 8192.0)

    def test_an_unresolved_peak_is_doubled_up_to(self):
        """With nothing to measure there is nothing to solve for.

        A series that leaves the peak on about one sample cannot say how
        much finer it needs to be, so the rate doubles until it can. Asked
        for an accuracy so loose that only the two sample floor binds, the
        rate must still climb off the rate of the data to reach it.
        """
        model = self.model(accuracy=1.0)
        self.assertGreaterEqual(model.resolved_samples(model.ref_snr), 2.0)
        self.assertGreater(model.sample_rate, SRATE)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargResolution))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
