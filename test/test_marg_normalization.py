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
"""Tests that time marginalization integrates rather than averages.

The marginalized likelihood is an integral of the likelihood against the
prior, so it has to depend on the prior it was integrated against and not
on the grid it happened to be computed on. Both are checked here against
direct integration of the unmarginalized model, which needs no
marginalization machinery to be correct.
"""

import copy
import unittest

import numpy
from scipy.special import logsumexp
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
)
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
INJ = dict(mass1=1.42, mass2=1.36, distance=50., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)


class TestMargNormalization(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        # the generator puts coalescence at t=0, so this puts it at tc
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        seed = 5
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 97
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='TaylorF2',
                          ra=INJ['ra'], dec=INJ['dec'],
                          polarization=INJ['polarization'])
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination']}

    def prior(self, halfwidth):
        variable = ['distance', 'inclination', 'tc']
        return variable, JointDistribution(
            list(variable), SinAngle(inclination=None),
            Uniform(distance=(10, 200)),
            Uniform(tc=(TC - halfwidth, TC + halfwidth)))

    def sky_marginalized(self, npoint, nseed=6):
        """Mean and spread of the sky and time marginalized likelihood."""
        variable = ['distance', 'inclination', 'tc', 'polarization',
                    'ra', 'dec']
        static = {k: v for k, v in self.static.items()
                  if k not in ('ra', 'dec', 'polarization')}
        dists = [SinAngle(inclination=None), Uniform(distance=(10, 200)),
                 Uniform(tc=(TC - 0.1, TC + 0.1)),
                 UniformAngle(polarization=None), UniformAngle(ra=None),
                 CosAngle(dec=None)]
        values = []
        for s in range(nseed):
            numpy.random.seed(700 + s)
            model = models.MarginalizedTime(
                list(variable), copy.deepcopy(self.data),
                low_frequency_cutoff=self.flow, psds=self.psds,
                static_params=static,
                prior=JointDistribution(list(variable), *dists),
                marginalize_phase=True,
                marginalize_vector_params='tc,ra,dec,polarization',
                marginalize_vector_samples=npoint, sample_rate=4096,
                marginalize_sky_initial_samples=2e5)
            model.update(**self.point)
            values.append(model.loglr)
        return numpy.mean(values), numpy.std(values) / nseed ** 0.5

    def test_sky_marginalization_does_not_depend_on_the_point_count(self):
        """The number of points must buy precision, not a bigger answer.

        The times are drawn from each detector's signal to noise series,
        so what they carry is their probability under that series. Taking
        it over the drawn times rather than over the series made it one in
        however many were drawn, and the answer grew by the log of that:
        over a range of sixty four in the count it climbed by four and a
        half, far outside its own scatter.
        """
        seen = [(n,) + self.sky_marginalized(n) for n in (128, 2048)]
        (_, coarse, coarse_err), (_, fine, fine_err) = seen
        self.assertLess(abs(coarse - fine),
                        4 * (coarse_err ** 2 + fine_err ** 2) ** 0.5,
                        "%s" % seen)

    def test_more_sky_points_are_more_precise(self):
        """What the count does buy is precision."""
        _, coarse = self.sky_marginalized(128)
        _, fine = self.sky_marginalized(2048)
        self.assertLess(fine, coarse,
                        "error on the mean went %.4f to %.4f" % (coarse, fine))

    def integrated(self, halfwidth, npoint=4001):
        """The marginal computed by summing the likelihood over the prior."""
        variable, prior = self.prior(halfwidth)
        model = models.Relative(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=prior,
            fiducial_params={'mass1': INJ['mass1'], 'tc': TC}, epsilon=0.1)
        values = []
        for tc in numpy.linspace(TC - halfwidth, TC + halfwidth, npoint):
            model.update(tc=tc, **self.point)
            values.append(model.loglr)
        # a uniform prior makes the integral the mean over the window
        return logsumexp(values) - numpy.log(npoint)

    def marginalized(self, halfwidth, sample_rate=8192, nseed=8):
        variable, prior = self.prior(halfwidth)
        values = []
        for s in range(nseed):
            numpy.random.seed(20 + s)
            model = models.RelativeTime(
                list(variable), copy.deepcopy(self.data),
                low_frequency_cutoff=self.flow, psds=self.psds,
                static_params=self.static, prior=prior,
                fiducial_params={'mass1': INJ['mass1'], 'tc': TC},
                epsilon=0.1, marginalize_vector_params='tc',
                sample_rate=sample_rate, marginalize_vector_samples=2000)
            model.update(**self.point)
            values.append(model.loglr)
        # the error on the mean, which is what the tests compare against
        return numpy.mean(values), numpy.std(values) / len(values) ** 0.5

    def test_matches_direct_integration(self):
        """The whole point: it must be the integral it claims to be."""
        for halfwidth in (0.001, 0.004):
            mean, error = self.marginalized(halfwidth)
            self.assertAlmostEqual(
                mean, self.integrated(halfwidth), delta=4 * error,
                msg="half-width %s, error %.3f" % (halfwidth, error))

    def test_widening_the_prior_lowers_the_evidence(self):
        """Doubling the window must cost log(2), not nothing.

        Averaging instead of integrating gives an answer that barely moves,
        which is what this is here to catch.
        """
        narrow, _ = self.marginalized(0.002)
        wide, _ = self.marginalized(0.008)
        self.assertAlmostEqual(narrow - wide, 2 * numpy.log(2), delta=0.15)

    def test_the_answer_settles_as_the_time_grid_is_refined(self):
        """sample_rate must be a resolution knob and nothing more.

        It sets the spacing of the time series the times are drawn from.
        The peak of that series is well under a millisecond wide for a
        signal like this one, so a coarse grid does not resolve it and the
        answer is off; what matters is that refining the grid settles,
        which is what lets the rate be chosen to meet an accuracy. Before
        the weights were fixed the answer drifted instead, and no choice of
        rate was any more right than another.
        """
        values = [self.marginalized(0.004, sample_rate=rate)[0]
                  for rate in (4096, 8192, 16384)]
        steps = numpy.abs(numpy.diff(values))
        self.assertLess(steps[-1], steps[0] / 3.,
                        "refining the grid should settle: %s" % values)
        self.assertLess(steps[-1], 0.05, "%s" % values)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestMargNormalization))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
