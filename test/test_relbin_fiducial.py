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
"""Tests that a fiducial waveform can be found rather than supplied.

Relative binning is only accurate near its fiducial waveform, and the
masses are what move the waveform, so these vary a mass and check the
answer against the same likelihood computed at full resolution, which
needs no fiducial at all.

A signal is injected into simulated noise so that there is a peak to
find, and the models are built directly, so the file runs in seconds.
"""

import copy
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
INJ = dict(mass1=1.42, mass2=1.36, distance=50., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)


class TestAutoFiducial(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        peak = float(hp.sample_times[numpy.argmax(abs(hp.data) ** 2
                                                  + abs(hc.data) ** 2)])
        cls.data, cls.psds = {}, {}
        seed = 5
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 97
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            signal.start_time += TC - peak
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass2=INJ['mass2'], f_lower=FLOW,
                          approximant='TaylorF2', tc=TC, ra=INJ['ra'],
                          dec=INJ['dec'], polarization=INJ['polarization'])
        cls.variable = ['distance', 'inclination', 'mass1']
        cls.prior = JointDistribution(
            cls.variable, SinAngle(inclination=None),
            Uniform(distance=(10, 200)), Uniform(mass1=(1.3, 1.55)))
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination'],
                     'mass1': INJ['mass1']}

        # the same likelihood without any relative binning at all
        exact = models.MarginalizedPhaseGaussianNoise(
            list(cls.variable), copy.deepcopy(cls.data),
            low_frequency_cutoff=cls.flow, psds=cls.psds,
            static_params=cls.static, prior=cls.prior)
        exact.update(**cls.point)
        cls.exact_loglr = exact.loglr

    def model(self, fiducial, optimize=False, epsilon=0.1):
        return models.Relative(
            list(self.variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior,
            fiducial_params=fiducial, epsilon=epsilon,
            optimize_fiducial=optimize)

    def test_a_bad_fiducial_really_does_hurt(self):
        """Establishes that there is something here worth fixing."""
        model = self.model({'mass1': 1.30})
        model.update(**self.point)
        self.assertGreater(abs(model.loglr - self.exact_loglr), 1.,
                           "a fiducial this far off should be inaccurate")

    def test_optimizing_recovers_the_exact_likelihood(self):
        """Starting from the same bad guess, it must find its way back."""
        model = self.model({'mass1': 1.30}, optimize=True, epsilon=0.03)
        model.update(**self.point)
        self.assertAlmostEqual(model.loglr, self.exact_loglr, delta=0.1)

    def test_the_fiducial_found_is_near_the_peak(self):
        """The point it settles on should be a signal, not noise."""
        model = self.model({'mass1': 1.30}, optimize=True)
        self.assertAlmostEqual(model.fid_params['mass1'], INJ['mass1'],
                               delta=0.02)

    def test_what_is_left_is_only_the_binning(self):
        """The fiducial found must leave a properly converging error.

        It does not land exactly on the test point, so some error remains,
        but it has to be the ordinary bin resolution error that goes away
        as bins are added, not a fixed offset from a fiducial in the wrong
        place.
        """
        errors = []
        for epsilon in (0.1, 0.03, 0.01):
            model = self.model({'mass1': 1.30}, optimize=True,
                               epsilon=epsilon)
            model.update(**self.point)
            errors.append(abs(model.loglr - self.exact_loglr))

        for coarse, fine in zip(errors, errors[1:], strict=False):
            self.assertGreater(coarse / fine, 3.,
                               "error should fall with the bins: %s" % errors)
        self.assertLess(errors[-1], 0.01, "%s" % errors)

    def test_only_shape_changing_parameters_are_searched(self):
        """Distance and inclination scale the waveform, they do not shape it.

        Searching over them walks along a degeneracy that says nothing
        about where the bins should go, and ends wherever the prior stops.
        """
        model = self.model({'mass1': INJ['mass1']})
        draws = model.prior_distribution.rvs(size=6)
        found = model.shape_params(list(self.variable), draws)
        self.assertEqual(found, ['mass1'])

    def test_the_fiducial_does_not_end_up_against_a_prior_bound(self):
        low, high = 10, 200
        model = self.model({'mass1': 1.30}, optimize=True)
        distance = model.fid_params['distance']
        self.assertGreater(distance - low, 1e-3 * (high - low))
        self.assertGreater(high - distance, 1e-3 * (high - low))

    def test_it_does_not_spoil_a_good_fiducial(self):
        """Optimizing when the fiducial was already fine must be safe."""
        given = self.model({'mass1': INJ['mass1']}, epsilon=0.03)
        given.update(**self.point)
        found = self.model({'mass1': INJ['mass1']}, optimize=True,
                           epsilon=0.03)
        found.update(**self.point)
        self.assertAlmostEqual(found.loglr, given.loglr, delta=0.1)

    def test_it_leaves_the_model_ready_to_use(self):
        """The search must not leave the model sitting somewhere odd."""
        model = self.model({'mass1': 1.30}, optimize=True)
        model.update(**self.point)
        self.assertEqual(model.current_params['mass1'], self.point['mass1'])


    def test_optimizing_without_a_prior_is_skipped(self):
        """There is nothing to draw a starting point from."""
        model = models.Relative(
            list(self.variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, fiducial_params={'mass1': 1.3756},
            epsilon=0.5, optimize_fiducial=True)
        model.update(**self.point)
        self.assertTrue(numpy.isfinite(model.loglr))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAutoFiducial))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
