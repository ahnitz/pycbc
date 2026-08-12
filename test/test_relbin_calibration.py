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
"""Tests calibration error support in the relative binning model.

The correction is a spline through a few control points, so it varies far
more slowly with frequency than the waveform does. That is what makes it
affordable here: it can be evaluated at the bin edges alone. These check
that doing so gives the same answer as correcting the waveform at every
frequency of the data, which is what the full resolution model does.
"""

import copy
import time
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.strain.recalibrate import CubicSpline
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
INJ = dict(mass1=1.42, mass2=1.36, distance=50., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)
NPOINT = 4


class TestRelbinCalibration(unittest.TestCase):

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
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='TaylorF2', tc=TC,
                          ra=INJ['ra'], dec=INJ['dec'],
                          polarization=INJ['polarization'])
        cls.variable = ['distance', 'inclination']
        cls.prior = JointDistribution(
            cls.variable, SinAngle(inclination=None),
            Uniform(distance=(10, 200)))
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination']}

    def recalibration(self):
        return {ifo: CubicSpline(minimum_frequency=FLOW,
                                 maximum_frequency=1000.,
                                 n_points=NPOINT, ifo_name=ifo)
                for ifo in self.data}

    def calibration_params(self, amplitude=0.05, phase=0.05):
        """A correction big enough to matter, and not flat with frequency."""
        params = {}
        for ifo in self.data:
            for i in range(NPOINT):
                sign = 1 if i % 2 == 0 else -1
                params['recalib_amplitude_%s_%d' % (ifo, i)] = sign * amplitude
                params['recalib_phase_%s_%d' % (ifo, i)] = sign * phase
        return params

    def relbin(self, recalibration=None, epsilon=0.03):
        return models.Relative(
            list(self.variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior,
            fiducial_params={}, epsilon=epsilon,
            recalibration=recalibration)

    def exact(self, recalibration=None):
        return models.MarginalizedPhaseGaussianNoise(
            list(self.variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=self.prior,
            recalibration=recalibration)

    def test_calibration_actually_changes_the_answer(self):
        """Otherwise the rest of these tests would prove nothing."""
        without = self.exact()
        without.update(**self.point)

        with_cal = self.exact(self.recalibration())
        with_cal.update(**self.point, **self.calibration_params())

        self.assertGreater(abs(with_cal.loglr - without.loglr), 1.)

    def test_relbin_matches_full_resolution_with_calibration(self):
        """Correcting at the bin edges must match correcting everywhere."""
        params = dict(self.point, **self.calibration_params())

        exact = self.exact(self.recalibration())
        exact.update(**params)

        model = self.relbin(self.recalibration())
        model.update(**params)

        self.assertAlmostEqual(model.loglr, exact.loglr, delta=0.05)

    def test_it_holds_for_a_larger_correction(self):
        """A correction well beyond what is expected must still agree."""
        params = dict(self.point,
                      **self.calibration_params(amplitude=0.2, phase=0.2))

        exact = self.exact(self.recalibration())
        exact.update(**params)

        model = self.relbin(self.recalibration())
        model.update(**params)

        self.assertAlmostEqual(model.loglr, exact.loglr, delta=0.05)

    def test_no_calibration_is_unchanged(self):
        """Models without calibration must give exactly what they did."""
        plain = self.relbin()
        plain.update(**self.point)

        model = self.relbin(self.recalibration())
        model.update(**self.point, **self.calibration_params(0., 0.))

        self.assertAlmostEqual(model.loglr, plain.loglr, delta=1e-6)

    def test_calibration_is_nearly_free(self):
        """The whole point is that it costs bins, not data samples."""
        params = dict(self.point, **self.calibration_params())

        plain = self.relbin()
        model = self.relbin(self.recalibration())

        def timed(m, point):
            m.update(**point)
            m.loglr
            start = time.time()
            for i in range(100):
                m.update(**dict(point, distance=40. + 0.01 * i))
                m.loglr
            return (time.time() - start) / 100

        cost = timed(model, params) / timed(plain, self.point)
        self.assertLess(cost, 2.0,
                        "calibration made the likelihood %.1f times slower"
                        % cost)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinCalibration))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
