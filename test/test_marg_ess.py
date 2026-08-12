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
"""Tests the effective sample size the vector marginalization reports.

The marginalized likelihood is a weighted average over drawn points. The
effective sample size is how many of them the answer really rests on, and
the error of the answer is about one over its square root. This checks
that the reported number means that: that it rises with the points drawn,
and that one over its root really is the spread of the answer.
"""

import copy
import unittest

import numpy

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

FLOW, SEGLEN, SRATE, TC = 25., 32, 2048, 1187008882.42840


class TestMargESS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # a fixed binary neutron star signal injected into simulated noise;
        # self-contained so the test does not depend on a shared fixture
        inj = dict(mass1=1.4, mass2=1.35, distance=60., inclination=0.5,
                   ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: inj[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        hp.start_time += TC
        hc.start_time += TC
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        cls.data, cls.psds = {}, {}
        seed = 3
        for ifo in ['H1', 'L1']:
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=seed)
            noise._epoch = TC - SEGLEN / 2
            seed += 101
            signal = Detector(ifo).project_wave(
                hp, hc, inj['ra'], inj['dec'], inj['polarization'])
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=inj['mass1'], mass2=inj['mass2'], f_lower=FLOW,
                          approximant='TaylorF2', ra=inj['ra'], dec=inj['dec'],
                          polarization=inj['polarization'])
        cls.point = {'distance': inj['distance'],
                     'inclination': inj['inclination']}

    def build(self, npoint):
        # time marginalization draws a genuine vector of times from the
        # signal-to-noise series, which is what an effective sample size
        # is about; polarization is done analytically, with no vector
        variable = ['distance', 'inclination', 'tc']
        return models.MarginalizedTime(
            variable, copy.deepcopy(self.data), low_frequency_cutoff=self.flow,
            psds=self.psds, static_params=self.static,
            prior=JointDistribution(variable, SinAngle(inclination=None),
                                    Uniform(distance=(10, 200)),
                                    Uniform(tc=(TC - 0.1, TC + 0.1))),
            marginalize_phase=True, marginalize_vector_params='tc',
            sample_rate=4096, marginalize_vector_samples=npoint)

    def ess_and_spread(self, npoint, nseed=12):
        ess, loglr = [], []
        for s in range(nseed):
            numpy.random.seed(500 + s)
            model = self.build(npoint)
            model.update(**self.point)
            loglr.append(model.loglr)
            ess.append(model.vector_ess)
        return numpy.mean(ess), numpy.std(loglr)

    def test_ess_is_reported(self):
        numpy.random.seed(0)
        model = self.build(1000)
        model.update(**self.point)
        model.loglr  # the marginalization, which fills in the ess
        self.assertTrue(numpy.isfinite(model.vector_ess))
        self.assertGreater(model.vector_ess, 0)
        self.assertLessEqual(model.vector_ess, 1000 + 1)

    def test_ess_rises_with_points(self):
        e256, _ = self.ess_and_spread(256)
        e4096, _ = self.ess_and_spread(4096)
        self.assertGreater(e4096, e256)

    def test_one_over_root_ess_predicts_the_spread(self):
        """The point of the number: it forecasts the error.

        One over the root of the effective sample size should track the
        spread of the log likelihood across independent draws, to within
        a factor. It is a proxy, not an identity: the exact constant
        depends on the integrand, and both quantities are themselves
        estimated from a finite number of runs, so the constant is not
        pinned tightly and should not be asserted to be. What must hold is
        that the proxy stays the same order of magnitude as the spread it
        stands in for, at every count.
        """
        for npoint in (256, 1024, 4096):
            ess, spread = self.ess_and_spread(npoint)
            predicted = 1.0 / ess ** 0.5
            ratio = spread / predicted
            self.assertGreater(ratio, 0.1,
                               "1/sqrt(ess)=%.4f far above spread=%.4f at "
                               "%d points" % (predicted, spread, npoint))
            self.assertLess(ratio, 10.0,
                            "1/sqrt(ess)=%.4f far below spread=%.4f at "
                            "%d points" % (predicted, spread, npoint))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargESS))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
