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
from validation import FLOW, TC, get_seed, make_data

from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models


class TestMargESS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed = get_seed(cls.__name__)
        cls.data, cls.psds, inj = make_data(cls.seed, ifos=['H1', 'L1'])
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

        The spread of the log likelihood across independent draws should
        be about one over the root of the effective sample size. Checked
        at a few counts; the constant of proportionality is order one and
        must not drift with the count, or the number would not be an
        effective sample size.
        """
        ratios = []
        for npoint in (256, 1024, 4096):
            ess, spread = self.ess_and_spread(npoint)
            predicted = 1.0 / ess ** 0.5
            ratios.append(spread / predicted)
        # the ratio must be order one and roughly constant across counts,
        # not growing or shrinking with them
        ratios = numpy.array(ratios)
        self.assertTrue((ratios > 0.2).all() and (ratios < 5.0).all(),
                        "spread / (1/sqrt(ess)) = %s" % ratios)
        self.assertLess(ratios.max() / ratios.min(), 4.0,
                        "the relation drifts with count: %s" % ratios)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargESS))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
