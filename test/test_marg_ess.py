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

It also pins down the arithmetic that produces the number, which is done
in linear space for speed rather than through logsumexp: the value has to
agree with the log-space route, has to survive log weights far outside the
range of exp, and has to give the obvious answers in the obvious cases.
"""

import copy
import unittest

import numpy
from scipy.special import logsumexp

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.inference.models.tools import marginalize_likelihood
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


def reported_ess(vloglr, logw):
    """The effective sample size marginalize_likelihood reports for a
    prescribed vector of log likelihood ratios.

    Feeding the wanted vector in as sh with hh zero makes the internal
    loglr exactly that vector, so the test drives the real code path
    instead of a copy of it.
    """
    vloglr = numpy.asarray(vloglr, dtype=float)
    _, ess = marginalize_likelihood(vloglr, numpy.zeros(len(vloglr)),
                                    logw=logw, return_ess=True)
    return ess


def logspace_ess(vloglr, logw):
    """The same quantity via logsumexp, kept as the reference."""
    lw = numpy.asarray(vloglr, dtype=float) + logw
    return float(numpy.exp(2.0 * logsumexp(lw) - logsumexp(2.0 * lw)))


class TestESSArithmetic(unittest.TestCase):
    """The effective sample size is now formed in linear space.

    The weights are exponentiated once, after their largest log is
    subtracted, and the ratio (sum w)^2 / sum w^2 is taken directly. That
    is much cheaper than two more logsumexp calls on an array already in
    hand, but it is not bit-identical to them: rounding differs at the
    last couple of digits. These tests fix how far it is allowed to drift
    and check the cases where a naive exponential would break.
    """

    def test_matches_the_log_space_form(self):
        """Agreement with logsumexp, over weights of every character.

        The offsets push the log weights far above and far below the range
        where exp is representable, which is the whole reason the maximum
        is subtracted first. The spread controls how peaked the weights
        are, from nearly uniform to dominated by a few points.
        """
        rng = numpy.random.RandomState(9)
        for n in (2, 17, 512, 4096):
            for spread in (1e-3, 1.0, 5.0, 40.0):
                for offset in (-900.0, 0.0, 700.0):
                    v = rng.normal(size=n) * spread + offset
                    logw = -numpy.log(n)
                    fast, ref = reported_ess(v, logw), logspace_ess(v, logw)
                    self.assertLess(
                        abs(fast / ref - 1.0), 1e-10,
                        "n=%d spread=%g offset=%g: %.17g vs %.17g"
                        % (n, spread, offset, fast, ref))

    def test_matches_with_uneven_weights(self):
        """The log weights need not be a single number per point.

        The marginalization passes a whole array of importance weights
        once the draws stop being uniform, so the two routes have to agree
        when both terms vary point to point.
        """
        rng = numpy.random.RandomState(11)
        n = 1000
        v = rng.normal(size=n) * 8.0 + 300.0
        logw = rng.normal(size=n) * 4.0
        logw -= logsumexp(logw)
        self.assertLess(abs(reported_ess(v, logw)
                            / logspace_ess(v, logw) - 1.0), 1e-10)

    def test_known_cases_come_out_right(self):
        """The two extremes and the ceiling between them."""
        # every point counting the same: the answer is the count itself
        for n in (1, 2, 37, 1024):
            ess = reported_ess(numpy.full(n, 12.5), -numpy.log(n))
            self.assertLess(abs(ess / n - 1.0), 1e-12)

        # one point carrying the integral: the answer rests on one sample
        v = numpy.full(64, -400.0)
        v[13] = 0.0
        self.assertLess(abs(reported_ess(v, -numpy.log(64)) - 1.0), 1e-12)

        # and between them the count is the ceiling, whatever the weights
        rng = numpy.random.RandomState(5)
        for n in (3, 64, 2048):
            for spread in (0.1, 10.0, 100.0):
                ess = reported_ess(rng.normal(size=n) * spread, -numpy.log(n))
                self.assertGreaterEqual(ess, 1.0 - 1e-12)
                self.assertLessEqual(ess, n + 1e-9)

    def test_weightless_points_are_not_points(self):
        """Excluded draws arrive as -inf and must contribute nothing."""
        rng = numpy.random.RandomState(21)
        v = rng.normal(size=40) * 3.0
        keep = numpy.zeros(40, dtype=bool)
        keep[:9] = True
        logw = numpy.where(keep, -numpy.log(9), -numpy.inf)
        self.assertLess(abs(reported_ess(v, logw)
                            - logspace_ess(v[keep], -numpy.log(9))), 1e-9)

        # with none of them finite there is nothing to count
        self.assertTrue(numpy.isnan(
            reported_ess(numpy.full(8, -numpy.inf), -numpy.log(8))))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargESS))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestESSArithmetic))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
