# Copyright (C) 2026 Alexander Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
"""Unit tests for the analytic delay-to-sky extrinsic draw.

These cover the pieces that can be checked against a closed form without
building a model or reading data, which is deliberately where the load is put:
a per-injection likelihood comparison is NOT a reliable regression test here,
because at some injections the marginalised loglr is not reproducible run to
run (see test_analytic_sky.py notes in the class docstrings).
"""
import unittest
import numpy
from scipy import stats

from pycbc.inference.models.tools import DistMarg
from utils import simple_exit


class TestAliasDraw(unittest.TestCase):
    """The alias table is the sampling primitive the candidate draw rests on.

    It must reproduce the law it is given exactly, and must never return a
    cell of zero weight -- underflowed SNR cells rely on that.
    """

    def setUp(self):
        self.rng = numpy.random.default_rng(7)
        self.ndraw = 400000

    def _check(self, weights, label):
        weights = numpy.asarray(weights, dtype=float)
        n = len(weights)
        prob = weights / weights.sum()
        idx = DistMarg._alias_draw(weights.copy(), self.ndraw)
        self.assertGreaterEqual(idx.min(), 0, label)
        self.assertLess(idx.max(), n, label)
        obs = numpy.bincount(idx, minlength=n).astype(float)
        # nothing may come back from a zero-probability cell
        self.assertEqual(obs[prob <= 0].sum(), 0,
                         "%s: drew a zero-weight cell" % label)
        exp = prob * self.ndraw
        keep = exp >= 5
        if keep.sum() > 1:
            chi2 = (((obs[keep] - exp[keep]) ** 2) / exp[keep]).sum()
            pval = stats.chi2.sf(chi2, int(keep.sum()) - 1)
            self.assertGreater(pval, 1e-4, "%s: chi2 p=%.3g" % (label, pval))

    def test_single_cell(self):
        self._check([3.7], "single cell")

    def test_uniform(self):
        # every table entry is exactly 1.0, which is the degenerate case for
        # the Vose build's two worklists
        self._check(numpy.ones(85), "uniform")

    def test_all_mass_on_one_cell(self):
        w = numpy.zeros(85)
        w[40] = 1.0
        self._check(w, "all mass on one cell")

    def test_mostly_underflowed(self):
        # the realistic shape: exp(loglr - max) underflows in most cells
        w = numpy.zeros(1065)
        w[500:640] = numpy.exp(self.rng.normal(0, 3, 140))
        self._check(w, "mostly underflowed")

    def test_extreme_dynamic_range(self):
        w = numpy.array([1e-300] * 84 + [1.0])
        self._check(w, "1e-300 dynamic range")

    def test_large(self):
        self._check(self.rng.random(16384) ** 4, "n=16384")


class TestMirrorRootSelection(unittest.TestCase):
    """Drawing the mirror root among those inside the tc prior must not change
    the expectation.

    The two roots of ``nhat = npar +- sgeo * e3`` are mirror images through the
    detector plane. Drawing one blind and rejecting it if its coalescence time
    falls outside the prior estimates ``(w+ v+ + w- v-)/2``. Drawing only among
    the valid roots, with the weight scaled by ``(v+ + v-)/2``, must estimate
    the same thing while discarding nothing.
    """

    def test_expectation_is_unchanged(self):
        rng = numpy.random.default_rng(3)
        n = 2000000
        okp = rng.random(n) < 0.55
        okd = rng.random(n) < 0.45
        lp = numpy.exp(rng.normal(0, 1.0, n))
        ld = numpy.exp(rng.normal(0, 1.0, n))
        w = numpy.exp(rng.normal(0, 0.5, n))
        target = 0.5 * w * (okp * lp + okd * ld)

        coin = rng.random(n) < 0.5
        old = numpy.where(numpy.where(coin, okp, okd),
                          w * numpy.where(coin, lp, ld), 0.0)

        nroot = okp.astype(numpy.int8) + okd
        up = numpy.where(okp & okd, coin, okp)
        new = numpy.where(nroot > 0,
                          w * (0.5 * numpy.maximum(nroot, 1))
                          * numpy.where(up, lp, ld), 0.0)

        for label, est in (("reject-after", old), ("draw-among-valid", new)):
            diff = est - target
            sem = diff.std(ddof=1) / numpy.sqrt(n)
            self.assertLess(abs(diff.mean() / sem), 5.0,
                            "%s biased: %+.3g +- %.3g" % (label, diff.mean(),
                                                          sem))
        # and the new scheme must never throw a sample away that has a root
        self.assertEqual((new[nroot > 0] == 0).sum(), 0)


class TestEpochOffsetPrecision(unittest.TestCase):
    """Delays must be formed as offsets from an epoch, not from absolute GPS.

    Each delay is a difference of two coalescence-time-like quantities. From
    absolute GPS those are order 1e9 s, which a double holds to about 2.4e-7 s
    -- a four-thousandth of a sample at 4096 Hz -- and that error lands on the
    cell boundaries the candidate weights are read from, so a drawn delay can
    be attributed to the wrong cell and weighted with its SNR.
    """

    def test_offset_form_is_far_more_accurate(self):
        from fractions import Fraction
        gps = 1000000000.235
        delta = 1.0 / 4096
        rng = numpy.random.default_rng(1)
        exact_delta = Fraction(1, 4096)
        abs_err = []
        off_err = []
        for _ in range(500):
            i0 = int(rng.integers(0, 87))
            cell = int(rng.integers(0, 87))
            dither = float(rng.uniform(-delta / 2, delta / 2))
            exact = Fraction(i0 - cell) * exact_delta + Fraction(dither)
            absolute = (gps + i0 * delta + dither) - (gps + cell * delta)
            offset = (i0 * delta + dither) - cell * delta
            abs_err.append(abs(Fraction(absolute) - exact))
            off_err.append(abs(Fraction(offset) - exact))
        worst_abs = float(max(abs_err))
        worst_off = float(max(off_err))
        # absolute GPS should be around 1e-4 of a cell, the offset form
        # negligible; require at least six orders of magnitude between them
        self.assertGreater(worst_abs / max(worst_off, 1e-30), 1e6,
                           "offset form is not better: %.3g vs %.3g"
                           % (worst_off, worst_abs))
        self.assertLess(worst_off / delta, 1e-10,
                        "offset form still loses a cell fraction: %.3g"
                        % (worst_off / delta))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAliasDraw))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMirrorRootSelection))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestEpochOffsetPrecision))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
