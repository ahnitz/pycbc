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
"""Tests the warning when the distance interpolator range is exceeded.

The distance-marginalized likelihood is interpolated over a fixed range
of signal to noise ratio; outside it the value is dropped to zero, which
biases the result. That used to happen silently. This checks that the
first time it happens the run says so, once, and that staying inside the
range says nothing.
"""

import logging
import unittest

import numpy
from utils import simple_exit

from pycbc.inference.models.tools import setup_distance_marg_interpolant


class TestDistanceBoundsWarning(unittest.TestCase):

    def interpolant(self, snr_range=(5, 40)):
        dist_rescale = numpy.linspace(0.5, 2.0, 60)
        weights = numpy.ones(60) / 60
        return setup_distance_marg_interpolant(
            (dist_rescale, weights), snr_range=snr_range, density=(60, 60))

    def test_in_range_is_silent(self):
        interp = self.interpolant()
        with self.assertNoLogs(level='WARNING'):
            # a signal to noise ratio of about 20, well inside (5, 40)
            interp(200.0, 100.0)
            interp(numpy.array([150.0, 300.0]), numpy.array([80.0, 120.0]))

    def test_out_of_range_warns_once(self):
        interp = self.interpolant()
        with self.assertLogs(level='WARNING') as caught:
            interp(1e9, 1e12)
            interp(1e9, 1e12)
            interp(numpy.array([1e9]), numpy.array([1e12]))
        # exactly one warning despite three excursions
        self.assertEqual(len(caught.records), 1)
        self.assertIn('snr_range', caught.records[0].getMessage().lower()
                      .replace('signal to noise ratio', 'snr_range'))

    def test_out_of_range_still_returns_minus_inf(self):
        """The warning must not change the value: still -inf out of range."""
        interp = self.interpolant()
        logging.disable(logging.CRITICAL)
        try:
            self.assertEqual(interp(1e9, 1e12), -numpy.inf)
            v = interp(numpy.array([200.0, 1e9]), numpy.array([100.0, 1e12]))
        finally:
            logging.disable(logging.NOTSET)
        self.assertTrue(numpy.isfinite(v[0]))
        self.assertEqual(v[1], -numpy.inf)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestDistanceBoundsWarning))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
