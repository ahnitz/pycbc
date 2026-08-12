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
drawn from. What it has to be depends on how wide the peak of that series
is, which depends on the signal, so no fixed value suits every analysis.
These check that the peak width is measured and the rate follows from it,
against direct integration of the unmarginalized model.
"""

import copy
import unittest

from test_marg_normalization import TC, TestMargNormalization
from utils import simple_exit

from pycbc.inference import models

TARGET = 2.5


class TestMargResolution(TestMargNormalization):
    """Reuses the injection and the direct integration from its base."""

    def relative_time(self, sample_rate, halfwidth=0.004):
        variable, prior = self.prior(halfwidth)
        return models.RelativeTime(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=prior,
            fiducial_params={'mass1': self.static['mass1'], 'tc': TC},
            epsilon=0.1, marginalize_vector_params='tc',
            sample_rate=sample_rate, marginalize_vector_samples=2000)

    def test_auto_resolves_the_peak(self):
        model = self.relative_time('auto')
        self.assertGreaterEqual(
            model.resolved_samples(model.ref_snr), TARGET,
            "chose %s Hz" % model.sample_rate)

    def test_auto_does_not_pay_for_more_than_it_needs(self):
        """Half the chosen rate must not already have been enough.

        The series is rebuilt at every likelihood evaluation, so a rate
        higher than the peak calls for is a running cost.
        """
        model = self.relative_time('auto')
        coarser = self.relative_time(model.sample_rate / 2)
        self.assertLess(coarser.resolved_samples(coarser.ref_snr), TARGET,
                        "%s Hz would have done" % coarser.sample_rate)

    def test_auto_matches_direct_integration(self):
        mean, error = self.marginalized(0.004, sample_rate='auto')
        self.assertAlmostEqual(mean, self.integrated(0.004),
                               delta=max(0.02, 4 * error))

    def test_an_unresolved_peak_is_what_goes_wrong(self):
        """The measurement has to be measuring the thing that matters.

        A rate that leaves the peak on about one sample is unreliable: the
        answer depends on where the peak falls between samples. This checks
        that such a rate exists below the chosen one and really is worse,
        so the criterion is not just a number that happens to pass.
        """
        chosen = self.relative_time('auto').sample_rate
        reference = self.integrated(0.004)

        errors = []
        rate = chosen / 2
        while rate >= 2048:
            model = self.relative_time(rate)
            if model.resolved_samples(model.ref_snr) < 1.5:
                mean, _ = self.marginalized(0.004, sample_rate=rate)
                errors.append(abs(mean - reference))
            rate /= 2

        self.assertTrue(errors, "no unresolved rate to compare against")
        self.assertGreater(max(errors), 0.05,
                           "an unresolved peak should hurt: %s" % errors)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargResolution))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
