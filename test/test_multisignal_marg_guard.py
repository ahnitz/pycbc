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
"""Tests that marginalization inside multi_signal is refused, not ignored.

The multi-signal likelihood adds the components' inner products and their
cross terms. A marginalized likelihood is an integral, not a sum of those
terms, so the two cannot simply be combined: the result is a wrong answer
rather than a slow one. The models knew this and said so in a log message,
but reported support anyway, which left the wrong answer reachable with no
sign that anything was off. It must be refused where it is asked for.
"""

import unittest

from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower

FLOW, SEGLEN, SRATE, TC = 25., 8, 2048, 1187008882.42840


class TestMultiSignalMargGuard(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # pure simulated noise is enough: the guard is about which
        # configurations are allowed, not about the value of the likelihood
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        cls.data, cls.psds = {}, {}
        seed = 13
        for ifo in ['H1', 'L1']:
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=seed)
            noise._epoch = TC - SEGLEN / 2
            seed += 41
            cls.data[ifo] = noise.to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}

    def model(self, **marg):
        """A single-template model with the marginalizations asked for."""
        return models.SingleTemplate(
            ['distance'], {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                               approximant='TaylorF2', tc=TC,
                               ra=1.7, dec=-0.4, polarization=0.3,
                               inclination=0.5),
            sample_rate=4096, **marg)

    def test_unmarginalized_model_is_supported(self):
        """The combination that does work must keep working."""
        model = self.model(marginalize_phase=False)
        self.assertIn(type(model), model.multi_signal_support)

    def test_marginalized_model_is_refused(self):
        """Asking for the broken combination must raise, not log and go on."""
        model = self.model(marginalize_phase=True)
        with self.assertRaises(ValueError) as caught:
            model.multi_signal_support
        self.assertIn('multi_signal', str(caught.exception))

    def test_the_refusal_survives_hasattr(self):
        """The refusal has to reach the caller.

        MultiSignalModel probes the property with hasattr, which swallows an
        AttributeError. A ValueError is not swallowed, so the refusal is
        seen rather than being read as 'this model has no support'.
        """
        model = self.model(marginalize_phase=True)
        with self.assertRaises(ValueError):
            hasattr(model, 'multi_signal_support')


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestMultiSignalMargGuard))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
