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
"""Tests locking the time marginalization around the reference peak.

The region is chosen once, from a reference; the peak it was chosen for
moves as the parameters vary. These cover asking for a region of a given
width, and noticing when the peak has moved to the edge of it.
"""

import copy
import unittest

import numpy

from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.inference.models.tools import peak_offset
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42
FLOW, SEGLEN, SRATE = 25., 32, 2048
SAMPLE_RATE = 4096
VARIABLE = ['distance', 'inclination', 'tc']


class TestPeakLock(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=40., inclination=0.5, coa_phase=1.1)
        hp.start_time += TC
        hc.start_time += TC
        psd = aLIGOZeroDetHighPower(int(SRATE * SEGLEN / 2) + 1,
                                    1. / SEGLEN, FLOW)
        cls.data, cls.psds = {}, {}
        seed = 7
        for ifo in ['H1', 'L1']:
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=seed)
            noise._epoch = TC - SEGLEN / 2
            seed += 11
            signal = Detector(ifo).project_wave(hp, hc, 1.7, -0.4, 0.3)
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.static = dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                          approximant='TaylorF2', ra=1.7, dec=-0.4,
                          polarization=0.3)

    def model(self, **kwargs):
        numpy.random.seed(5)
        prior = JointDistribution(
            list(VARIABLE), SinAngle(inclination=None),
            Uniform(distance=(10, 200)), Uniform(tc=(TC - 0.1, TC + 0.1)))
        return models.SingleTemplate(
            list(VARIABLE), copy.deepcopy(self.data),
            low_frequency_cutoff={ifo: FLOW for ifo in self.data},
            psds=self.psds, static_params=self.static, prior=prior,
            marginalize_phase=True, marginalize_vector_params='tc',
            sample_rate=SAMPLE_RATE, marginalize_vector_samples=1000,
            peak_lock_snr=4.0, **kwargs)

    def test_the_region_is_the_width_asked_for(self):
        """peak_lock_time sets the width of the region directly."""
        for half in (0.002, 0.01, 0.05):
            model = self.model(peak_lock_time=half)
            for ifo in model.data:
                width = model.tend[ifo] - model.tstart[ifo]
                self.assertAlmostEqual(width, 2 * half,
                                       delta=2. / SAMPLE_RATE)

    def test_it_is_narrower_than_the_prior(self):
        """Locking has to restrict something, or it is doing nothing."""
        locked = self.model(peak_lock_time=0.002)
        for ifo in locked.data:
            self.assertLess(locked.tend[ifo] - locked.tstart[ifo], 0.2)

    def test_the_peak_sits_in_the_middle(self):
        """Locked on the reference peak, the weight is centred."""
        model = self.model(peak_lock_time=0.01)
        with self.assertNoLogs(level='WARNING'):
            model.update(distance=40., inclination=0.5)
            model.loglr

    def test_a_peak_at_the_edge_is_reported(self):
        """Move the region off the peak and it has to say so.

        This is what happens in a real run when the reference the region
        was locked around stops describing the current parameters.
        """
        model = self.model(peak_lock_time=0.01)
        width = {ifo: model.tend[ifo] - model.tstart[ifo]
                 for ifo in model.data}
        for ifo in model.data:
            model.tstart[ifo] += 0.45 * width[ifo]
            model.tend[ifo] = model.tstart[ifo] + width[ifo]

        with self.assertLogs(level='WARNING') as caught:
            model.update(distance=40., inclination=0.5)
            model.loglr
        messages = [r.getMessage() for r in caught.records]
        self.assertTrue(any('locked region' in m for m in messages),
                        "no warning for a peak at the edge: %s" % messages)

    def test_the_measure_tracks_where_the_weight_is(self):
        """Zero in the middle, one against an edge."""
        n = 101
        for position, expected in ((50, 0.0), (75, 0.5), (100, 1.0)):
            logweight = numpy.full(n, -numpy.inf)
            logweight[position] = 0.0
            self.assertAlmostEqual(peak_offset(logweight), expected, places=6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPeakLock))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
