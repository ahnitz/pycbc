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
"""Tests that the SNR series is built only when the draw will use it.

Once the marginalization points are precalculated, the per-evaluation
draw takes them from that pool and never looks at the signal-to-noise
series, so building it every call is wasted. The inner products it is
built alongside feed the likelihood and must still be there. This checks
both: that the series is skipped when the pool exists and built when it
does not, and that the inner products are present either way.
"""

import copy
import unittest

import numpy

from pycbc.detector import Detector
from pycbc.distributions import (CosAngle, JointDistribution, SinAngle,
                                  Uniform, UniformAngle)
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

FLOW, SEGLEN, SRATE, TC = 25., 32, 2048, 1187008882.42840
INJ = dict(mass1=1.4, mass2=1.35, distance=60., inclination=0.5,
           ra=1.7, dec=-0.4, polarization=0.3, coa_phase=1.1)


class TestRelbinSNRSkip(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # a fixed binary neutron star signal in simulated noise; the test
        # is self-contained and does not depend on a shared fixture
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
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
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}

    def build(self):
        variable = ['distance', 'inclination', 'tc', 'polarization',
                    'ra', 'dec']
        dists = [SinAngle(inclination=None), Uniform(distance=(10, 200)),
                 Uniform(tc=(TC - 0.1, TC + 0.1)),
                 UniformAngle(polarization=None), UniformAngle(ra=None),
                 CosAngle(dec=None)]
        fid = {'mass1': INJ['mass1'], 'tc': TC, 'ra': INJ['ra'],
               'dec': INJ['dec'], 'polarization': INJ['polarization']}
        numpy.random.seed(11)
        return models.RelativeTimeDom(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                               f_lower=FLOW, approximant='TaylorF2'),
            prior=JointDistribution(list(variable), *dists),
            fiducial_params=fid, epsilon=0.1, marginalize_phase=True,
            marginalize_vector_params='tc,ra,dec,polarization',
            marginalize_vector_samples=1000, sample_rate=4096,
            marginalize_sky_initial_samples=3e4,
            precalculate_marginalization_points=20000)

    def test_series_built_only_without_the_precalc_pool(self):
        model = self.build()
        model.update(distance=INJ['distance'], inclination=INJ['inclination'])
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertTrue(hasattr(model, 'premarg'))

        wfs = model.get_waveforms(model.current_params)

        # with the pool the draw takes precalculated points, so the series
        # is not built; the inner products it is built with still are
        self.assertEqual(len(model.get_snr(wfs)), 0)
        self.assertEqual(len(model.sh), len(self.data))
        self.assertEqual(len(model.hh), len(self.data))

        # without a pool (setup, or drawing fresh each call) the series is
        # drawn from and must be built for every detector
        saved = model.premarg
        del model.premarg
        try:
            self.assertEqual(len(model.get_snr(wfs)), len(self.data))
            self.assertEqual(len(model.sh), len(self.data))
        finally:
            model.premarg = saved


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRelbinSNRSkip))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
