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
"""Every way of drawing the sky has to describe the same sky.

The pregenerated map and the analytic draw propose differently and carry
different weights, and the tilt proposes differently again. All three are
importance sampling the same posterior, so the sky they reconstruct --
proposals weighted by what they earn -- has to be the same one. This is
the test that says the proposal and its weight belong together; a draw
that placed samples well but weighted them wrongly would pass a test of
where the samples are and fail here.

Deliberately compares the *shape*, normalized. The three do not agree on
the total, and that difference is a known one in the map estimator's
normalization rather than anything about the draw -- see the branches
marg-sky-normalization and marg-sky-absolute. Asserting on the total here
would tie this test to that bug.
"""

import copy
import unittest

import numpy

from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import (JointDistribution, SinAngle, CosAngle,
                                 Uniform)
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42
# the signal has to fit inside the segment: from 45 Hz this
# waveform runs 64 s, and a 32 s segment would hold a third of
# it while the template was charged for all of it
FLOW, SEGLEN, SRATE = 45., 128, 2048
IFOS = ['H1', 'L1']
# V1 is built too, but only the inert-case test uses it: the
# comparison between draws is a two-detector statement
EXTRA = ['V1']
SKY = dict(ra=6.0062, dec=1.1195, polarization=0.3)
POINT = dict(distance=45., inclination=0.5)
VARIABLE = ['distance', 'inclination', 'tc', 'ra', 'dec']
SEEDS, VSAMPLES = 6, 6000


class TestSkyDrawsAgree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=45., inclination=0.5,
                                 coa_phase=1.1)
        hp.start_time += TC
        hc.start_time += TC
        psd = aLIGOZeroDetHighPower(int(SRATE * SEGLEN / 2) + 1,
                                    1. / SEGLEN, FLOW)
        cls.data, cls.psds = {}, {}
        for n, ifo in enumerate(IFOS + EXTRA):
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=11 + n)
            noise._epoch = TC - SEGLEN / 2
            sig = Detector(ifo).project_wave(hp, hc, SKY['ra'], SKY['dec'],
                                             SKY['polarization'])
            cls.data[ifo] = noise.add_into(sig).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.static = dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                          approximant='TaylorF2',
                          polarization=SKY['polarization'])
        cls.fiducial = dict(mass1=1.4, tc=TC, **SKY)
        cls.drawn = {name: cls.sky(name, **kw) for name, kw in (
            ('map', {}),
            ('analytic', dict(marginalize_sky_analytic=True)),
            ('tilted', dict(marginalize_sky_analytic=True,
                            marginalize_sky_amplitude=True)))}

    @classmethod
    def sky(cls, _name, **kwargs):
        """The sky one draw reconstructs: its points, weighted."""
        ra, sdec, wt = [], [], []
        for seed in range(SEEDS):
            numpy.random.seed(900 + seed)
            prior = JointDistribution(
                list(VARIABLE), SinAngle(inclination=None),
                Uniform(distance=(10, 400)), Uniform(tc=(TC - .1, TC + .1)),
                Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None))
            model = models.RelativeTimeDom(
                list(VARIABLE),
                {i: copy.deepcopy(cls.data[i]) for i in IFOS},
                low_frequency_cutoff={i: FLOW for i in IFOS},
                psds={i: cls.psds[i] for i in IFOS},
                static_params=cls.static, prior=prior,
                fiducial_params=cls.fiducial, epsilon=0.1,
                marginalize_phase=True,
                marginalize_vector_params='tc,ra,dec',
                marginalize_vector_samples=VSAMPLES, sample_rate=4096,
                **kwargs)
            model.update(**POINT)
            # what a sample contributes to the marginalization is its own
            # likelihood times the correction its draw carries, and the
            # reconstruction path is where the model hands back the first
            # of those per sample rather than already summed
            model.reconstruct_vector = True
            try:
                like = numpy.asarray(model.loglr)
            finally:
                model.reconstruct_vector = False
            vp = model.marginalize_vector_params
            total = like + numpy.asarray(model.marginalize_vector_weights)
            live = numpy.isfinite(total)
            ra.append(numpy.mod(numpy.asarray(vp['ra'])[live],
                                2 * numpy.pi))
            sdec.append(numpy.sin(numpy.asarray(vp['dec'])[live]))
            w = numpy.exp(total[live] - total[live].max())
            wt.append(w / w.sum())
        return (numpy.concatenate(ra), numpy.concatenate(sdec),
                numpy.concatenate(wt))

    def moments(self, name):
        ra, sdec, w = self.drawn[name]
        s = w.sum()
        mra, msd = (ra * w).sum() / s, (sdec * w).sum() / s
        vra = (w * (ra - mra) ** 2).sum() / s
        vsd = (w * (sdec - msd) ** 2).sum() / s
        return mra, msd, numpy.sqrt(vra), numpy.sqrt(vsd)

    def test_the_analytic_draw_finds_the_map_draw_s_sky(self):
        """The two propose completely differently. If the weights are
        right they still describe one sky."""
        a = self.moments('map')
        b = self.moments('analytic')
        for k, what in enumerate(('mean ra', 'mean sin(dec)',
                                  'width in ra', 'width in sin(dec)')):
            self.assertLess(
                abs(a[k] - b[k]), 0.25 * a[2 + (k % 2)],
                "%s: map gives %s, the analytic draw %s" % (what, a[k], b[k]))

    def test_tilting_does_not_move_the_sky(self):
        """The tilt changes where samples are proposed and takes the
        change back out of the weight, so the sky must not move."""
        a = self.moments('analytic')
        b = self.moments('tilted')
        for k, what in enumerate(('mean ra', 'mean sin(dec)',
                                  'width in ra', 'width in sin(dec)')):
            self.assertLess(
                abs(a[k] - b[k]), 0.25 * a[2 + (k % 2)],
                "%s: untilted %s, tilted %s" % (what, a[k], b[k]))

    def test_the_option_is_inert_where_there_is_no_azimuth(self):
        """One detector has no delay and three have no freedom left after
        them, so the tilt has nothing to weight in either. Turning it on
        must leave those draws exactly as they were -- not close, the
        same -- or it is reaching somewhere it does not belong.
        """
        for ifos in (['H1'], ['H1', 'L1', 'V1']):
            got = {}
            for tilt in (False, True):
                numpy.random.seed(900)
                prior = JointDistribution(
                    list(VARIABLE), SinAngle(inclination=None),
                    Uniform(distance=(10, 400)),
                    Uniform(tc=(TC - .1, TC + .1)),
                    Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None))
                data = {i: copy.deepcopy(self.data[i]) for i in ifos}
                model = models.RelativeTimeDom(
                    list(VARIABLE), data,
                    low_frequency_cutoff={i: FLOW for i in data},
                    psds={i: self.psds[i] for i in data},
                    static_params=self.static, prior=prior,
                    fiducial_params=self.fiducial, epsilon=0.1,
                    marginalize_phase=True,
                    marginalize_vector_params='tc,ra,dec',
                    marginalize_vector_samples=VSAMPLES, sample_rate=4096,
                    marginalize_sky_analytic=True,
                    marginalize_sky_amplitude=tilt)
                model.update(**POINT)
                got[tilt] = (model.loglr,
                             numpy.asarray(
                                 model.marginalize_vector_params['ra']).copy())
            self.assertEqual(got[False][0], got[True][0],
                             "%d detectors: the estimate moved" % len(ifos))
            self.assertTrue(
                numpy.array_equal(got[False][1], got[True][1]),
                "%d detectors: the sky drawn is not the same" % len(ifos))

    def test_the_draws_are_actually_different(self):
        """Without this the two tests above would pass on a draw that
        quietly fell back to the map."""
        same = 0
        for k in ('analytic', 'tilted'):
            if numpy.array_equal(self.drawn['map'][0][:50],
                                 self.drawn[k][0][:50]):
                same += 1
        self.assertEqual(same, 0,
                         "a draw produced the same points as the map")


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSkyDrawsAgree))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
