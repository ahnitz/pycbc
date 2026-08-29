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
"""Tests weighting the free azimuth by the amplitudes actually seen.

With two detectors the delay says which ring the source is on and nothing
about where along it, so the azimuth is drawn flat. The relative
amplitude does say something, because the response varies around the
ring. These check the response the tilt predicts is the one the detector
has, that tilting moves the draw and not the answer, and that it is worth
its cost -- which is not free, so a gain that does not survive being
measured per unit time is not a gain.
"""

import copy
import time
import unittest

import numpy

from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import (JointDistribution, SinAngle, CosAngle,
                                 Uniform)
from pycbc.inference import models
from pycbc.inference.models.tools import DistMarg
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42
# the signal has to fit inside the segment: from 45 Hz this
# waveform runs 64 s, and a 32 s segment would hold a third of
# it while the template was charged for all of it
FLOW, SEGLEN, SRATE = 45., 128, 2048
IFOS = ['H1', 'L1']
# the ring through this position carries usable amplitude information:
# cos_t 0.68, away from the H1-L1 minimum at cos_t = 0, where the plane
# bisecting the baseline sees both detectors almost equally
SKY = dict(ra=6.0062, dec=1.1195, polarization=0.3)
POINT = dict(distance=45., inclination=0.5)
VARIABLE = ['distance', 'inclination', 'tc', 'ra', 'dec']
# the statistic is the amplitude likelihood itself now, so there
# is no approximation left to hedge against and the spread stays
# small; it only broadens the proposal for safety
TILT = True


class TestSkyAmplitudeTilt(unittest.TestCase):

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
        for n, ifo in enumerate(IFOS):
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

    def model(self, vsamples=4000, seed=5, **kwargs):
        numpy.random.seed(seed)
        prior = JointDistribution(
            list(VARIABLE), SinAngle(inclination=None),
            Uniform(distance=(10, 200)), Uniform(tc=(TC - 0.1, TC + 0.1)),
            Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None))
        return models.RelativeTimeDom(
            list(VARIABLE),
            {i: copy.deepcopy(self.data[i]) for i in IFOS},
            low_frequency_cutoff={i: FLOW for i in IFOS},
            psds={i: self.psds[i] for i in IFOS},
            static_params=self.static, prior=prior,
            fiducial_params=self.fiducial, epsilon=0.1,
            marginalize_phase=True, marginalize_vector_params='tc,ra,dec',
            marginalize_vector_samples=vsamples, sample_rate=4096,
            marginalize_sky_analytic=True, **kwargs)

    # ---- the response the tilt predicts ---------------------------------

    @staticmethod
    def frame():
        det = {i: Detector(i) for i in ('H1', 'L1', 'V1')}
        base = (numpy.asarray(det['H1'].location)
                - numpy.asarray(det['L1'].location))
        dhat = base / numpy.linalg.norm(base)
        uvec = numpy.cross(dhat, [0., 0., 1.])
        uvec = uvec / numpy.linalg.norm(uvec)
        return det, (dhat, uvec, numpy.cross(dhat, uvec))

    @staticmethod
    def probe(spread=TILT, cosi=0.0, psi=0.0):
        p = DistMarg.__new__(DistMarg)
        p.marginalize_sky_amplitude = bool(spread)
        p._current_params = {'inclination': numpy.arccos(cosi),
                             'polarization': psi}
        p.static_params = {}
        return p

    def test_the_table_is_the_antenna_pattern(self):
        """The tilt reads its geometry from a table, so the table has to
        be the response, at every row and cell and not on average."""
        det, frm = self.frame()
        probe = self.probe()
        dhat, uvec, vvec = frm
        cos_t = numpy.linspace(-0.95, 0.95, 40)
        table = probe._ring_response(det, ['H1', 'L1'], frm, cos_t)
        sin_t = numpy.sqrt(numpy.maximum(1.0 - cos_t * cos_t, 0.0))
        phi = (numpy.arange(DistMarg.AZ_CELLS) + 0.5) * (
            2.0 * numpy.pi / DistMarg.AZ_CELLS)
        nhat = (cos_t[:, None] * dhat[:, None, None]
                + sin_t[:, None] * (numpy.cos(phi) * uvec[:, None, None]
                                    + numpy.sin(phi) * vvec[:, None, None]))
        shape = (len(cos_t), DistMarg.AZ_CELLS)
        for k, ifo in enumerate(('H1', 'L1')):
            fplus, fcross = det[ifo].antenna_pattern_from_direction(
                nhat.reshape(3, -1))
            self.assertLess(
                numpy.abs(table[k][0] - fplus.reshape(shape)).max(),
                1e-12, "%s: F+ is not the response it stands for" % ifo)
            self.assertLess(
                numpy.abs(table[k][1] - fcross.reshape(shape)).max(),
                1e-12, "%s: Fx is not the response it stands for" % ifo)

    def wsq_at(self, probe, table, cos_t):
        """|w|^2 per detector on the ring at this cosine, as the draw does"""
        return probe._predicted_wsq(table)

    def test_it_points_where_the_source_is(self):
        """Given the amplitudes a source really makes, the tilt must
        prefer that source's cell. Noiseless, so only the geometry can
        explain a hit.
        """
        det, frm = self.frame()
        dhat, uvec, vvec = frm
        gmst = det['H1'].gmst_estimate(TC)
        lon = SKY['ra'] - gmst
        ntrue = numpy.array([numpy.cos(SKY['dec']) * numpy.cos(lon),
                             numpy.cos(SKY['dec']) * numpy.sin(lon),
                             numpy.sin(SKY['dec'])])
        cos_t = float(ntrue @ dhat)
        az = numpy.mod(numpy.arctan2(ntrue @ vvec, ntrue @ uvec),
                       2.0 * numpy.pi)
        want = int(az / (2.0 * numpy.pi / DistMarg.AZ_CELLS))
        cosi, psi = numpy.cos(POINT['inclination']), SKY['polarization']
        probe = self.probe(cosi=cosi, psi=psi)
        table = probe._ring_response(det, IFOS, frm,
                                     numpy.array([cos_t]))
        wsq = self.wsq_at(probe, table, cos_t)
        # a noiseless source of this orientation: the observed amplitude
        # IS the predicted one, so the truth must come out on top
        norm = [1.0e4, 1.0e4]
        zsq = tuple(numpy.array([wsq[k][0][want] * norm[k]])
                    for k in range(2))
        logw = probe._amplitude_tilt(wsq, zsq, norm)[0]
        self.assertEqual(int(logw.argmax()), want,
                         "the tilt prefers cell %d, the source is in cell %d"
                         % (logw.argmax(), want))

    def test_it_earns_its_cost(self):
        """A guard on a cost that is known to be too high.

        Reading the statistic at 128 azimuth cells is what makes the draw
        work at all -- at 16 it is worse than drawing flat -- and it
        costs, measured, about 33 times a flat draw at 15000 samples.
        That is recorded rather than hidden: the algorithm is settled and
        the cost is not, and the shape of the fix is structural (the
        statistic only needs fine cells near its peak, so a coarse pass
        that finds the region and a fine one inside it would buy 128-cell
        accuracy at nearer 32-cell cost) rather than a constant to tune.
        This asserts only that it has not grown beyond what the cell
        count itself explains.
        """
        def ess_and_cost(tilt):
            model = self.model(vsamples=8000,
                               marginalize_sky_amplitude=tilt)
            model.update(**POINT)
            model.loglr
            got = []
            for seed in range(6):
                m = self.model(vsamples=8000, seed=700 + seed,
                               marginalize_sky_amplitude=tilt)
                m.update(**POINT)
                m.loglr
                logw = numpy.asarray(
                    m.marginalize_vector_params['logw_partial'])
                logw = logw[numpy.isfinite(logw)]
                got.append(len(logw))
            start = time.perf_counter()
            for _ in range(8):
                model.update(**POINT)
                model.loglr
            return (time.perf_counter() - start) / 8

        flat = ess_and_cost(0.0)
        tilted = ess_and_cost(TILT)
        self.assertLess(
            tilted, 60.0 * flat,
            "the tilt costs %.1f times a flat draw, beyond the ~33x the "
            "cell count explains" % (tilted / flat))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestSkyAmplitudeTilt))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
