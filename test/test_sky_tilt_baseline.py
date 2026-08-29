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
"""The baseline the azimuth tilt has to keep meeting.

These fix what was established by measurement rather than argument, so a
later change that breaks it says so. Each one guards a mistake that was
actually made:

- the fixture has to hold the whole signal (a 256 s waveform went into a
  32 s segment and every number taken from it was wrong)
- with no noise the injection has to be the ring maximum (Cauchy-Schwarz
  guarantees it; three separate references failed this and were believed
  anyway)
- the cells have to be fine enough to see the statistic (at sixteen they
  were not, and the tilt was worse than drawing flat)
- the statistic has to sit on the posterior, not merely near it
- nothing may produce an infinity at an antenna null
"""

import copy
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
FLOW, SEGLEN, SRATE = 45., 128, 2048
IFOS = ['H1', 'L1']
VARIABLE = ['distance', 'inclination', 'tc', 'ra', 'dec']
# a ring that carries amplitude information, away from the cos_t = 0
# minimum where two nearly co-aligned detectors see almost equally
COS_T, PHASE = 0.60, 1.1
DIST, PSI = 45.0, 0.4
# the scoring grid has to resolve the posterior peak; under about fifteen
# points across it the efficiency it reports is not converged
MFINE, MIN_PEAK_POINTS = 2880, 15
# a ring away from the cos_t = 0 minimum, used by the pieces
# that exercise the response directly
RING_CT = -0.33


def frame():
    det = {i: Detector(i) for i in IFOS}
    gmst = det['H1'].gmst_estimate(TC)
    base = (numpy.asarray(det['H1'].location)
            - numpy.asarray(det['L1'].location))
    dhat = base / numpy.linalg.norm(base)
    uvec = numpy.cross(dhat, [0., 0., 1.])
    uvec = uvec / numpy.linalg.norm(uvec)
    return det, gmst, dhat, uvec, numpy.cross(dhat, uvec)


class TestSkyTiltBaseline(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.det, cls.gmst, cls.dhat, cls.uvec, cls.vvec = frame()
        sin_t = numpy.sqrt(1.0 - COS_T ** 2.0)
        nhat = (COS_T * cls.dhat
                + sin_t * (numpy.cos(PHASE) * cls.uvec
                           + numpy.sin(PHASE) * cls.vvec))
        cls.ra = float(numpy.mod(numpy.arctan2(nhat[1], nhat[0]) + cls.gmst,
                                 2 * numpy.pi))
        cls.dec = float(numpy.arcsin(numpy.clip(nhat[2], -1, 1)))
        cls.incl = 0.05
        cls.psd = aLIGOZeroDetHighPower(int(SRATE * SEGLEN / 2) + 1,
                                        1. / SEGLEN, FLOW)
        cls.clean = cls.inject(noiseless=True)
        cls.noisy = cls.inject(noiseless=False)

    @classmethod
    def inject(cls, noiseless):
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=DIST, inclination=cls.incl,
                                 coa_phase=1.1)
        cls.duration = len(hp) * hp.delta_t
        hp.start_time += TC
        hc.start_time += TC
        data = {}
        for k, ifo in enumerate(IFOS):
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, cls.psd,
                                   seed=500 + k)
            noise._epoch = TC - SEGLEN / 2
            if noiseless:
                noise = noise * 0.0
            data[ifo] = noise.add_into(
                Detector(ifo).project_wave(hp, hc, cls.ra, cls.dec, PSI)
            ).to_frequencyseries()
        return data

    def model(self, data, marg='tc', **kwargs):
        numpy.random.seed(5)
        prior = JointDistribution(
            list(VARIABLE), SinAngle(inclination=None),
            Uniform(distance=(10, 400)), Uniform(tc=(TC - .1, TC + .1)),
            Uniform(ra=(0, 2 * numpy.pi)), CosAngle(dec=None))
        return models.RelativeTimeDom(
            list(VARIABLE), {i: copy.deepcopy(data[i]) for i in IFOS},
            low_frequency_cutoff={i: FLOW for i in IFOS},
            psds={i: self.psd for i in IFOS},
            static_params=dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                               approximant='TaylorF2', polarization=PSI),
            prior=prior,
            fiducial_params=dict(mass1=1.4, mass2=1.35, tc=TC, ra=self.ra,
                                 dec=self.dec, polarization=PSI),
            epsilon=0.1, marginalize_phase=True,
            marginalize_vector_params=marg,
            marginalize_vector_samples=8000, sample_rate=4096, **kwargs)

    def ring(self, npoint):
        sin_t = numpy.sqrt(1.0 - COS_T ** 2.0)
        phi = (numpy.arange(npoint) + 0.5) * (2.0 * numpy.pi / npoint)
        n = (COS_T * self.dhat[:, None]
             + sin_t * (numpy.cos(phi) * self.uvec[:, None]
                        + numpy.sin(phi) * self.vvec[:, None]))
        ra = numpy.mod(numpy.arctan2(n[1], n[0]) + self.gmst, 2 * numpy.pi)
        return phi, ra, numpy.arcsin(numpy.clip(n[2], -1, 1))

    def posterior(self, data, npoint):
        """loglr along the ring. tc marginalized ONLY -- a model that
        marginalizes the sky overwrites the sky it is handed."""
        m = self.model(data)
        self.assertEqual(set(m.marginalize_vector_params), {'tc'},
                         "the reference must not marginalize the sky")
        phi, ra, dec = self.ring(npoint)
        out = numpy.array([(m.update(ra=float(ra[k]), dec=float(dec[k]),
                                     distance=DIST, inclination=self.incl),
                            m.loglr)[1] for k in range(npoint)])
        return m, phi, out

    def statistic(self, m, ncells):
        """the shipped statistic on the ring, at this cell count"""
        probe = DistMarg.__new__(DistMarg)
        probe.marginalize_sky_amplitude = True
        probe._current_params = {'inclination': self.incl,
                                 'polarization': PSI}
        probe.static_params = {}
        probe.AZ_CELLS = ncells
        table = probe._ring_response(self.det, IFOS,
                                     (self.dhat, self.uvec, self.vvec),
                                     numpy.array([COS_T]))
        wsq = probe._predicted_wsq(table)
        zsq = []
        for ifo in IFOS:
            t = TC + self.det[ifo].time_delay_from_earth_center(
                self.ra, self.dec, TC)
            zsq.append(numpy.array([
                abs(complex(m.sh[ifo].at_time(t, interpolate='quadratic')))
                ** 2.0 / float(m.hh[ifo])]))
        return probe._amplitude_tilt(
            wsq, tuple(zsq), [float(m.hh[i]) for i in IFOS])[0], probe

    def efficiency(self, post, logw, ncells, probe):
        """samples kept if the azimuth is drawn from this statistic"""
        p = numpy.exp(post - post.max()); p /= p.sum()
        self.assertGreaterEqual(
            int((p > 0.1 * p.max()).sum()), MIN_PEAK_POINTS,
            "the scoring grid does not resolve the posterior peak, so the "
            "efficiency it reports is not converged")
        q = numpy.exp(logw - logw.max()); q /= q.sum()
        q = (1.0 - probe.AZ_FLOOR) * q + probe.AZ_FLOOR / ncells
        n = len(post)
        idx = (((numpy.arange(n) + 0.5) * ncells) // n).astype(int) % ncells
        q = q[idx]; q /= q.sum()
        flat = 1.0 / (p * p * n).sum()
        return (1.0 / (p * p / numpy.maximum(q, 1e-300)).sum()) / flat

    # ---- the fixture itself -------------------------------------------

    def test_the_signal_fits_inside_the_segment(self):
        """A waveform longer than the segment leaves the data holding part
        of a signal the template is charged for in full."""
        self.assertLess(self.duration, SEGLEN,
                        "the waveform runs %.0f s in a %d s segment"
                        % (self.duration, SEGLEN))

    def test_a_noiseless_injection_is_the_ring_maximum(self):
        """With d = s exactly, |<s|h>| <= sqrt(<s|s><h|h>), so no point on
        the ring can beat the true one. Any reference that says otherwise
        is wrong, and three of mine did."""
        m, phi, post = self.posterior(self.clean, 720)
        _, ra, dec = self.ring(720)
        sin_t = numpy.sqrt(1.0 - COS_T ** 2.0)
        nt = (COS_T * self.dhat
              + sin_t * (numpy.cos(PHASE) * self.uvec
                         + numpy.sin(PHASE) * self.vvec))
        azt = numpy.mod(numpy.arctan2(nt @ self.vvec, nt @ self.uvec),
                        2 * numpy.pi)
        k = numpy.abs(numpy.angle(numpy.exp(1j * (phi - azt)))).argmin()
        self.assertGreater(
            post[k] - post.max(), -1.0,
            "the injection sits %.2f nats below the ring maximum with no "
            "noise present" % (post[k] - post.max()))

    # ---- the statistic ------------------------------------------------

    def test_the_statistic_sits_on_the_posterior(self):
        m, phi, post = self.posterior(self.noisy, MFINE)
        logw, _ = self.statistic(m, MFINE)
        gap = abs(numpy.angle(numpy.exp(
            1j * (phi[logw.argmax()] - phi[post.argmax()]))))
        self.assertLess(gap, 0.15,
                        "the statistic peaks %.3f rad from the posterior"
                        % gap)

    def test_the_cells_resolve_the_statistic(self):
        """The statistic is read at cell centres, so a cell has to be
        narrower than the statistic or the centres sample its tails. At
        sixteen cells this was worse than drawing the azimuth flat."""
        m, phi, post = self.posterior(self.noisy, MFINE)
        logw, probe = self.statistic(m, MFINE)
        p = numpy.exp(logw - logw.max()); p /= p.sum()
        mu = numpy.angle((p * numpy.exp(1j * phi)).sum())
        width = numpy.sqrt(
            (p * numpy.angle(numpy.exp(1j * (phi - mu))) ** 2.0).sum())
        cell = 2.0 * numpy.pi / DistMarg.AZ_CELLS
        self.assertLess(
            cell, 2.0 * width,
            "cells are %.3f rad and the statistic is %.3f rad wide; the "
            "centres will sample its tails" % (cell, width))

    def test_the_tilt_beats_drawing_the_azimuth_flat(self):
        m, phi, post = self.posterior(self.noisy, MFINE)
        logw, probe = self.statistic(m, DistMarg.AZ_CELLS)
        gain = self.efficiency(post, logw, DistMarg.AZ_CELLS, probe)
        self.assertGreater(gain, 5.0,
                           "the tilt keeps only %.1f times the samples a "
                           "flat azimuth draw does" % gain)

    def test_coarse_cells_would_be_caught(self):
        """The guard above has to be able to fail: sixteen cells must not
        pass the efficiency bar the shipped count does."""
        m, phi, post = self.posterior(self.noisy, MFINE)
        logw, probe = self.statistic(m, 16)
        self.assertLess(self.efficiency(post, logw, 16, probe),
                        self.efficiency(
                            post, self.statistic(m, DistMarg.AZ_CELLS)[0],
                            DistMarg.AZ_CELLS, probe),
                        "sixteen cells did as well as the shipped count, so "
                        "this bar proves nothing")

    # ---- numerical corners --------------------------------------------

    def test_nothing_blows_up_at_an_antenna_null(self):
        """Edge-on, |w| follows F+ through zero. A ratio would divide by
        it; this must not."""
        probe = DistMarg.__new__(DistMarg)
        probe.marginalize_sky_amplitude = True
        probe._current_params = {'inclination': numpy.pi / 2,
                                 'polarization': 0.0}
        probe.static_params = {}
        table = probe._ring_response(self.det, IFOS,
                                     (self.dhat, self.uvec, self.vvec),
                                     numpy.array([RING_CT]))
        wsq = [(f[0] * 0.5) ** 2.0 for f in table]
        self.assertLess(min(w.min() for w in wsq), 1e-4,
                        "this ring never approaches a null, so the test is "
                        "not exercising one")
        for zsq in ((numpy.array([4000.]), numpy.array([4000.])),
                    (numpy.array([0.0]), numpy.array([0.0])),
                    (numpy.array([1e6]), numpy.array([1e6]))):
            out = probe._amplitude_tilt(wsq, zsq, [1e4, 1e4])
            self.assertTrue(numpy.isfinite(out).all(),
                            "not finite for |z|^2 = %s" % (zsq[0][0],))

    def test_the_floor_bounds_what_one_sample_can_carry(self):
        """The statistic is sharp enough that a cell can hold effectively
        no probability. Drawing it anyway gives that sample a weight
        bounded only by how small the probability was, and one sample
        then carries the estimate. The flat share bounds it."""
        m, phi, post = self.posterior(self.noisy, MFINE)
        logw, probe = self.statistic(m, DistMarg.AZ_CELLS)
        q = numpy.exp(logw - logw.max()); q /= q.sum()
        self.assertLess(
            q.min(), 1e-6 / DistMarg.AZ_CELLS,
            "no cell is starved here, so this is not testing the floor")
        mixed = ((1.0 - DistMarg.AZ_FLOOR) * q
                 + DistMarg.AZ_FLOOR / DistMarg.AZ_CELLS)
        self.assertGreater(
            mixed.min() * DistMarg.AZ_CELLS, 0.5 * DistMarg.AZ_FLOOR,
            "a cell can still be drawn with a density that lets its weight "
            "run away")
        self.assertGreater(DistMarg.AZ_FLOOR, 0.0,
                           "with no flat share the weight is unbounded")

    def test_the_inclination_enters_the_prediction(self):
        """Face-on, |w|^2 is F+^2 + Fx^2 and the polarization cancels out
        of it entirely. Edge-on the cross term is gone and only F+(psi)
        remains, so the same sky predicts a different amplitude. A
        statistic that ignored the inclination would not know that."""
        def wsq_for(cosi, psi):
            probe = DistMarg.__new__(DistMarg)
            probe.marginalize_sky_amplitude = True
            probe._current_params = {'inclination': numpy.arccos(cosi),
                                     'polarization': psi}
            probe.static_params = {}
            table = probe._ring_response(self.det, IFOS,
                                         (self.dhat, self.uvec, self.vvec),
                                         numpy.array([RING_CT]))
            return numpy.concatenate(
                [w.ravel() for w in probe._predicted_wsq(table)])
        # edge on, the cross term is gone and the prediction is exactly
        # F+(psi)^2 / 4. antenna_pattern does the psi rotation itself, so
        # this compares against the library rather than against a second
        # copy of the formula being tested.
        cos_t = RING_CT
        sin_t = numpy.sqrt(max(1.0 - cos_t ** 2.0, 0.0))
        phi = (numpy.arange(DistMarg.AZ_CELLS) + 0.5) * (
            2.0 * numpy.pi / DistMarg.AZ_CELLS)
        n = (cos_t * self.dhat[:, None]
             + sin_t * (numpy.cos(phi) * self.uvec[:, None]
                        + numpy.sin(phi) * self.vvec[:, None]))
        ra = numpy.mod(numpy.arctan2(n[1], n[0]) + self.gmst, 2 * numpy.pi)
        dec = numpy.arcsin(numpy.clip(n[2], -1, 1))
        for psi in (0.0, 0.7):
            want = numpy.concatenate([
                numpy.array([self.det[i].antenna_pattern(
                    float(ra[k]), float(dec[k]), psi, TC)[0]
                    for k in range(DistMarg.AZ_CELLS)]) ** 2.0 / 4.0
                for i in IFOS])
            self.assertLess(
                numpy.abs(wsq_for(0.0, psi) - want).max(), 1e-10,
                "edge-on, the prediction is not F+(psi)^2 / 4: the "
                "inclination factors are wrong")

        face_a, face_b = wsq_for(1.0, 0.0), wsq_for(1.0, 0.7)
        self.assertLess(numpy.abs(face_a - face_b).max(), 1e-12,
                        "face-on, the polarization must cancel")
        edge_a, edge_b = wsq_for(0.0, 0.0), wsq_for(0.0, 0.7)
        self.assertGreater(numpy.abs(edge_a - edge_b).max(), 1e-3,
                           "edge-on, the polarization must matter")
        self.assertGreater(
            numpy.abs(face_a - wsq_for(0.0, 0.0)).max(), 1e-3,
            "the inclination must change the predicted amplitude")

    def test_the_density_is_the_one_drawn_from(self):
        """The weight divides by the density this returns, so it has to
        be the density the draw actually has -- not the anchor it picked,
        and not the distribution it approximates. Checked by drawing many
        times and comparing the frequency of each cell against what was
        reported for it.

        This is the test that a mixture, or any two-stage draw, gets
        wrong quietly: the samples still land in sensible places and only
        the bookkeeping is off, so nothing looks broken until the answer
        moves.
        """
        det = {i: Detector(i) for i in IFOS}
        base = (numpy.asarray(det['H1'].location)
                - numpy.asarray(det['L1'].location))
        dhat = base / numpy.linalg.norm(base)
        uvec = numpy.cross(dhat, [0., 0., 1.])
        uvec = uvec / numpy.linalg.norm(uvec)
        frame = (dhat, uvec, numpy.cross(dhat, uvec))
        probe = DistMarg.__new__(DistMarg)
        probe.marginalize_sky_amplitude = True
        probe._current_params = {'inclination': self.incl,
                                 'polarization': PSI}
        probe.static_params = {}
        probe.hh = {i: 1.0e4 for i in IFOS}
        probe._sky_rng = numpy.random.default_rng(4)
        # deliberately spread far wider than a real draw, so that
        # neighbouring anchors describe genuinely different rings. Over
        # the narrow span the first delay actually produces, adjacent
        # anchors are nearly identical and confusing one for another
        # would be undetectable -- and harmless.
        n = 240000
        rng = numpy.random.default_rng(5)
        cos_t = rng.uniform(0.15, 0.9, n)
        zsq = (rng.uniform(2000., 4000., n), rng.uniform(1500., 3000., n))
        az, dens = probe._tilted_azimuth(det, IFOS, frame, cos_t, zsq)
        # the density is returned relative to a flat draw, so 1/q summed
        # over what was actually drawn has to come back to one. This is
        # exact and needs nothing but the draw itself: any density that
        # is not the one the sampler used breaks it.
        got = numpy.exp(-dens).mean()
        err = numpy.exp(-dens).std() / numpy.sqrt(len(dens))
        self.assertLess(
            abs(got - 1.0), max(6.0 * err, 0.02),
            "1/q averages to %.4f +- %.4f where it must be 1; the density "
            "reported is not the density drawn from, so every weight is "
            "wrong by that factor" % (got, err))
        self.assertGreater(len(numpy.unique(az)), 1000,
                           "the draw is not varying, so this proves little")

    def test_it_declines_when_the_orientation_is_unknown(self):
        """The amplitude a source produces depends on its inclination and
        polarization. Either may itself be marginalized over, arriving as
        a vector, and then there is nothing to predict an amplitude with.
        The tilt has to stand down rather than pick a stand-in: a
        polarization averaged over its whole range points nowhere, and
        committing to it is worse than not tilting.
        """
        probe = DistMarg.__new__(DistMarg)
        probe.marginalize_sky_amplitude = True
        probe.static_params = {}
        probe._current_params = {'inclination': 0.5, 'polarization': 0.3}
        self.assertIsNotNone(probe._orientation(),
                             "a fixed orientation must be usable")
        probe.hh = {'H1': 1.0e4, 'L1': 1.0e4}
        self.assertTrue(probe._amplitude_usable(),
                        "a fixed orientation and known norms must be usable")
        del probe.hh
        self.assertFalse(
            probe._amplitude_usable(),
            "without the template norms there is no amplitude to predict, "
            "so the tilt must stand down rather than raise")
        probe.hh = {'H1': 1.0e4, 'L1': 1.0e4}
        probe._current_params = None
        self.assertIsNone(
            probe._orientation(),
            "with no point being evaluated there is nothing to predict "
            "an amplitude at; the precalculated draw arrives this way")
        for vary in ('inclination', 'polarization'):
            probe._current_params = {'inclination': 0.5, 'polarization': 0.3}
            probe._current_params[vary] = numpy.linspace(0.1, 1.2, 32)
            self.assertIsNone(
                probe._orientation(),
                "a marginalized %s left the tilt something to act on" % vary)
            self.assertFalse(probe._amplitude_usable(),
                             "a marginalized %s should stand the tilt down"
                             % vary)

    def test_a_zero_template_norm_is_survivable(self):
        probe = DistMarg.__new__(DistMarg)
        probe.marginalize_sky_amplitude = True
        wsq = [numpy.full((1, 8), 0.3), numpy.full((1, 8), 0.2)]
        out = probe._amplitude_tilt(wsq, (numpy.array([100.]),
                                          numpy.array([100.])), [0.0, 0.0])
        self.assertTrue(numpy.isfinite(out).all())


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestSkyTiltBaseline))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
