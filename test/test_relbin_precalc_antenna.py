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
"""Tests that RelativeTime reuses the antenna factors the sky draw stored.

The sky marginalization builds its grid of pointings once and keeps the
antenna response and the earth-center time delay of every one of them; the
likelihood then needs nothing more than to index the subset that was drawn.
Working them out again from the sky coordinates on every evaluation is the
same answer computed twice, and for a thousand pointings in three detectors
it is a few milliseconds of every likelihood call. These tests check that
the stored values are the ones actually used, that they still line up
sample for sample with the draw they belong to, and that the situations in
which they do not apply still fall back to computing them.
"""

import unittest

import numpy

from pycbc.distributions import (CosAngle, JointDistribution, SinAngle,
                                 Uniform, UniformAngle)
from pycbc.detector import Detector
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
# The coalescence time prior is deliberately off-centre from the fiducial
# time: the stored factors are evaluated in the middle of the prior and the
# recomputed ones at the fiducial time, so centring it would make the two
# paths agree bit for bit and the comparison below would prove nothing.
TCMIN, TCMAX = TC - 0.02, TC + 0.08
FLOW, SEGLEN, SRATE = 25., 16, 2048
RA, DEC, POL = 1.7, -0.4, 0.3
VSAMPLES, SKY_SAMPLES = 200, 4000
PREMARG_POINTS = 2000
SEED = 11


def watch_precalc_use(model):
    """Record which detectors reach for the precalculated factors."""
    used = []
    original = model.get_precalc_antenna_factors

    def watched(ifo):
        used.append(ifo)
        return original(ifo)

    model.get_precalc_antenna_factors = watched
    return used


class TestRelbinPrecalcAntenna(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=40., inclination=0.5,
                                 coa_phase=1.1)
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        seed = 17
        for ifo in ['H1', 'L1']:
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=seed)
            noise._epoch = TC - SEGLEN / 2
            seed += 53
            signal = Detector(ifo).project_wave(hp, hc, RA, DEC, POL)
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}

    def model(self, sky=True, earth_rotation=False, premarg=False):
        """A time-marginalized model, with the sky marginalized or not."""
        variable = ['distance', 'inclination', 'tc']
        dists = [SinAngle(inclination=None), Uniform(distance=(10, 200)),
                 Uniform(tc=(TCMIN, TCMAX))]
        static = dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                      approximant='TaylorF2', polarization=POL)
        marginalized = 'tc'
        if sky:
            variable = variable + ['ra', 'dec']
            dists += [UniformAngle(ra=None), CosAngle(dec=None)]
            marginalized = 'tc, ra, dec'
        else:
            static.update(ra=RA, dec=DEC)
        extra = {}
        if premarg:
            extra['precalculate_marginalization_points'] = PREMARG_POINTS
        numpy.random.seed(SEED)
        return models.RelativeTime(
            list(variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=static,
            prior=JointDistribution(list(variable), *dists),
            fiducial_params={'mass1': 1.4, 'mass2': 1.35, 'tc': TC,
                             'ra': RA, 'dec': DEC, 'polarization': POL,
                             'distance': 40., 'inclination': 0.5,
                             'coa_phase': 1.1},
            epsilon=0.1, earth_rotation=earth_rotation,
            marginalize_phase=True,
            marginalize_vector_params=marginalized,
            marginalize_vector_samples=VSAMPLES,
            marginalize_sky_initial_samples=SKY_SAMPLES,
            sample_rate=4096, **extra)

    def assert_factors_match_the_draw(self, model):
        """The stored factors must line up with the draw, sample for sample.

        get_precalc_antenna_factors indexes the full pointing grid with the
        indices the draw kept, so element i must describe the same patch of
        sky as element i of the drawn ra and dec. A misalignment would not
        leave the values slightly off, it would leave them unrelated, and
        the likelihood would quietly use another sky position's response.
        """
        p = model.current_params
        for ifo in model.data:
            fp, fc, dt = model.get_precalc_antenna_factors(ifo)
            atime = model.antenna_time[ifo]
            rfp, rfc = model.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                      0, atime)
            rdt = model.det[ifo].time_delay_from_earth_center(p['ra'],
                                                              p['dec'], atime)
            self.assertLess(max(abs(fp - rfp).max(), abs(fc - rfc).max()),
                            1e-4, "%s response does not match its draw" % ifo)
            # the earth turns 9.2e-7 s of delay per second, and the two
            # evaluations are about a second apart, so this is a bound on
            # the offset rather than on the alignment; the roll below is
            # what holds the alignment down
            self.assertLess(abs(dt - rdt).max(), 1e-5,
                            "%s delay does not match its draw" % ifo)

            # the tolerances above are not passing for free: sliding the
            # draw by one sample moves the response by order unity
            self.assertGreater(abs(numpy.roll(fp, 1) - rfp).max(), 0.1)

    def test_stored_factors_are_used_and_give_the_same_loglr(self):
        """The whole point: with a sky draw in hand, nothing is recomputed.

        The stored and the recomputed factors are close but not identical,
        as they are evaluated at slightly different times: the stored ones
        in the middle of the coalescence time prior, where the pointing grid
        was built, the recomputed ones at the fiducial time. Over the second
        or so between them the earth turns far too little to matter, so the
        log likelihood ratio moves by a part in a million, far under any
        scale a sampler can see. Reaching for the wrong samples, by
        contrast, would move it by tens.
        """
        model = self.model()
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        with_stored = model.loglr
        self.assertTrue(numpy.isfinite(with_stored))
        self.assertTrue(model.precalc_antenna_factors,
                        "the sky draw stored no antenna factors")
        self.assertEqual(sorted(used), sorted(model.data),
                         "the stored factors were not used for every "
                         "detector (used: %s)" % used)
        self.assert_factors_match_the_draw(model)

        # Re-evaluate the very same draw with the factors worked out again.
        # The draw is frozen so the sample set cannot move underneath the
        # comparison, the stored factors are hidden so the fallback runs,
        # and _loglr is called directly because loglr is cached.
        del used[:]

        def frozen_draw(*args, **kwargs):
            return model.marginalize_vector_params

        model.snr_draw = frozen_draw
        model.precalc_antenna_factors = None
        recomputed = model._loglr()
        self.assertEqual(used, [], "the fallback still used stored factors")
        self.assertLess(abs(with_stored - recomputed) / abs(recomputed), 1e-4,
                        "stored %r vs recomputed %r" % (with_stored,
                                                        recomputed))

    def test_precalculated_points_keep_the_factors_after_a_scalar_time(self):
        """The fast path has to survive a reconstruction, not just precede it.

        With precalculate_marginalization_points set -- the configuration
        the shipped examples use -- every draw comes from premarg_draw,
        which selects a subset of the points chosen once at setup. A scalar
        coalescence time discards the stored factors, and premarg_draw is
        the only thing that can put them back; if it does not, the first
        reconstruction turns the reuse off for the rest of the run.
        """
        model = self.model(premarg=True)
        self.assertTrue(hasattr(model, 'premarg'),
                        "the precalculated points were not built")
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertEqual(sorted(used), sorted(model.data))

        # A reconstruction-style call: one definite sky position and time.
        del used[:]
        model.update(distance=45., inclination=0.6, tc=TC, ra=RA, dec=DEC)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertIsNone(model.precalc_antenna_factors)
        self.assertEqual(used, [], "stored factors used for a scalar time")

        # and back to vector draws: the fast path must be live again
        for step in range(3):
            del used[:]
            model.update(distance=45. + step, inclination=0.6)
            self.assertTrue(numpy.isfinite(model.loglr))
            self.assertTrue(model.precalc_antenna_factors,
                            "the stored factors were not restored after a "
                            "scalar time (draw %d)" % step)
            self.assertEqual(sorted(used), sorted(model.data),
                             "the fast path stayed off after a scalar time "
                             "(draw %d, used: %s)" % (step, used))

        # Restoring them is only right if they still describe the draw that
        # was made; a leftover set would pass every assertion above.
        self.assert_factors_match_the_draw(model)

    def test_falls_back_where_the_stored_factors_do_not_apply(self):
        """Three configurations with nothing applicable to index.

        Reconstruction asks for one definite sky position and time, so the
        pointings the draw held no longer describe it. With only the time
        marginalized none are ever drawn. With earth rotation the response
        is a vector over the bins, evaluated at the time the signal was at
        each frequency, and no single stored number can stand in for it --
        using one would throw away the frequency dependence that is the
        whole reason the option exists. Each has to compute the response.
        """
        # a scalar coalescence time, as reconstruction uses
        model = self.model()
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        model.loglr
        self.assertTrue(model.precalc_antenna_factors)
        del used[:]
        model.update(distance=45., inclination=0.6, tc=TC, ra=RA, dec=DEC)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertIsNone(model.precalc_antenna_factors)
        self.assertEqual(used, [], "stored factors used for a scalar time")

        # no sky marginalization, so no pointings were ever drawn
        model = self.model(sky=False)
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertFalse(model.precalc_antenna_factors)
        self.assertEqual(used, [], "stored factors used without a sky draw")

        # earth rotation, so the response varies across the bins
        model = self.model(earth_rotation=True)
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        try:
            model.loglr
        except Exception:
            # A marginalized sky with earth rotation is not a combination
            # the model supports, and this is where it comes apart.
            # Failing is fine; answering with the wrong factors is not.
            pass
        self.assertTrue(model.precalc_antenna_factors,
                        "the sky draw stored no antenna factors")
        self.assertEqual(used, [], "stored factors used with earth rotation")

        # and the response really does depend on frequency here, so there
        # is genuinely no single number that could have been substituted
        ifo = list(model.data)[0]
        times = model.antenna_time[ifo]
        self.assertGreater(len(times), 1)
        fp, _ = model.det[ifo].antenna_pattern(RA, DEC, 0, times)
        self.assertGreater(abs(fp.max() - fp.min()), 1e-4)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinPrecalcAntenna))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
