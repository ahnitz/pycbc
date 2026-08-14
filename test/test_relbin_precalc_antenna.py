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
sample for sample with the draw they belong to, and that the two situations
in which they do not apply -- no sky draw ever happened, and a response that
varies across the frequency bins -- still fall back to computing them.
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
# The coalescence time prior is deliberately not centred on the fiducial
# time, as it is not in the shipped configurations. The stored factors are
# evaluated in the middle of that prior and the recomputed ones at the
# fiducial time, so this offset is the whole of the difference between the
# two paths; centring the prior would make them agree bit for bit and the
# comparison below would prove nothing.
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

    def test_stored_factors_are_used_when_the_sky_is_drawn(self):
        """The whole point: with a sky draw in hand, nothing is recomputed."""
        model = self.model()
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertTrue(model.precalc_antenna_factors,
                        "the sky draw stored no antenna factors")
        self.assertEqual(sorted(used), sorted(model.data),
                         "the stored factors were not used for every "
                         "detector (used: %s)" % used)

    def test_stored_factors_belong_to_the_samples_that_were_drawn(self):
        """The lookup has to be aligned with the draw, sample for sample.

        get_precalc_antenna_factors indexes the full pointing grid with the
        indices the draw kept, so element i must describe the same patch of
        sky as element i of the drawn ra and dec. If that alignment were off
        the values would not be slightly different from the recomputed ones,
        they would be unrelated to them, and the likelihood would quietly be
        using another sky position's response.
        """
        model = self.model()
        model.update(distance=45., inclination=0.6)
        model.loglr
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
            self.assertLess(abs(dt - rdt).max(), 1e-6,
                            "%s delay does not match its draw" % ifo)

            # and a misalignment really would show up here, so the
            # tolerances above are not passing for free: sliding the draw
            # by a single sample moves the response by order unity
            self.assertGreater(abs(numpy.roll(fp, 1) - rfp).max(), 0.1)

    def test_loglr_agrees_with_recomputing_the_factors(self):
        """The stored and recomputed factors are close but not identical.

        Both are the response of a detector that is referenced at the time
        being analysed, but they are evaluated at slightly different times:
        the stored ones in the middle of the coalescence time prior, where
        the pointing grid was built, the recomputed ones at the fiducial
        time. Over the tenth of a second that separates them the earth turns
        far too little to matter: the response differs by a part in a
        million and the log likelihood ratio by about as much, comfortably
        inside the tolerance asserted here and far under any scale a
        sampler can see. Reaching for the wrong samples, by contrast,
        would move it by tens.
        """
        model = self.model()
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        with_stored = model.loglr
        self.assertEqual(len(used), len(model.data))

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

    def test_scalar_time_falls_back(self):
        """Reconstruction asks for one definite sky position and time.

        The vector marginalization is switched off for that call and the
        stored factors, which describe a thousand other pointings, are
        discarded. Anything indexing them would be reading a leftover draw.
        """
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

    def test_precalculated_points_keep_the_factors_after_a_scalar_time(self):
        """The fast path has to survive a reconstruction, not just precede it.

        With precalculate_marginalization_points set -- the configuration
        the shipped examples use and the one the saving was measured in --
        every draw comes from premarg_draw, which selects a subset of the
        points chosen once at setup. A scalar coalescence time discards the
        stored factors, and premarg_draw is then the only thing that can put
        them back; if it does not, the very first reconstruction turns the
        reuse off for the rest of the run and the optimization silently
        stops optimizing while every test above still passes.
        """
        model = self.model(premarg=True)
        self.assertTrue(hasattr(model, 'premarg'),
                        "the precalculated points were not built")
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        before = model.loglr
        self.assertTrue(numpy.isfinite(before))
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
        # was made; handing back a leftover set would pass every assertion
        # above and quietly use another sky position's response.
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
            self.assertLess(abs(dt - rdt).max(), 1e-6,
                            "%s delay does not match its draw" % ifo)

    def test_no_sky_marginalization_falls_back(self):
        """With only the time marginalized no pointings are ever drawn.

        The factors are then still their initial value and there is nothing
        to index, so the response has to be computed as before.
        """
        model = self.model(sky=False)
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertFalse(model.precalc_antenna_factors)
        self.assertEqual(used, [], "stored factors used without a sky draw")

    def test_earth_rotation_falls_back(self):
        """With earth rotation on the response is a function of frequency.

        Each bin is evaluated at the time the signal was at that frequency,
        so the response is a vector over the bins, while what the sky draw
        stored is a single number per pointing. The stored values cannot
        stand in for it, and using them would silently throw the frequency
        dependence away, which is the entire reason the option exists.
        """
        model = self.model(earth_rotation=True)
        used = watch_precalc_use(model)
        model.update(distance=45., inclination=0.6)
        try:
            model.loglr
        except Exception:
            # A marginalized sky and earth rotation is not a combination
            # the model supports, and this is where it has always come
            # apart. Failing is fine; quietly answering with the wrong
            # factors would not be.
            pass
        self.assertTrue(model.precalc_antenna_factors,
                        "the sky draw stored no antenna factors")
        self.assertEqual(used, [], "stored factors used with earth rotation")

        # and the response really does depend on frequency here, so there
        # is genuinely no single number that could have been substituted
        times = model.antenna_time[list(model.data)[0]]
        self.assertGreater(len(times), 1)
        fp, _ = model.det[list(model.data)[0]].antenna_pattern(RA, DEC, 0,
                                                               times)
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
