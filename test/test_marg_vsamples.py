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
"""Tests choosing the number of marginalization samples while running.

How many points the marginalization needs depends on where the sampler is,
so the number is not something setup can settle: it is followed instead,
from the effective sample size the marginalization reports, against the
scatter in nats the model was asked to keep to.

The number is allowed to go down as well as up, which is what makes these
tests necessary rather than nice to have: a controller that can move both
ways can chase its own noise for a whole run and cost more than it saves.
So as well as checking that it moves in the right direction and stops
inside its bounds, these check that it stops -- that the number converges
and then stays put over a long run of evaluations, and that where it
converges to does not depend much on where it started.

The signal is a light binary black hole in its own noise realisation, not
one the other marginalization tests use, so passing means the rule carried
over rather than that a threshold was fitted to a fixture.
"""

import copy
import time
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import (CosAngle, JointDistribution, SinAngle,
                                 Uniform)
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
HALFWIDTH = 0.05
# the resolution is fixed throughout: what is being tested is the other
# knob, and a rate that leaves the peak on a few samples is where the
# number of draws has the most to do
RATE = 8192
INJ = dict(mass1=2.2, mass2=1.6, distance=90., inclination=1.2,
           ra=2.5, dec=-0.9, polarization=1.4, coa_phase=2.1)


class TestMargVSamples(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant='IMRPhenomD', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: INJ[k] for k in ('mass1', 'mass2', 'distance',
                                   'inclination', 'coa_phase')})
        # the generator puts coalescence at t=0, so this puts it at tc
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        seed = 17
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                seed=seed)
            seed += 41
            ts._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(
                hp, hc, INJ['ra'], INJ['dec'], INJ['polarization'])
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

        cls.flow = {ifo: FLOW for ifo in cls.data}
        cls.static = dict(mass1=INJ['mass1'], mass2=INJ['mass2'],
                          f_lower=FLOW, approximant='IMRPhenomD',
                          ra=INJ['ra'], dec=INJ['dec'],
                          polarization=INJ['polarization'])
        cls.point = {'distance': INJ['distance'],
                     'inclination': INJ['inclination']}

    def model(self, vsamples='auto', accuracy=0.05):
        variable = ['distance', 'inclination', 'tc']
        prior = JointDistribution(
            list(variable), SinAngle(inclination=None),
            Uniform(distance=(10, 300)),
            Uniform(tc=(TC - HALFWIDTH, TC + HALFWIDTH)))
        return models.RelativeTimeDom(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=prior,
            fiducial_params={'mass1': INJ['mass1'], 'mass2': INJ['mass2'],
                             'tc': TC},
            epsilon=0.1, marginalize_vector_params='tc', sample_rate=RATE,
            marginalization_accuracy=accuracy,
            marginalize_vector_samples=vsamples)

    def sky_model(self, pool, vsamples='auto', accuracy=0.005):
        """A model that draws its points from a precalculated pool.

        The pool is what makes changing the number of points cheap -- it is
        the same drawing, subsetted -- and it is also the ceiling, since the
        subset is taken without replacement and the antenna factors stored
        alongside it belong to those points and no others.
        """
        variable = ['distance', 'inclination', 'tc', 'ra', 'dec']
        prior = JointDistribution(
            list(variable), SinAngle(inclination=None),
            Uniform(distance=(10, 300)), Uniform(ra=(0, 2 * numpy.pi)),
            CosAngle(dec=None),
            Uniform(tc=(TC - HALFWIDTH, TC + HALFWIDTH)))
        static = {k: v for k, v in self.static.items()
                  if k not in ('ra', 'dec')}
        return models.RelativeTimeDom(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=static, prior=prior,
            fiducial_params={'mass1': INJ['mass1'], 'mass2': INJ['mass2'],
                             'tc': TC, 'ra': INJ['ra'], 'dec': INJ['dec']},
            epsilon=0.1, marginalize_vector_params='tc, ra, dec',
            sample_rate=RATE, marginalization_accuracy=accuracy,
            marginalize_vector_samples=vsamples,
            marginalize_sky_initial_samples=2e4,
            precalculate_marginalization_points=pool)

    def drive(self, model, npoint, point=None, seed=13):
        """Evaluate the same point over and over, as a sampler sitting on
        the peak would, and record what the number of samples does."""
        state = numpy.random.get_state()
        try:
            numpy.random.seed(seed)
            trail = []
            for _ in range(npoint):
                model.update(**(point if point else self.point))
                model.loglr
                trail.append(model.vsamples)
        finally:
            numpy.random.set_state(state)
        return trail

    def changes(self, trail):
        """Where the number of samples moved, and which way."""
        return [(i, numpy.sign(trail[i] - trail[i - 1]))
                for i in range(1, len(trail)) if trail[i] != trail[i - 1]]

    def scatter(self, model, ndraw=64):
        """How far apart the same point lands on repeated evaluation."""
        state = numpy.random.get_state()
        try:
            numpy.random.seed(29)
            values = []
            for _ in range(ndraw):
                model.update(**self.point)
                values.append(float(model.loglr))
        finally:
            numpy.random.set_state(state)
        return numpy.std(values, ddof=1)

    def test_it_draws_more_when_the_measured_size_is_short(self):
        """The direction that matters.

        Too few points is noise in the likelihood that the sampler cannot
        tell from structure, so a run started well below what the accuracy
        needs has to climb, and to climb until the size it measures is the
        size that was asked for.
        """
        model = self.model()
        model.set_vsamples(32)
        self.drive(model, 300)
        self.assertGreater(model.vsamples, 32)
        self.assertGreater(model.vector_ess, model.wanted_ess())

    def test_it_draws_fewer_when_the_measured_size_is_past_the_margin(self):
        """The direction the resolution and the bins are not allowed.

        Points are only a subset of one drawing, so taking fewer in an easy
        region is a real saving and nothing has to be rebuilt to have it.
        What must not happen is going below what the accuracy asks for, so
        that is checked alongside.
        """
        model = self.model()
        settled = self.drive(model, 300)[-1]
        model.set_vsamples(settled * 20)
        trail = self.drive(model, 900)
        self.assertLess(trail[-1], settled * 20)
        self.assertGreater(model.vector_ess, model.wanted_ess())

    def test_it_settles_and_stays_settled(self):
        """The risk the design takes on by allowing decreases.

        A controller that can move both ways can sit on the boundary and
        cross it back and forth for a whole run. What stops that here is
        that the two decisions are separated by a margin, so this asserts
        what the margin is for: over a long run the number of samples stops
        changing, it gets there in a handful of steps, and it never turns
        around more than once on the way.
        """
        for start in ('auto', 32):
            model = self.model()
            if start != 'auto':
                model.set_vsamples(start)
            trail = self.drive(model, 1200)
            changes = self.changes(trail)
            turns = sum(1 for i in range(1, len(changes))
                        if changes[i][1] != changes[i - 1][1])
            self.assertLessEqual(len(changes), 8,
                                 "from %s: %s" % (start, changes))
            self.assertLessEqual(turns, 1, "from %s: %s" % (start, changes))
            self.assertEqual(trail[-500:], [trail[-1]] * 500,
                             "from %s: still moving at the end, %s"
                             % (start, changes))

    def test_where_it_settles_hardly_depends_on_where_it_started(self):
        """The margin costs something and it should be only that.

        Coming down stops as soon as the measured size is inside the margin
        rather than on the budget itself, so starting far above settles
        higher than starting below. That gap is the hysteresis, and it must
        be no wider than the margin that produced it -- otherwise the number
        of samples is a memory of the starting guess rather than a property
        of the signal. Coming up stops a tenth past the budget and coming
        down within the margin of four of it, so the width to expect is
        four over that tenth, a little under four.
        """
        below = self.model()
        below.set_vsamples(32)
        low = self.drive(below, 1200)[-1]

        above = self.model()
        above.set_vsamples(low * 20)
        high = self.drive(above, 1200)[-1]

        self.assertGreaterEqual(high, low * 0.9)
        self.assertLessEqual(high, low * 4.0, "%s against %s" % (high, low))

    def test_the_settled_number_delivers_the_accuracy(self):
        """Settling is not the point; settling in the right place is.

        The number of samples is followed from a promise about the scatter
        of the marginalized likelihood, so the scatter it ends up with has
        to be about the size promised. Checked loosely by half: the law it
        rests on is a fit good to some tens of a percent, and 64 draws pin
        a standard deviation to about a tenth.
        """
        accuracy = 0.05
        model = self.model(accuracy=accuracy)
        self.drive(model, 300)
        self.assertLess(self.scatter(model), 1.5 * accuracy,
                        "settled at %s samples" % model.vsamples)

    def test_a_tighter_accuracy_settles_at_more_samples(self):
        """The knob has to be a knob, and to point the right way."""
        settled = []
        for accuracy in (0.08, 0.04, 0.02):
            model = self.model(accuracy=accuracy)
            settled.append(self.drive(model, 400)[-1])
        self.assertTrue(settled[0] < settled[1] < settled[2], "%s" % settled)

    def test_it_stays_above_the_floor(self):
        """Below a few tens of samples the average being marginalized has a
        bias in its logarithm larger than any accuracy anyone asks for, so
        an accuracy loose enough to want fewer does not get fewer."""
        model = self.model(accuracy=5.0)
        trail = self.drive(model, 600)
        self.assertGreaterEqual(min(trail), 32)
        self.assertEqual(trail[-1], min(trail))

    def test_the_precalculated_pool_is_the_ceiling(self):
        """Growing the draw within the pool is nearly free; growing past it
        would mean building the pool again, and pairing the stored antenna
        factors up with points that were not the ones they belong to. An
        accuracy that asks for more than the pool holds therefore gets the
        pool, and is told so rather than being quietly under-served."""
        pool = 400
        model = self.sky_model(pool, accuracy=0.002)
        self.assertLessEqual(model.vsamples, pool)
        with self.assertLogs(level='WARNING') as logs:
            trail = self.drive(model, 300)
        self.assertLessEqual(max(trail), pool)
        self.assertTrue(any('no more samples' in m for m in logs.output),
                        "%s" % logs.output)

    def test_points_far_from_the_best_do_not_move_it(self):
        """What the marginalization needs in the tails is not what it needs
        at the peak, and a sampler spends much of its time passing through
        points whose answer does not matter. Their effective sample size is
        small for reasons that have nothing to do with the accuracy of the
        answer, so acting on them would run the number of samples up for
        nothing."""
        model = self.model(accuracy=0.02)
        self.drive(model, 300)
        settled = model.vsamples
        far = dict(self.point)
        far['distance'] = 295.
        self.drive(model, 900, point=far)
        self.assertLess(float(model.loglr), model.marg_best_loglr - 20.)
        self.assertEqual(model.vsamples, settled)

    def auto_rate_model(self, vsamples='auto', accuracy=0.01):
        """The same model with the rate left to the model as well."""
        variable = ['distance', 'inclination', 'tc']
        prior = JointDistribution(
            list(variable), SinAngle(inclination=None),
            Uniform(distance=(10, 300)),
            Uniform(tc=(TC - HALFWIDTH, TC + HALFWIDTH)))
        return models.RelativeTimeDom(
            list(variable), copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, prior=prior,
            fiducial_params={'mass1': INJ['mass1'], 'mass2': INJ['mass2'],
                             'tc': TC},
            epsilon=0.1, marginalize_vector_params='tc', sample_rate='auto',
            marginalization_accuracy=accuracy,
            marginalize_vector_samples=vsamples)

    def cost(self, model, ncall=40):
        """Seconds a likelihood evaluation takes."""
        numpy.random.seed(31)
        start = time.time()
        for _ in range(ncall):
            model.update(**self.point)
            model.loglr
        return (time.time() - start) / ncall

    def test_asking_for_both_does_not_buy_the_expensive_one(self):
        """The accuracy depends only on the product of resolution and
        draws, and the two are not the same price: the signal to noise
        series is rebuilt every call, so the rate is close to linear in
        cost, while quadrupling the draws costs tens of a percent. Asking
        for both must therefore not cost more than choosing a low rate by
        hand and asking only for the draws -- which is what happens if the
        rate is allowed to spend the whole budget before the draws are
        sized. A factor of two is allowed on a timing taken in one process
        on one signal; the point is that it is not the fourfold or more
        that leaning on the rate costs.
        """
        both = self.auto_rate_model()
        byhand = self.model(accuracy=0.01)
        self.drive(both, 200)
        self.drive(byhand, 200)
        together, apart = self.cost(both), self.cost(byhand)
        self.assertLess(together, 2.0 * apart,
                        "%s Hz and %s samples took %.2f ms against %.2f for "
                        "%s Hz and %s samples"
                        % (both.sample_rate, both.vsamples, together * 1e3,
                           apart * 1e3, byhand.sample_rate, byhand.vsamples))

    def test_asking_for_both_still_reaches_the_accuracy(self):
        """Spending less on the rate is only worth having if the accuracy
        still arrives, so the scatter is measured rather than predicted.

        This is also the configuration someone would actually write: the
        rate settled once at setup, the draws sized for it and then
        followed, and the two together delivering what was asked.
        """
        accuracy = 0.01
        model = self.auto_rate_model(accuracy=accuracy)
        trail = self.drive(model, 400)
        self.assertEqual(trail[-100:], [trail[-1]] * 100)
        self.assertLess(self.scatter(model), 1.5 * accuracy,
                        "%s Hz and %s samples"
                        % (model.sample_rate, model.vsamples))

    def test_an_explicit_count_leaves_the_rate_choice_alone(self):
        """A number of draws given explicitly is a number the model may not
        change, so the rate is the only thing left to buy the accuracy with
        and it has to buy all of it, exactly as before this option existed.
        Asking for the draws as well must then cost strictly less rate.
        """
        accuracy = 0.01
        explicit = self.auto_rate_model(vsamples=1000, accuracy=accuracy)
        # the rate alone, against the resolution the accuracy law asks of it
        self.assertGreaterEqual(
            explicit.resolved_samples(explicit.ref_snr),
            1.0 / (1000 * accuracy ** 2))
        chosen = self.auto_rate_model(accuracy=accuracy)
        self.assertLess(chosen.sample_rate, explicit.sample_rate)

    def test_the_rate_climbs_only_when_the_draws_run_out_of_room(self):
        """The rate is the lever of last resort.

        While there is room for the draws the resolution asked for is the
        least that can be trusted; when the budget needs more draws than
        there can be, the rate is the only thing left and it climbs. Which
        way round it went is visible in the draws: they are up against their
        ceiling in the second case and nowhere near it in the first.
        """
        from pycbc.inference.models.relbin import (MOST_VSAMPLES,
                                                   SANE_RESOLVED)
        loose = self.auto_rate_model(accuracy=0.01)
        self.assertLess(loose.resolved, 2.0 * SANE_RESOLVED)
        self.assertLess(loose.vsamples, MOST_VSAMPLES / 2)

        tight = self.auto_rate_model(accuracy=0.002)
        self.assertGreater(tight.resolved, SANE_RESOLVED)
        self.assertGreater(tight.sample_rate, loose.sample_rate)
        self.assertGreater(tight.vsamples, MOST_VSAMPLES / 2)

    def test_an_explicit_number_is_left_alone(self):
        """Asking for a number of samples must still get exactly that."""
        model = self.model(vsamples=777)
        self.assertFalse(model.adapt_vsamples)
        self.assertEqual(set(self.drive(model, 300)), {777})


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargVSamples))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
