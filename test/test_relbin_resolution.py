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
"""Tests that the relative binning resolution check tracks the actual error.

These build models directly with parameters chosen to stress the
approximation, rather than running a sampler, and on simulated noise
rather than downloaded data, so the whole file runs in a few seconds.
"""

import unittest
from itertools import pairwise

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.inference.models import relbin
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform


class TestRelbinResolution(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # a short simulated signal is enough, and avoids downloading data
        tc = 1187008882.42840
        flow, seglen, srate = 25., 16, 2048
        delta_f = 1. / seglen
        flen = int(srate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        cls.flow = {'H1': flow, 'L1': flow}
        cls.static = {'mass1': 1.3757, 'mass2': 1.3757, 'f_lower': flow,
                      'approximant': 'TaylorF2', 'polarization': 0.,
                      'ra': 3.44615914, 'dec': -0.40808407, 'tc': tc}
        cls.variable = ['distance', 'inclination']
        cls.prior = JointDistribution(cls.variable, SinAngle(inclination=None),
                                      Uniform(distance=(10, 100)))
        cls.q = {'distance': 42.0, 'inclination': 2.5}

        # a signal at the test point, injected loud enough that the
        # resolution matters: the check that refines only where the answer
        # matters can only be exercised where there is a signal to get
        # wrong. Pure noise would leave the likelihood flat and nothing to
        # resolve.
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=flow, delta_t=1. / srate,
            mass1=cls.static['mass1'], mass2=cls.static['mass2'],
            distance=cls.q['distance'], inclination=cls.q['inclination'],
            coa_phase=0.)
        hp.start_time += tc
        hc.start_time += tc
        cls.psds, cls.data = {}, {}
        seed = 7
        for ifo in ['H1', 'L1']:
            ts = noise_from_psd(int(seglen * srate), 1. / srate, psd,
                                seed=seed)
            ts._epoch = tc - seglen / 2
            seed += 11
            signal = Detector(ifo).project_wave(
                hp, hc, cls.static['ra'], cls.static['dec'],
                cls.static['polarization'])
            cls.data[ifo] = ts.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd

    def model(self, epsilon, static=None, fiducial=None):
        return models.Relative(
            list(self.variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=static or self.static, prior=self.prior,
            fiducial_params=fiducial or {'mass1': 1.3756}, epsilon=epsilon)

    def test_error_falls_with_resolution(self):
        """Adding bins must resolve the ratio better, and at second order.

        The bin width scales with epsilon and the error of a linear
        interpolation goes as the width squared, so halving epsilon should
        cut the reported error by roughly four.
        """
        values = [self.model(eps).check_bin_resolution()
                  for eps in (1.0, 0.5, 0.25)]
        self.assertTrue(values[0] > values[1] > values[2],
                        "error should fall as bins are added: %s" % values)
        for coarse, fine in pairwise(values):
            self.assertGreater(coarse / fine, 2.,
                               "expected close to second order: %s" % values)

    def test_error_tracks_the_likelihood_error(self):
        """The reported error must predict the error it stands in for.

        Both are measured against a much finer binning of the same model,
        so this needs no external reference.
        """
        exact = self.model(0.005)
        exact.update(**self.q)
        reference = exact.loglr

        seen = []
        for epsilon in (1.0, 0.25, 0.05):
            model = self.model(epsilon)
            model.update(**self.q)
            seen.append((model.check_bin_resolution(),
                         abs(model.loglr - reference)))

        # the diagnostic is only useful if the two fall together
        for (ecoarse, lcoarse), (efine, lfine) in pairwise(seen):
            self.assertGreater(ecoarse, efine, "diagnostic: %s" % seen)
            self.assertGreater(lcoarse, lfine, "likelihood: %s" % seen)

    def test_flags_a_poorly_placed_fiducial(self):
        """Bins that suit the fiducial need not suit the signal."""
        offset = dict(self.static)
        offset['mass1'], offset['mass2'] = 1.39, 1.36
        good = self.model(0.5).check_bin_resolution()
        bad = self.model(0.5, static=offset).check_bin_resolution()
        self.assertGreater(bad, good * 10)
        self.assertGreater(bad, 1e-3)

    def test_check_leaves_the_model_alone(self):
        """A diagnostic must not move the point the model is sitting at."""
        model = self.model(0.5)
        model.update(**self.q)
        before = model.loglr
        model.check_bin_resolution()
        self.assertEqual(model.current_params['distance'], self.q['distance'])
        self.assertEqual(before, model.loglr)

    def test_check_is_cheap(self):
        """It must cost a handful of waveforms, not a full analysis.

        The point of measuring the error at the bin midpoints rather than
        on a fine grid is that everything is generated at the resolution
        of the bins. Counting the generator calls and the samples asked of
        them is what says so; a check that quietly reached for the data's
        own frequencies would still return a sensible number.
        """
        model = self.model(0.5)
        asked = []
        original = relbin.get_fd_waveform_sequence

        def record(sample_points=None, **kwargs):
            asked.append(len(sample_points))
            return original(sample_points=sample_points, **kwargs)

        relbin.get_fd_waveform_sequence = record
        try:
            value = model.check_bin_resolution(ndraw=3)
        finally:
            relbin.get_fd_waveform_sequence = original

        self.assertTrue(numpy.isfinite(value))
        self.assertGreaterEqual(value, 0.)

        # one waveform at the edges per distinct edge set, from the
        # likelihood's own call, plus one at the midpoints per detector
        per_draw = len(model.edge_unique) + len(model.data)
        self.assertEqual(len(asked), 3 * per_draw)

        # and every one of them at the resolution of the bins, which is
        # far short of the data the analysis works at
        nbins = max(len(model.fedges[ifo]) for ifo in model.data)
        ndata = min(len(model.f[ifo]) for ifo in model.data)
        self.assertLessEqual(max(asked), nbins)
        self.assertLess(sum(asked), ndata)

    # 'auto' keeps the bins good enough while the sampler runs

    def wander(self, model, n=4000, spread=0.05):
        """Evaluate near the peak, the way a settled sampler would."""
        for i in range(n):
            model.update(distance=self.q['distance'] + spread * (i % 7 - 3),
                         inclination=self.q['inclination'])
            assert numpy.isfinite(model.loglr)
        return model

    def test_auto_leaves_a_good_setup_alone(self):
        """Nothing should be rebuilt when nothing is wrong."""
        model = self.wander(self.model('auto'))
        self.assertEqual(model.rebins, 0)
        self.assertEqual(model.epsilon, 0.5)

    # a fiducial off enough from the template that the coarse bins cost
    # real accuracy (interpolation error 0.3, worth 5 in the log
    # likelihood at this signal's strength), while the template still
    # matches the injected signal so the likelihood stays loud there.
    # Offsetting the template instead would mismatch the data and leave
    # nothing to resolve; a much larger offset would break the relative
    # binning outright and drive the likelihood negative.
    BAD_FIDUCIAL = {'mass1': 1.37}

    def test_auto_refines_when_the_bins_are_not_good_enough(self):
        """A fiducial away from the signal should buy bins, while running.

        The setup cannot know this: it depends on where the sampler goes.
        """
        fixed = self.model(0.5, fiducial=self.BAD_FIDUCIAL)
        fixed.update(**self.q)
        auto = self.wander(self.model('auto', fiducial=self.BAD_FIDUCIAL))
        auto.update(**self.q)

        self.assertGreater(auto.rebins, 0)
        self.assertLess(auto.epsilon, 0.5)

        fine = self.model(0.002, fiducial=self.BAD_FIDUCIAL)
        fine.update(**self.q)
        self.assertLess(abs(auto.loglr - fine.loglr),
                        abs(fixed.loglr - fine.loglr))

    def test_auto_does_not_rebin_away_from_the_peak(self):
        """Bins are only worth adding where the answer matters.

        A point the sampler is passing through on its way somewhere better
        does not need the bins to describe it, and paying for it would
        mean paying for the worst place the sampler ever wandered.
        """
        model = self.model('auto', fiducial=self.BAD_FIDUCIAL)
        model.update(**self.q)
        # as if a far better point had already been seen, so that
        # everything from here on is out in the tails
        model.best_loglr = model.loglr + 1000.

        self.wander(model)
        self.assertEqual(model.rebins, 0)
        self.assertEqual(model.epsilon, 0.5)

    def test_auto_settles_and_only_ever_refines(self):
        """It must stop refining rather than keep going for as long as it
        runs, and never go back and forth, which would never settle."""
        model = self.model('auto', fiducial=self.BAD_FIDUCIAL)
        seen = []
        for _ in range(12):
            self.wander(model, n=600)
            seen.append(model.epsilon)
        self.assertEqual(seen, sorted(seen, reverse=True), "%s" % seen)

        after = model.rebins
        self.wander(model)
        self.assertEqual(model.rebins, after,
                         "still refining after settling")

    def test_the_screen_weighs_the_error_by_the_signal(self):
        """A ratio error small on its own is not small at a loud signal.

        The screen weighs the interpolation error by the signal to noise
        ratio, because that is what turns it into an error in the log
        likelihood. A well-matched fiducial with coarse bins leaves a ratio
        error near 6e-4: worth less than the accuracy at this fixture's own
        strength, and worth more than it at a loud one, so the ratio alone
        cannot decide. Checked on the screen's own quantities, so it does
        not need an injection louder than the fixture can hold.
        """
        model = self.model(0.5, fiducial={'mass1': 1.3756})
        model.update(**self.q)
        model.get_waveforms(model.current_params)
        error = model.interpolation_error_from_reference()

        # small as a ratio error
        self.assertLess(error, 1e-3)
        # worth skipping at this signal's own strength
        here_snr = (2.0 * max(model.loglr, 0.0)) ** 0.5
        self.assertLess(here_snr * error, model.accuracy)
        # not worth skipping at a loud one: the same bins cost more than
        # the accuracy asked for
        self.assertGreater(100.0 * error, model.accuracy)

    def test_diagnostic_survives_vector_valued_parameters(self):
        """A marginalized or batched parameter is a vector, not a number.

        The model can be sitting at a whole set of points at once: a
        marginalized parameter is held as the vector of samples drawn for
        it, and some samplers hand a batch of points in together. The
        diagnostic generates a waveform of its own, and the generators take
        one point at a time, so it has to reduce them first.
        """
        model = self.model(0.5)
        model.update(**self.q)
        model.get_waveforms(model.current_params)
        expected = model.interpolation_error_from_reference()

        # what a marginalized run or a batched sampler leaves behind:
        # the same values, held as vectors rather than as numbers
        for name in ('distance', 'inclination', 'tc'):
            if name in model.current_params:
                model.current_params[name] = numpy.full(
                    64, model.current_params[name])

        value = model.interpolation_error_from_reference()
        self.assertTrue(numpy.isscalar(value) or numpy.ndim(value) == 0)
        self.assertGreaterEqual(value, 0.)
        # the parameters made vectors here either scale the ratio, which
        # cancels in a relative error, or never reach the generator
        self.assertAlmostEqual(value, expected, places=10)


    def test_a_model_without_a_prior_can_be_built(self):
        """The check runs at construction and must not need a prior."""
        model = models.Relative(
            list(self.variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=self.static, fiducial_params={'mass1': 1.3756},
            epsilon=0.5)
        model.update(**self.q)
        self.assertTrue(numpy.isfinite(model.loglr))
        self.assertEqual(model.check_bin_resolution(), 0.)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinResolution))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
