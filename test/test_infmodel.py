# Copyright (C) 2021 Alex Nitz
#
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are the unittests for pycbc.inference.models
"""

import unittest
import copy
from utils import simple_exit
import numpy
from scipy import special
from pycbc.catalog import Merger
from pycbc.psd import interpolate, inverse_spectrum_truncation, aLIGOZeroDetHighPower
from pycbc.noise import noise_from_psd
from pycbc.frame import read_frame
from pycbc.filter import highpass, resample_to_delta_t
from astropy.utils.data import download_file
from pycbc.inference import models
from pycbc.distributions import Uniform, JointDistribution, SinAngle, UniformAngle
from pycbc.waveform.waveform import FailedWaveformError

class TestModels(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ###### Get data for references analysis of 170817
        m = Merger("GW170817")
        ifos = ['H1', 'V1', 'L1']
        cls.psds = {}
        cls.data = {}

        for ifo in ifos:
            print("Processing {} data".format(ifo))

            # Download the gravitational wave data for GW170817
            # Not using dcc.ligo.org because it sometimes hangs in the test
            # suite. But original URL is commented below
            # url = "https://dcc.ligo.org/public/0146/P1700349/001/"
            url = "https://media.githubusercontent.com/media/gwastro/pycbc_data/master/"
            url += "{}-{}1_LOSC_CLN_4_V1-1187007040-2048.gwf"
            fname = download_file(url.format(ifo[0], ifo[0]), cache=True)
            ts = read_frame(fname, "{}:LOSC-STRAIN".format(ifo),
                            start_time=int(m.time - 260),
                            end_time=int(m.time + 40))
            ts = highpass(ts, 15.0)
            ts = resample_to_delta_t(ts, 1.0/2048)
            ts = ts.time_slice(m.time-112, m.time + 16)
            cls.data[ifo] = ts.to_frequencyseries()

            psd = interpolate(ts.psd(4), ts.delta_f)
            psd = inverse_spectrum_truncation(psd, int(4 * psd.sample_rate),
                                              trunc_method='hann',
                                              low_frequency_cutoff=20.0)
            cls.psds[ifo] = psd

        cls.static = {'mass1':1.3757,
                  'mass2':1.3757,
                  'f_lower':20.0,
                  'approximant':"TaylorF2",
                  'polarization':0,
                  'ra': 3.44615914,
                  'dec': -0.40808407,
                  'tc':  1187008882.42840,
                 }

        cls.variable = ['distance',
                    'inclination',
                      ]
        cls.flow = {'H1':25, 'L1':25, 'V1':25}
        inclination_prior = SinAngle(inclination=None)
        distance_prior = Uniform(distance=(10, 100))
        tc_prior = Uniform(tc=(m.time-0.1, m.time+0.1))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, inclination_prior,
                                       distance_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop('polarization')
        cls.variable2 = cls.variable + ['polarization']
        cls.prior2 = JointDistribution(cls.variable2, inclination_prior,
                                       distance_prior, pol)

        # set up for gated model tests
        cls.static3 = cls.static.copy()
        cls.static3['t_gate_start'] = cls.static3['tc']
        cls.static3['t_gate_end'] = cls.static3['tc'] + 0.05
        cls.variable3 = cls.variable
        cls.prior3 = JointDistribution(cls.variable3, inclination_prior,
                                       distance_prior)

        ###### Expected answers
        # Answer taken from marginalized gaussian model
        cls.q1 = {'distance':42.0, 'inclination':2.5}
        cls.a1 = 541.8235746138382

        # answer taken from brute marginize pol + phase
        cls.a2 = 542.581
        cls.pol_samples = 200

        # answer with gate applied, no normalization
        cls.a3 = -1246.1948739646468

    def test_base_phase_marg(self):
        model = models.MarginalizedPhaseGaussianNoise(
                                self.variable, copy.deepcopy(self.data),
                                low_frequency_cutoff=self.flow,
                                psds=self.psds,
                                static_params=self.static,
                                prior=self.prior,
                                )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=1e-3)

    def test_relative_phase_marg(self):
        model = models.Relative(self.variable, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds = self.psds,
                                 static_params = self.static,
                                 prior = self.prior,
                                 fiducial_params = {'mass1':1.3756},
                                 epsilon = .1,
                                )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.002)

    def test_single_phase_marg(self):
        model = models.SingleTemplate(
                        self.variable, copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static,
                        prior = self.prior,
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a1, model.loglr, delta=0.02)

    def test_single_pol_phase_marg(self):
        model = models.SingleTemplate(
                        self.variable2, copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior2,
                        marginalize_vector_samples = 1000,
                        marginalize_vector_params = 'polarization',
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a2, model.loglr, delta=0.04)

    def test_gated_gaussian_psd_opts(self):
        model_toeplitz = models.GatedGaussianNoise(
                                 self.variable3, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds=self.psds,
                                 static_params=self.static3,
                                 prior=self.prior3,)
        model_matmul = models.GatedGaussianNoise(
                                 self.variable3, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds=self.psds,
                                 static_params=self.static3,
                                 prior=self.prior3, paint_method='matmul',)
        model_toeplitz.update(**self.q1)
        model_matmul.update(**self.q1)
        # check likelihoods match calculated
        self.assertAlmostEqual(self.a3, model_toeplitz.loglr, delta=0.01)
        self.assertAlmostEqual(self.a3, model_matmul.loglr, delta=0.01)
        # check paint method is being set correctly
        self.assertEqual('toeplitz', model_toeplitz.paint_method)
        self.assertEqual('matmul', model_matmul.paint_method)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=400,
                        cores=1,
                        base_model='marginalized_polarization',
                        )
        model.update(**self.q1)
        self.assertAlmostEqual(self.a2, model.loglr, delta=0.002)


class TestWaveformErrors(unittest.TestCase):
    """Tests that models handle no waveform errors correctly."""

    @classmethod
    def setUpClass(cls):
        cls.psds = {}
        cls.data = {}
        # static params for the test
        tc = 1187008882.42840
        flow = 20
        cls.static = {
            'approximant':"IMRPhenomD",
            'mass1': 40.,
            'mass2': 40.,
            'polarization': 0,
            'ra': 3.44615914,
            'dec': -0.40808407,
            'tc': tc,
            'distance': 100.,
            'inclination': 2.5
            }
        cls.variable = ['spin1z', 'f_lower']
        ifos = ['H1', 'L1', 'V1']
        # generate the reference psd
        seglen = 8
        delta_f = 1./seglen
        sample_rate = 4096
        delta_t = 1./sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        # put non-zero values in the beginning and end of the psd
        # so the gating models will work
        psd[0:int(flow/delta_f+1)] = psd[int(flow/delta_f+1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {'H1': flow, 'L1': flow, 'V1': flow}
        # generate the noise
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static['tc'] - seglen/2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # setup priors
        spin_prior = Uniform(spin1z=(-1., 2.))
        flowbad = 4000.
        flower_prior = Uniform(f_lower=(flow, flowbad+100.))
        pol = UniformAngle(polarization=None)
        cls.prior = JointDistribution(cls.variable, spin_prior, flower_prior)

        # set up for marginalized polarization tests
        cls.static2 = cls.static.copy()
        cls.static2.pop('polarization')
        cls.variable2 = cls.variable + ['polarization']
        cls.prior2 = JointDistribution(cls.variable2, spin_prior, flower_prior,
                                       pol)
        # set up gated parameters
        staticgate = cls.static.copy()
        staticgate['t_gate_start'] = tc - 0.05
        staticgate['t_gate_end'] = tc
        cls.staticgate = staticgate
        # margpol
        staticgate2 = cls.static2.copy()
        staticgate2['t_gate_start'] = tc - 0.05
        staticgate2['t_gate_end'] = tc
        cls.staticgate2 = staticgate2
        # the parameters to test:
        # these parameters should pass
        cls.pass_params = {'spin1z': 0., 'f_lower': flow}
        # these parameters should trigger a NoWaveformError
        cls.nowf_params = {'spin1z': 0., 'f_lower': flowbad}
        # these parameters should cause a FailedWaveformError
        cls.fail_params = {'spin1z': 2., 'f_lower': flow}

    def _run_tests(self, model, check_pass=True, check_nowf=True,
                   check_failed=True, check_raises=True):
        # check that the model works
        if check_pass:
            model.update(**self.pass_params)
            self.assertTrue(numpy.isfinite(model.loglr))
        # check that a no waveform error is caught correctly
        if check_nowf:
            model.update(**self.nowf_params)
            self.assertEqual(model.loglr, -numpy.inf)
        # check that a failed waveform is caught correctly
        if check_failed:
            model.update(**self.fail_params)
            self.assertEqual(model.loglr, -numpy.inf)
        # check that an error is raised if ignore_failed_waveforms is False
        if check_raises:
            model.ignore_failed_waveforms = False
            model.update(**self.fail_params)
            with self.assertRaises((FailedWaveformError, RuntimeError)):
                model.loglr

    def test_base_phase_marg(self):
        model = models.MarginalizedPhaseGaussianNoise(
                                self.variable, copy.deepcopy(self.data),
                                low_frequency_cutoff=self.flow,
                                psds=self.psds,
                                static_params=self.static,
                                prior=self.prior,
                                ignore_failed_waveforms=True)
        self._run_tests(model)

    def test_relative_phase_marg(self):
        model = models.Relative(self.variable, copy.deepcopy(self.data),
                                 low_frequency_cutoff=self.flow,
                                 psds = self.psds,
                                 static_params = self.static,
                                 prior = self.prior,
                                 fiducial_params = {},
                                 #fiducial_params = {'mass1': 40.},
                                 epsilon = .1,
                                 ignore_failed_waveforms=True)
        # relative model doesn't respect flower, so no point in testing nowf
        self._run_tests(model, check_nowf=False)

    def test_brute_pol_phase_marg(self):
        # Uses the old polarization syntax untill we decide to remove it.
        # Untill then, this also tests that that interface stays working.
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=4,
                        cores=1,
                        base_model='marginalized_polarization',
                        ignore_failed_waveforms=True
                        )
        # we need to do the check raises test separately because the underlying
        # base model's ignore_failed_waveforms needs to be set
        self._run_tests(model, check_raises=False)
        model = models.BruteParallelGaussianMarginalize(
                        self.variable, data=copy.deepcopy(self.data),
                        low_frequency_cutoff=self.flow,
                        psds = self.psds,
                        static_params = self.static2,
                        prior = self.prior,
                        marginalize_phase=4,
                        cores=1,
                        base_model='marginalized_polarization',
                        ignore_failed_waveforms=False
                        )
        self._run_tests(model, check_pass=False, check_nowf=False,
                        check_failed=False, check_raises=True)

    def test_gated_gaussian_noise(self):
        model = models.GatedGaussianNoise(
            self.variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            prior=self.prior,
            ignore_failed_waveforms=True)
        self._run_tests(model)

    def test_gated_gaussian_margpol(self):
        model = models.GatedGaussianMargPol(
            self.variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate2,
            prior=self.prior,
            ignore_failed_waveforms=True)
        self._run_tests(model)


class TestMarginalizedPolModels(unittest.TestCase):
    """Tests that marginalized polarization models return the same
    loglikelihood as the unmarginalized model numerically integrated over the
    same set of polarization samples.
    """

    @classmethod
    def setUpClass(cls):
        # Use the same setup as TestWaveformErrors
        cls.psds = {}
        cls.data = {}
        tc = 1187008882.42840
        flow = 20
        cls.static = {
            'approximant': "IMRPhenomD",
            'mass1': 40.,
            'mass2': 40.,
            'ra': 3.44615914,
            'dec': -0.40808407,
            'tc': tc,
            'distance': 100.,
            'inclination': 2.5,
            'f_lower': flow
        }
        cls.marg_variable = ['spin1z']
        cls.orig_variable = cls.marg_variable + ['polarization']
        ifos = ['H1', 'L1', 'V1']
        seglen = 8
        delta_f = 1. / seglen
        sample_rate = 4096
        delta_t = 1. / sample_rate
        flen = int(sample_rate * seglen / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, delta_f, flow)
        psd[0:int(flow / delta_f + 1)] = psd[int(flow / delta_f + 1)]
        psd[-2:] = psd[-2]
        seed = 1000
        cls.flow = {'H1': flow, 'L1': flow, 'V1': flow}
        for ifo in ifos:
            tsamples = int(seglen * sample_rate)
            ts = noise_from_psd(tsamples, delta_t, psd, seed=seed)
            ts._epoch = cls.static['tc'] - seglen / 2
            seed += 1027
            cls.data[ifo] = ts.to_frequencyseries()
            cls.psds[ifo] = psd
        # Set up the static parameters for the gated models
        staticgate = cls.static.copy()
        staticgate['t_gate_start'] = tc - 0.5
        staticgate['t_gate_end'] = tc
        cls.staticgate = staticgate

    def _test_models(self, margpol_model, orig_model, polarization_samples):
        """Tests that manual integration of the loglikelihood over
        polarizations gives the same result as the marginalized model."""
        # Get the loglikelihood from the margpol model
        params = {'spin1z': 0.}
        margpol_model.update(**params)
        margpol_logl = margpol_model.loglikelihood
        # Numerically integrate the loglikelihood over polarizations
        logls = numpy.zeros(len(polarization_samples))
        norm = numpy.log(len(logls))
        for ii,pol in enumerate(polarization_samples):
            params['polarization'] = pol
            orig_model.update(**params)
            logls[ii] = orig_model.loglikelihood
        integrated_logl = special.logsumexp(logls) - norm
        # Assert that the loglikelihoods are close
        self.assertAlmostEqual(
            margpol_logl, integrated_logl, places=2,
            msg=f"Loglikelihood mismatch: MargPol={margpol_logl}, "
                f"Integrated={integrated_logl}"
        ) 

    def test_gaussian_models(self):
        """Tests that the GaussianNoise models are consistent."""
        margpol_model = models.MarginalizedPolarization(
            self.marg_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True
        )
        orig_model = models.GaussianNoise(
            self.orig_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.static,
            ignore_failed_waveforms=True
        )
        # now test
        polsamples = margpol_model.marginalize_vector_params['polarization']
        self._test_models(margpol_model, orig_model, polsamples)

    def test_gated_models(self):
        """Tests that the Gated models are consistent."""
        # Initialize the GatedGaussianMargPol model
        margpol_model = models.GatedGaussianMargPol(
            self.marg_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True
        )
        # Initialize the GatedGaussianNoise model
        orig_model = models.GatedGaussianNoise(
            self.orig_variable, data=copy.deepcopy(self.data),
            low_frequency_cutoff=self.flow,
            psds=self.psds,
            static_params=self.staticgate,
            ignore_failed_waveforms=True
        )
        # now test
        polsamples = margpol_model.pol
        self._test_models(margpol_model, orig_model, polsamples)

class TestDistanceDrawRefinement(unittest.TestCase):
    """A drawn distance should not carry the spacing of the grid it came
    from. The marginalization integrates the grid away; a draw does not."""

    DMIN, DMAX, S, H = 10., 1000., 300., 300.

    def grid(self, n):
        """The same three arrays setup_marginalization builds."""
        locs = numpy.linspace(self.DMIN, self.DMAX, n)
        w = locs ** 2                      # uniform in volume
        ref = 0.5 * (self.DMAX + self.DMIN)
        r = ref / locs
        return locs, numpy.log(w / w.sum()) + self.S * r - 0.5 * self.H * r ** 2

    def exact(self, m=200001):
        d = numpy.linspace(self.DMIN, self.DMAX, m)
        r = 0.5 * (self.DMAX + self.DMIN) / d
        lp = self.S * r - 0.5 * self.H * r ** 2 + 2 * numpy.log(d)
        p = numpy.exp(lp - lp.max())
        c = numpy.cumsum(p)
        return d, c / c[-1]

    def draws(self, n, ndraw=6000, locs=True, seed=3):
        from pycbc.inference.models.tools import draw_sample
        g, logp = self.grid(n)
        rng = numpy.random.default_rng(seed)
        if locs:
            return numpy.array([draw_sample(logp, rng=rng, locs=g)
                                for _ in range(ndraw)])
        return numpy.array([g[draw_sample(logp, rng=rng)]
                            for _ in range(ndraw)])

    def test_refined_draw_leaves_the_grid(self):
        """Without refinement every draw is one of 40 values; with it, not."""
        self.assertLessEqual(len(numpy.unique(self.draws(40, locs=False))), 40)
        self.assertGreater(len(numpy.unique(self.draws(40))), 2000)

    def test_refinement_beats_a_grid_four_times_finer(self):
        """80 refined points track the true distribution better than 320
        unrefined ones, so the spacing and not the count was the limit.
        Measured across five seeds; the gain runs out around eight times."""
        d, c = self.exact()

        def ks(x):
            x = numpy.sort(x)
            e = numpy.interp(x, d, c)
            n = len(x)
            return numpy.abs(e - numpy.arange(1, n + 1) / n).max()

        for seed in range(3):
            self.assertLess(ks(self.draws(80, seed=seed)),
                            ks(self.draws(320, locs=False, seed=seed)))

    def test_no_half_spacing_offset(self):
        """A point sits in the middle of the weight it carries, not at the
        end of it. Getting that wrong shifts every draw by half a spacing."""
        d, c = self.exact()
        truth = numpy.interp(0.5, c, d)
        n = 40
        spacing = (self.DMAX - self.DMIN) / (n - 1)
        err = numpy.median(self.draws(n, ndraw=20000)) - truth
        self.assertLess(abs(err), 0.25 * spacing)

    def test_phase_is_conditioned_on_the_drawn_distance(self):
        """The level below must see the distance that was drawn, not the
        grid point nearest it, or the two disagree about the amplitude."""
        from pycbc.inference.models.tools import DistMarg
        from pycbc.inference.models.base import ModelStats
        m = DistMarg.__new__(DistMarg)
        n = 40
        locs = numpy.linspace(self.DMIN, self.DMAX, n)
        w = locs ** 2
        m.dist_ref = 0.5 * (self.DMAX + self.DMIN)
        m.dist_locs = locs
        m.distance_marginalization = (m.dist_ref / locs, w / w.sum())
        m.reconstruct_inline = ['vector', 'distance', 'phase']
        m.marginalize_vector_params = {}
        m.marginalize_phase = True
        m.marginalize_rng = numpy.random.default_rng(4)
        m._current_stats = ModelStats()
        seen = {}

        def spy(sh, hh):
            seen['sh'] = sh
            return {'coa_phase': 0.0}, 0.0, 0
        m.draw_phase = spy

        sh_total, hh_total = 300.0 + 0.0j, 300.0
        m.draw_inline(sh_total, hh_total, None)
        drawn = m._current_stats.distance
        # the amplitude handed down is the one the drawn distance implies
        self.assertAlmostEqual(seen['sh'],
                               sh_total * (m.dist_ref / drawn), places=6)
        # and that distance is genuinely off-grid
        self.assertGreater(numpy.abs(locs - drawn).min(), 0.0)


    def test_index_draw_is_unchanged(self):
        """Asking for an index still returns one, for the callers that
        index the vector samples with it."""
        from pycbc.inference.models.tools import draw_sample
        _, logp = self.grid(40)
        rng = numpy.random.default_rng(1)
        i = draw_sample(logp, rng=rng)
        self.assertIsInstance(int(i), int)
        self.assertTrue(0 <= int(i) < 40)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestModels))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveformErrors))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMarginalizedPolModels))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDistanceDrawRefinement))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
