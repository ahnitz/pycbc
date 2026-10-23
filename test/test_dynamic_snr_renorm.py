import unittest

import numpy as np

import pycbc.filter
from pycbc.filter.matched_filter_ratio import (
    dynamic_snr_renormalize,
    get_dynamic_snr_renorm_factor,
)
from pycbc.types import Array, TimeSeries


class TestDynamicSNRRenormalize(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.dt = 1.0 / 2048.0
        self.sample_rate = 2048

    def test_export_availability(self):
        """Verify dynamic_snr_renormalize and get_dynamic_snr_renorm_factor are exported from pycbc.filter."""
        self.assertTrue(hasattr(pycbc.filter, "dynamic_snr_renormalize"))
        self.assertIs(pycbc.filter.dynamic_snr_renormalize, dynamic_snr_renormalize)
        self.assertTrue(hasattr(pycbc.filter, "get_dynamic_snr_renorm_factor"))
        self.assertIs(
            pycbc.filter.get_dynamic_snr_renorm_factor, get_dynamic_snr_renorm_factor
        )

    def test_empty_and_trivial(self):
        """Test empty arrays and trivial single-element inputs."""
        empty = np.array([], dtype=np.complex64)
        res = dynamic_snr_renormalize(empty, self.dt)
        self.assertEqual(len(res), 0)

        single = np.array([3.0 + 4.0j], dtype=np.complex64)
        res_single = dynamic_snr_renormalize(single, self.dt)
        self.assertEqual(len(res_single), 1)
        self.assertAlmostEqual(res_single[0], single[0])

    def test_invalid_dt(self):
        """Test non-positive dt returns series unaltered."""
        arr = np.array([1.0 + 1.0j], dtype=np.complex64)
        res = dynamic_snr_renormalize(arr, 0.0)
        np.testing.assert_array_equal(res, arr)

        res2 = dynamic_snr_renormalize(arr, -1.0)
        np.testing.assert_array_equal(res2, arr)

    def test_invalid_window_params(self):
        """Test that window_duration <= 2 * hollow_duration raises ValueError."""
        arr = np.array([1.0 + 1.0j], dtype=np.complex64)
        with self.assertRaises(ValueError):
            get_dynamic_snr_renorm_factor(
                arr, dt=self.dt, window_duration=1.0, hollow_duration=0.5
            )
        with self.assertRaises(ValueError):
            get_dynamic_snr_renorm_factor(
                arr, dt=self.dt, window_duration=0.8, hollow_duration=0.5
            )

    def test_types_and_dtypes(self):
        """Test that input type (numpy ndarray, PyCBC Array, TimeSeries) and dtypes are preserved."""
        N = 4096
        # complex64 numpy
        arr64 = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        out64 = dynamic_snr_renormalize(arr64, self.dt)
        self.assertIsInstance(out64, np.ndarray)
        self.assertEqual(out64.dtype, np.complex64)

        # complex128 numpy
        arr128 = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex128)
        out128 = dynamic_snr_renormalize(arr128, self.dt)
        self.assertIsInstance(out128, np.ndarray)
        self.assertEqual(out128.dtype, np.complex128)

        # PyCBC Array
        pycbc_arr = Array(arr64)
        out_pycbc_arr = dynamic_snr_renormalize(pycbc_arr, self.dt)
        self.assertIsInstance(out_pycbc_arr, Array)
        self.assertEqual(out_pycbc_arr.dtype, np.complex64)

        # PyCBC TimeSeries
        ts = TimeSeries(arr64, delta_t=self.dt)
        out_ts = dynamic_snr_renormalize(ts)
        self.assertIsInstance(out_ts, TimeSeries)
        self.assertEqual(out_ts.delta_t, self.dt)

    def test_stationary_gaussian_noise_unscaled(self):
        """Test that for unscaled complex SNR in stationary Gaussian noise, renorm factor is close to 1.0."""
        N = 16 * self.sample_rate
        # Real and imag each have variance 1.0 -> expected power 0.5 * (|real|^2 + |imag|^2) = 1.0
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        renorm = dynamic_snr_renormalize(
            noise, self.dt, window_duration=8.0, hollow_duration=0.5
        )

        factor = np.abs(renorm) / np.abs(noise)
        self.assertTrue(np.all(factor <= 1.0001))
        self.assertGreater(np.mean(factor), 0.97)

    def test_stationary_gaussian_noise_with_scale(self):
        """Test that an internally scaled reference buffer with scale=dt produces renorm factor close to 1.0."""
        N = 16 * self.sample_rate
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(
            np.complex64
        ) * self.dt
        renorm = dynamic_snr_renormalize(
            noise, self.dt, window_duration=8.0, hollow_duration=0.5, scale=self.dt
        )

        factor = np.abs(renorm) / np.abs(noise)
        self.assertTrue(np.all(factor <= 1.0001))
        self.assertGreater(np.mean(factor), 0.97)

    def test_signal_peak_preservation(self):
        """Test that a short high-SNR signal peak is preserved by the hollow exclusion."""
        N = 16 * self.sample_rate
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        center = N // 2
        # Inject high SNR signal of duration 0.2s (+/- 0.1s)
        half_width = int(0.1 / self.dt)
        noise[center - half_width : center + half_width] += 60.0 + 0j

        renorm = dynamic_snr_renormalize(
            noise, self.dt, window_duration=8.0, hollow_duration=0.5
        )
        factor_at_peak = (renorm[center] / noise[center]).real
        self.assertAlmostEqual(factor_at_peak, 1.0, places=2)

    def test_glitch_tail_downweighting(self):
        """Test that non-stationary glitch noise spanning multiple seconds is downweighted."""
        N = 16 * self.sample_rate
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        center = N // 2
        # Glitch from 1s to 3s after center (power 100x)
        g_start = center + int(1.0 / self.dt)
        g_end = center + int(3.0 / self.dt)
        noise[g_start:g_end] *= 10.0

        renorm = dynamic_snr_renormalize(
            noise, self.dt, window_duration=8.0, hollow_duration=0.5
        )
        factor_during_glitch = np.mean(
            np.abs(renorm[g_start:g_end]) / np.abs(noise[g_start:g_end])
        )
        self.assertLess(factor_during_glitch, 0.25)

    def test_variance_floor_no_amplification(self):
        """Test that quiet periods or zero arrays are never amplified (variance floor >= 1.0)."""
        zeros = np.zeros(2048, dtype=np.complex64)
        renorm = dynamic_snr_renormalize(zeros, self.dt, variance_floor=1.0)
        np.testing.assert_array_equal(renorm, zeros)

        # Very low amplitude noise
        quiet = (np.random.randn(2048) + 1j * np.random.randn(2048)).astype(
            np.complex64
        ) * 0.01
        renorm_quiet = dynamic_snr_renormalize(quiet, self.dt, variance_floor=1.0)
        factor = np.abs(renorm_quiet) / np.abs(quiet)
        self.assertTrue(np.allclose(factor, 1.0))

    def test_nan_inf_robustness(self):
        """Test that isolated NaN or Inf in the series does not poison the rest of the cumsum."""
        N = 8192
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        # Place isolated NaN and Inf
        noise[100] = np.nan
        noise[200] = np.inf

        rf = get_dynamic_snr_renorm_factor(noise, self.dt)
        # All outputs must be finite and in (0, 1]
        self.assertTrue(np.all(np.isfinite(rf)))
        self.assertTrue(np.all(rf > 0.0))
        self.assertTrue(np.all(rf <= 1.0))
        # Far from NaN/Inf points, mean factor remains close to 1.0
        self.assertGreater(np.mean(rf[1000:]), 0.95)

    def test_bbh_short_duration_preservation(self):
        """Verify that short-duration heavy BBH waveform SNR peaks are preserved (>= 99.5%)."""
        from pycbc.filter import matched_filter, sigmasq
        from pycbc.psd import aLIGOZeroDetHighPower
        from pycbc.waveform import get_td_waveform

        sample_rate = 2048
        dt = 1.0 / sample_rate
        N = 32 * sample_rate
        hp, _ = get_td_waveform(
            approximant="IMRPhenomD", mass1=30, mass2=30, f_lower=20.0, delta_t=dt
        )
        hp.resize(N)
        hp.start_time = 0.0

        psd = aLIGOZeroDetHighPower(N // 2 + 1, 1.0 / (N * dt), 20.0)
        s = sigmasq(hp, psd, low_frequency_cutoff=20.0)
        hp /= np.sqrt(s)

        # 1. Pure signal preservation (no background noise)
        target_snr = 25.0
        data_pure = hp * target_snr
        snr_pure = matched_filter(hp, data_pure, psd=psd, low_frequency_cutoff=20.0)
        peak_idx_pure = np.argmax(np.abs(snr_pure))
        rf_pure = get_dynamic_snr_renorm_factor(snr_pure, dt)
        self.assertGreaterEqual(rf_pure[peak_idx_pure], 0.995)

        # 2. Signal in stationary colored noise (signal does not bias local noise variance)
        noise = np.random.randn(N)
        tilde = np.fft.rfft(noise)
        tilde *= np.sqrt(psd.data * 0.5 * sample_rate)
        noise_colored = TimeSeries(np.fft.irfft(tilde, n=N), delta_t=dt, epoch=0.0)

        snr_noise = matched_filter(
            hp, noise_colored, psd=psd, low_frequency_cutoff=20.0
        )
        rf_noise = get_dynamic_snr_renorm_factor(snr_noise, dt)

        data = noise_colored + hp * target_snr
        snr = matched_filter(hp, data, psd=psd, low_frequency_cutoff=20.0)
        peak_idx = np.argmax(np.abs(snr))
        rf_data = get_dynamic_snr_renorm_factor(snr, dt)

        # Signal in noise must match noise background RF to within 0.5%
        self.assertAlmostEqual(rf_data[peak_idx], rf_noise[peak_idx], places=2)

    def test_two_sided_quiet_amplification(self):
        """Test Stage 5 two-sided renorm: quiet detector noise (< 1.0) is boosted up to max_boost_factor."""
        N = 16 * self.sample_rate
        # Local variance 0.25 (std = 0.5) -> unconstrained factor would be 1/sqrt(0.25) = 2.0
        quiet = (0.5 * np.random.randn(N) + 0.5j * np.random.randn(N)).astype(
            np.complex64
        )

        # Baseline variance_floor=1.0 prevents any boost
        rf_baseline = get_dynamic_snr_renorm_factor(quiet, self.dt, variance_floor=1.0)
        self.assertTrue(np.all(rf_baseline <= 1.0001))

        # Two-sided renorm: variance_floor=0.25 with max_boost_factor=1.35
        rf_boosted = get_dynamic_snr_renorm_factor(
            quiet, self.dt, variance_floor=0.25, max_boost_factor=1.35
        )
        # Factor should be strictly capped by safety ceiling 1.35
        self.assertTrue(np.all(rf_boosted <= 1.3501))
        # Interior samples should show boost above 1.0
        interior = rf_boosted[int(4.0 / self.dt) : int(12.0 / self.dt)]
        self.assertGreater(np.mean(interior), 1.25)

    def test_safety_boost_ceiling_clipping(self):
        """Test that safety ceiling strictly clips boost even in zero or near-zero data."""
        zeros = np.zeros(2048, dtype=np.complex64)
        rf_zeros = get_dynamic_snr_renorm_factor(
            zeros, self.dt, variance_floor=0.1, max_boost_factor=1.40
        )
        self.assertTrue(np.all(np.isfinite(rf_zeros)))
        self.assertTrue(np.allclose(rf_zeros, 1.40))

        # Renormalizing zeros must still produce exact zeros
        res_zeros = dynamic_snr_renormalize(
            zeros, self.dt, variance_floor=0.1, max_boost_factor=1.40
        )
        np.testing.assert_array_equal(res_zeros, zeros)

    def test_two_sided_glitch_downweighting_preserved(self):
        """Test that non-stationary glitches are still suppressed when two-sided boost is active."""
        N = 16 * self.sample_rate
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64)
        center = N // 2
        g_start = center + int(1.0 / self.dt)
        g_end = center + int(3.0 / self.dt)
        noise[g_start:g_end] *= 10.0

        rf = get_dynamic_snr_renorm_factor(
            noise, self.dt, variance_floor=0.25, max_boost_factor=1.50
        )
        factor_during_glitch = np.mean(rf[g_start:g_end])
        self.assertLess(factor_during_glitch, 0.25)

    def test_zero_or_negative_variance_floor_robustness(self):
        """Test that variance_floor <= 0 is safely clamped to positive floor without divide-by-zero."""
        zeros = np.zeros(2048, dtype=np.complex64)
        rf_zero_floor = get_dynamic_snr_renorm_factor(
            zeros, self.dt, variance_floor=0.0, max_boost_factor=1.50
        )
        self.assertTrue(np.all(np.isfinite(rf_zero_floor)))
        self.assertTrue(np.all(rf_zero_floor <= 1.50))


class TestFIRChisqAndSearchHelpers(unittest.TestCase):
    def test_validate_dynamic_snr_options(self):
        from pycbc.filter.matched_filter_ratio import validate_dynamic_snr_options

        validate_dynamic_snr_options(1.0, None, disable=False)
        validate_dynamic_snr_options(0.5, 2.0, disable=False)
        validate_dynamic_snr_options(-1.0, -5.0, disable=True)

        with self.assertRaises(ValueError):
            validate_dynamic_snr_options(0.0, None)
        with self.assertRaises(ValueError):
            validate_dynamic_snr_options(1e-8, None)
        with self.assertRaises(ValueError):
            validate_dynamic_snr_options(1.0, 0.0)
        with self.assertRaises(ValueError):
            validate_dynamic_snr_options(1.0, 2e3)

    def test_chisq_filters_indexing_and_prefetch(self):
        from pycbc.vetoes.chisq_fir import ChisqFilters

        rng = np.random.default_rng(123)
        taps = (
            rng.standard_normal((4, 32)) + 1j * rng.standard_normal((4, 32))
        ).astype(np.complex64)
        counts = np.array([16, 24, 32, 20], dtype=np.int64)
        cf_single = ChisqFilters(taps, counts, block_length=256, reach=8)
        cf_batch = ChisqFilters(taps, counts, block_length=256, reach=8)

        with self.assertRaises(TypeError):
            _ = cf_single[1:3]
        with self.assertRaises(TypeError):
            _ = cf_single[[0, 1]]
        with self.assertRaises(IndexError):
            _ = cf_single[4]
        with self.assertRaises(IndexError):
            _ = cf_single[-5]

        cf_batch.prefetch([3, 1, -4, 3])
        for i in range(4):
            np.testing.assert_allclose(cf_single[i], cf_batch[i], rtol=1e-6, atol=1e-6)
        np.testing.assert_array_equal(cf_single[-1], cf_single[3])

    def test_ratio_chisq_chain_and_sg_projection(self):
        import types

        from pycbc.types import FrequencySeries
        from pycbc.vetoes.chain import SignalConsistencyChain
        from pycbc.vetoes.chisq_fir import ChisqFilters, RatioPowerChisq
        from pycbc.vetoes.sgchisq import SingleDetSGChisq

        rng = np.random.default_rng(99)
        sample_rate = 2048.0
        block_len = 512
        n_seg = 4096
        df_seg = sample_rate / n_seg
        n_freq = n_seg // 2 + 1

        psd = FrequencySeries(np.ones(n_freq, dtype=np.float32), delta_f=df_seg)
        psd._cache_key = ("H1", 0)
        stilde_arr = (
            rng.standard_normal(n_freq) + 1j * rng.standard_normal(n_freq)
        ).astype(np.complex64)
        stilde = FrequencySeries(stilde_arr, delta_f=df_seg)
        stilde.psd = psd
        stilde._is_whitened = True
        stilde._cache_key = ("H1", 0)
        stilde.analyze = slice(256, n_seg - 256)
        stilde.cumulative_index = 0

        top_arr = np.ones(n_freq, dtype=np.complex64)
        top_template = FrequencySeries(top_arr, delta_f=df_seg)
        top_template.f_lower = 20.0

        taps = np.zeros((2, 32), dtype=np.complex64)
        taps[:, 16] = 1.0 + 0.0j
        counts = np.array([32, 32], dtype=np.int64)
        filters_f = ChisqFilters(taps, counts, block_length=block_len, reach=10)

        ref_snr = (rng.standard_normal(n_seg) + 1j * rng.standard_normal(n_seg)).astype(
            np.complex64
        )
        ref_snr[1024] = 12.0 + 0.0j
        ref_snr[2048] = 9.0 + 0.0j

        pchisq = RatioPowerChisq(
            num_bins=8,
            sample_rate=sample_rate,
            flow=20.0,
            f_final=800.0,
            top_template=top_template,
            upper_f_by_coarse={},
            snr_threshold=5.0,
        )
        args = types.SimpleNamespace(
            chisq_bins=8,
            chisq_snr_threshold=5.0,
            chisq_window=0.1,
            chisq_iterations=10,
            autochi_stride=2,
            autochi_number_points=4,
            autochi_onesided=None,
            autochi_two_phase=False,
            autochi_reverse_template=False,
            autochi_max_valued=False,
            autochi_max_valued_dof=None,
            sgchisq_snr_threshold=5.0,
            sgchisq_locations=[],
            bank_veto_bank_file=None,
            bank_veto_time_freq=False,
            signal_consistency_chain=None,
        )

        bank_stub = types.SimpleNamespace(
            table={
                "template_hash": np.array([101, 102], dtype=np.int64),
                "f_lower": np.array([20.0, 20.0], dtype=np.float32),
                0: {"mass1": 30.0, "mass2": 30.0, "spin1z": 0.0, "spin2z": 0.0},
                1: {"mass1": 25.0, "mass2": 25.0, "spin1z": 0.0, "spin2z": 0.0},
            }
        )
        pchisq.bank = bank_stub
        sg = SingleDetSGChisq(
            bank_stub,
            num_bins=8,
            snr_threshold=5.0,
            chisq_locations=[],
        )
        sg.do = True
        sg.params = {101: "20-0.5,20-0.6", 102: "20-0.5,20-0.6"}

        chain = SignalConsistencyChain.from_cli(
            args, bank_stub, sample_rate, pchisq, sg
        )

        ids_local = np.array([0, 1], dtype=np.int32)
        t_idxs = np.array([1024, 2048], dtype=np.int64)
        snrs = np.array([12.0 + 0.0j, 9.0 + 0.0j], dtype=np.complex64)
        fine_indices = np.array([0, 1], dtype=np.int32)

        pchisq.set_batch(
            stilde, 0, ids_local, fine_indices[ids_local], t_idxs, ref_snr, filters_f
        )
        kept, vals = chain.values(
            snrs,
            t_idxs,
            fine_indices[ids_local],
            pchisq.template,
            pchisq.correlation,
            stilde=stilde,
        )
        pchisq.ref_f_cache.clear()
        self.assertTrue(np.all(kept))
        self.assertEqual(len(vals["chisq"]), 2)
        self.assertEqual(len(vals["chisq_dof"]), 2)
        self.assertIn("cont_chisq", vals)
        self.assertIn("sg_chisq", vals)
        self.assertTrue(np.all(np.isfinite(vals["chisq"])))
        self.assertTrue(np.all(np.isfinite(vals["cont_chisq"])))
        self.assertEqual(len(pchisq.ref_f_cache), 0)
        self.assertTrue(np.all(np.isfinite(vals["sg_chisq"])))
        self.assertTrue(np.all(vals["sg_chisq"] > 0))

    def test_chain_defaults_and_subthreshold_chisq(self):
        import types

        from pycbc.events import ranking
        from pycbc.vetoes.chain import SignalConsistencyChain
        from pycbc.vetoes.chisq_fir import RatioPowerChisq

        pchisq = RatioPowerChisq(num_bins=64, sample_rate=2048.0, snr_threshold=5.25)
        args = types.SimpleNamespace(
            chisq_bins=64,
            chisq_snr_threshold=5.25,
            chisq_window=1.0,
            chisq_iterations=1,
            autochi_stride=0,
            autochi_number_points=0,
            autochi_onesided=None,
            autochi_two_phase=False,
            autochi_reverse_template=False,
            autochi_max_valued=False,
            autochi_max_valued_dof=None,
            sgchisq_snr_threshold=0,
            sgchisq_locations=None,
            bank_veto_bank_file=None,
            bank_veto_time_freq=False,
            signal_consistency_chain=None,
        )
        bank_stub = types.SimpleNamespace(table={})
        chain = SignalConsistencyChain.from_cli(
            args,
            bank_stub,
            2048.0,
            pchisq,
            defaults={"chisq": 126, "chisq_dof": 126, "sg_chisq": 1.0},
        )
        # Sub-threshold trigger (SNR = 4.2 < 5.25) should receive nominal rchisq = 1.0
        kept, vals = chain.values(
            np.array([4.2 + 0.0j], dtype=np.complex64),
            np.array([1024], dtype=np.int64),
            np.array([0], dtype=np.int32),
            lambda j: None,
            lambda j, reach=0: (None, None),
        )
        self.assertTrue(np.all(kept))
        self.assertEqual(vals["chisq"][0], 126.0)
        self.assertEqual(vals["chisq_dof"][0], 126)
        # Verify EventManager conversion (dof / 2 + 1) -> get_newsnr yields rchisq = 1.0
        em_trigs = {
            "snr": np.array([4.2], dtype=np.float32),
            "chisq": vals["chisq"].astype(np.float32),
            "chisq_dof": (vals["chisq_dof"] / 2 + 1).astype(np.int32),
        }
        np.testing.assert_allclose(ranking.get_newsnr(em_trigs), [4.2], rtol=1e-6)

    def test_ratio_filter_bank_handle_and_group_count(self):
        import os
        import tempfile

        import h5py

        from pycbc.waveform.bank import RatioFilterBank

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_ratio_bank.hdf")
            with h5py.File(path, "w") as f:
                f.attrs["parameters"] = np.array(
                    ["mass1", "mass2", "f_lower", "approximant"], dtype="S"
                )
                f.create_dataset(
                    "mass1", data=np.array([30.0, 29.0, 28.0], dtype=np.float32)
                )
                f.create_dataset(
                    "mass2", data=np.array([30.0, 29.0, 28.0], dtype=np.float32)
                )
                f.create_dataset(
                    "template_hash", data=np.array([1, 2, 3], dtype=np.int64)
                )
                f.create_dataset("approximant", data=np.array([b"IMRPhenomD"] * 3))
                f.create_dataset(
                    "f_lower", data=np.array([20.0, 20.0, 20.0], dtype=np.float32)
                )

                fd = f.create_group("fir_data")
                fd.attrs["num_levels"] = 1
                fd.attrs["sample_rate"] = 2048.0
                lvl1 = fd.create_group("level_1")
                grp0 = lvl1.create_group("0")
                grp0.create_dataset("taps", data=np.ones((3, 16), dtype=np.complex64))
                grp0.create_dataset(
                    "actual_tap_count", data=np.array([12, 8, 12], dtype=np.int32)
                )
                grp0.create_dataset(
                    "fine_bank_index", data=np.array([0, 1, 2], dtype=np.int32)
                )
                grp0.create_dataset(
                    "sample_shift", data=np.array([0, 0, 0], dtype=np.int32)
                )
                grp0.create_dataset("sigmas", data=np.ones((3, 3), dtype=np.float32))

            bank = RatioFilterBank(
                path,
                filter_length=2049,
                delta_f=0.5,
                dtype=np.complex64,
                load_compressed=False,
            )
            self.assertEqual(bank.group_template_count(0, 1), 3)
            bank.close_filehandler()
            self.assertIsNone(bank.filehandler)
            # Lazy reopen on access
            self.assertEqual(bank.group_template_count(0, 1), 3)
            taps, counts, indices = bank.get_firs(0, 1)
            # Equal tap counts (indices 0 and 2 both have 12 taps) preserve relative order (stable sort)
            np.testing.assert_array_equal(indices, np.array([1, 0, 2], dtype=np.int32))
            bank.close_filehandler()


if __name__ == "__main__":
    unittest.main()
