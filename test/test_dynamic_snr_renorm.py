import os
import unittest
import numpy as np

from pycbc.types import Array, TimeSeries
from pycbc.filter.matched_filter_ratio import dynamic_snr_renormalize, get_dynamic_snr_renorm_factor
import pycbc.filter


class TestDynamicSNRRenormalize(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.dt = 1.0 / 2048.0
        self.sample_rate = 2048

    def test_export_availability(self):
        """Verify dynamic_snr_renormalize and get_dynamic_snr_renorm_factor are exported from pycbc.filter."""
        self.assertTrue(hasattr(pycbc.filter, 'dynamic_snr_renormalize'))
        self.assertIs(pycbc.filter.dynamic_snr_renormalize, dynamic_snr_renormalize)
        self.assertTrue(hasattr(pycbc.filter, 'get_dynamic_snr_renorm_factor'))
        self.assertIs(pycbc.filter.get_dynamic_snr_renorm_factor, get_dynamic_snr_renorm_factor)

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
            get_dynamic_snr_renorm_factor(arr, dt=self.dt, window_duration=1.0, hollow_duration=0.5)
        with self.assertRaises(ValueError):
            get_dynamic_snr_renorm_factor(arr, dt=self.dt, window_duration=0.8, hollow_duration=0.5)

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
        renorm = dynamic_snr_renormalize(noise, self.dt, window_duration=8.0, hollow_duration=0.5)

        factor = np.abs(renorm) / np.abs(noise)
        self.assertTrue(np.all(factor <= 1.0001))
        self.assertGreater(np.mean(factor), 0.97)

    def test_stationary_gaussian_noise_with_scale(self):
        """Test that an internally scaled reference buffer with scale=dt produces renorm factor close to 1.0."""
        N = 16 * self.sample_rate
        noise = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64) * self.dt
        renorm = dynamic_snr_renormalize(noise, self.dt, window_duration=8.0, hollow_duration=0.5, scale=self.dt)

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
        noise[center - half_width:center + half_width] += (60.0 + 0j)

        renorm = dynamic_snr_renormalize(noise, self.dt, window_duration=8.0, hollow_duration=0.5)
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

        renorm = dynamic_snr_renormalize(noise, self.dt, window_duration=8.0, hollow_duration=0.5)
        factor_during_glitch = np.mean(np.abs(renorm[g_start:g_end]) / np.abs(noise[g_start:g_end]))
        self.assertLess(factor_during_glitch, 0.25)

    def test_variance_floor_no_amplification(self):
        """Test that quiet periods or zero arrays are never amplified (variance floor >= 1.0)."""
        zeros = np.zeros(2048, dtype=np.complex64)
        renorm = dynamic_snr_renormalize(zeros, self.dt, variance_floor=1.0)
        np.testing.assert_array_equal(renorm, zeros)

        # Very low amplitude noise
        quiet = (np.random.randn(2048) + 1j * np.random.randn(2048)).astype(np.complex64) * 0.01
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
        from pycbc.waveform import get_td_waveform
        from pycbc.psd import aLIGOZeroDetHighPower
        from pycbc.filter import matched_filter, sigmasq

        sample_rate = 2048
        dt = 1.0 / sample_rate
        N = 32 * sample_rate
        hp, _ = get_td_waveform(approximant='IMRPhenomD', mass1=30, mass2=30,
                                f_lower=20.0, delta_t=dt)
        hp.resize(N)
        hp.start_time = 0.0

        psd = aLIGOZeroDetHighPower(N//2 + 1, 1.0 / (N * dt), 20.0)
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

        snr_noise = matched_filter(hp, noise_colored, psd=psd, low_frequency_cutoff=20.0)
        rf_noise = get_dynamic_snr_renorm_factor(snr_noise, dt)

        data = noise_colored + hp * target_snr
        snr = matched_filter(hp, data, psd=psd, low_frequency_cutoff=20.0)
        peak_idx = np.argmax(np.abs(snr))
        rf_data = get_dynamic_snr_renorm_factor(snr, dt)

        # Signal in noise must match noise background RF to within 0.5%
        self.assertAlmostEqual(rf_data[peak_idx], rf_noise[peak_idx], places=2)

    def test_disable_env_var(self):
        """Test that setting PYCBC_DISABLE_DYNAMIC_SNR_RENORM=1 bypasses renormalization."""
        N = 4096
        arr = (np.random.randn(N) + 1j * np.random.randn(N)).astype(np.complex64) * 10.0
        os.environ['PYCBC_DISABLE_DYNAMIC_SNR_RENORM'] = '1'
        try:
            res = dynamic_snr_renormalize(arr, self.dt)
            np.testing.assert_array_equal(res, arr)
        finally:
            del os.environ['PYCBC_DISABLE_DYNAMIC_SNR_RENORM']

    def test_two_sided_quiet_amplification(self):
        """Test Stage 5 two-sided renorm: quiet detector noise (< 1.0) is boosted up to max_boost_factor."""
        N = 16 * self.sample_rate
        # Local variance 0.25 (std = 0.5) -> unconstrained factor would be 1/sqrt(0.25) = 2.0
        quiet = (0.5 * np.random.randn(N) + 0.5j * np.random.randn(N)).astype(np.complex64)

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
        interior = rf_boosted[int(4.0 / self.dt):int(12.0 / self.dt)]
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


if __name__ == '__main__':
    unittest.main()
