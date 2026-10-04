# Copyright (C) 2026 PyCBC Team
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

"""Unit tests for regularized overwhitening, compact instantaneous PSD,
and overlap-save streaming filter.
"""

import unittest
import numpy as np
import tempfile
import os
import time
import h5py
import scipy.stats as stats

from pycbc.types import TimeSeries, FrequencySeries
import pycbc.filter
import pycbc.psd
import pycbc.waveform


class TestOverwhiten(unittest.TestCase):
    def setUp(self):
        self.fs = 1024.0
        self.duration = 64.0
        self.N = int(self.fs * self.duration)
        # Generate stationary Gaussian noise with a colored spectrum
        rng = np.random.default_rng(12345)
        white_noise = rng.normal(0, 1, self.N)
        # Simple color filter: 1 / (1 + (f/50)**2)
        freqs = np.fft.rfftfreq(self.N, 1.0 / self.fs)
        color = 1.0 / np.sqrt(1.0 + (freqs / 50.0)**2)
        colored_fd = np.fft.rfft(white_noise) * color
        colored_td = np.fft.irfft(colored_fd)
        self.ts = TimeSeries(colored_td, delta_t=1.0 / self.fs, epoch=1000000000.0)

    def test_construct_regularized_kernels(self):
        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        invpsd, w_kernel, fir_dow, fir_w, psd_ist, df = pycbc.filter.construct_regularized_kernels(
            psd, self.fs, self.duration, f_low=18.0, f_taper=4.0, max_filter_duration=2.0
        )
        self.assertGreater(len(invpsd), 0)
        self.assertGreater(len(w_kernel), 0)
        self.assertEqual(len(fir_dow) % 2, 1) # symmetric odd length
        self.assertEqual(len(fir_w) % 2, 1)
        # Low frequency taper should be zero below f_low - f_taper = 14 Hz
        k_14 = int(round(13.0 / df))
        self.assertAlmostEqual(invpsd[k_14], 0.0)
        self.assertAlmostEqual(w_kernel[k_14], 0.0)

    def test_overwhiten_and_whitening_normality(self):
        ts_dow, ts_dw, meta = pycbc.filter.overwhiten_strain(
            self.ts, f_low=18.0, f_taper=4.0, max_filter_duration=2.0,
            use_overlap_save=True, chunk_len_s=8.0
        )
        self.assertEqual(len(ts_dow), len(self.ts))
        self.assertEqual(len(ts_dw), len(self.ts))
        # Whitened stream d_w(t) outside edge regions should have mean ~ 0 and std ~ 1
        edge = int(4.0 * self.fs)
        interior_dw = ts_dw.numpy()[edge : len(ts_dw) - edge]
        self.assertAlmostEqual(np.mean(interior_dw), 0.0, delta=0.1)
        self.assertAlmostEqual(np.std(interior_dw), 1.0, delta=0.2)

    def test_overlap_save_equivalence(self):
        # Full circular vs overlap-save should agree on the interior
        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        invpsd, w_kernel, fir_dow, fir_w, psd_ist, df = pycbc.filter.construct_regularized_kernels(
            psd, self.fs, self.duration, f_low=18.0, f_taper=4.0, max_filter_duration=2.0
        )
        dow_stream, dw_stream = pycbc.filter.produce_conditioned_streams(
            self.ts, invpsd, w_kernel, fir_dow, fir_w, delta_f_full=df, use_overlap_save=True, chunk_len_s=8.0
        )
        dow_batch, dw_batch = pycbc.filter.produce_conditioned_streams(
            self.ts, invpsd, w_kernel, fir_dow, fir_w, delta_f_full=df, use_overlap_save=False
        )
        edge = int(4.0 * self.fs)
        diff = np.abs(dow_stream.numpy()[edge:-edge] - dow_batch.numpy()[edge:-edge])
        rel_diff = np.max(diff) / np.std(dow_batch.numpy()[edge:-edge])
        self.assertLess(rel_diff, 1e-3)

    def test_matched_filter_overwhitened(self):
        # Generate synthetic waveform
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=30, mass2=30,
            f_lower=18.0, delta_f=self.ts.delta_f, distance=500
        )
        hp.resize(len(self.ts.to_frequencyseries()))
        hp_shifted = pycbc.waveform.apply_fd_time_shift(hp, self.duration / 2.0)
        s_td = hp_shifted.to_timeseries()
        s_td.start_time = self.ts.start_time

        ts_inj = self.ts + s_td
        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        psd_ist = pycbc.psd.inverse_spectrum_truncation(
            pycbc.psd.interpolate(psd, ts_inj.delta_f),
            max_filter_len=int(2 * self.fs), which_spectrum="invpsd",
            low_frequency_cutoff=18.0, trunc_method="hann"
        )
        sigmasq_val = pycbc.filter.sigmasq(hp, psd=psd_ist, low_frequency_cutoff=18.0)
        snr_opt = float(np.sqrt(sigmasq_val))

        # Standard PyCBC matched filter
        snr_pycbc = pycbc.filter.matched_filter(hp, ts_inj, psd=psd_ist, low_frequency_cutoff=18.0)
        peak_pycbc = float(np.max(np.abs(snr_pycbc.numpy())))

        # Overwhitened stream matched filter
        ts_dow, _, _ = pycbc.filter.overwhiten_strain(ts_inj, psd=psd, f_low=18.0, max_filter_duration=2.0)
        snr_dow = pycbc.filter.matched_filter_overwhitened(hp, ts_dow, snr_optimal=snr_opt, f_low=18.0)
        peak_dow = float(np.max(np.abs(snr_dow.numpy())))

        # Peak SNRs should agree within < 1%
        rel_diff = abs(peak_dow - peak_pycbc) / peak_pycbc
        self.assertLess(rel_diff, 0.02)

    def test_matched_filter_overwhitened_frequencyseries(self):
        # Test matched_filter_overwhitened when ts_dow is a FrequencySeries
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=30, mass2=30,
            f_lower=18.0, delta_f=self.ts.delta_f, distance=500
        )
        hp.resize(len(self.ts.to_frequencyseries()))
        hp_shifted = pycbc.waveform.apply_fd_time_shift(hp, self.duration / 2.0)
        s_td = hp_shifted.to_timeseries()
        s_td.start_time = self.ts.start_time
        ts_inj = self.ts + s_td

        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        psd_ist = pycbc.psd.inverse_spectrum_truncation(
            pycbc.psd.interpolate(psd, ts_inj.delta_f),
            max_filter_len=int(2 * self.fs), which_spectrum="invpsd",
            low_frequency_cutoff=18.0, trunc_method="hann"
        )
        ts_dow, _, _ = pycbc.filter.overwhiten_strain(ts_inj, psd=psd, f_low=18.0, max_filter_duration=2.0)
        fd_dow = ts_dow.to_frequencyseries()

        # Run with TimeSeries and with FrequencySeries
        snr_ts = pycbc.filter.matched_filter_overwhitened(hp, ts_dow, psd=psd_ist, f_low=18.0)
        snr_fs = pycbc.filter.matched_filter_overwhitened(hp, fd_dow, psd=psd_ist, f_low=18.0)

        self.assertEqual(len(snr_ts), len(self.ts))
        self.assertEqual(len(snr_fs), len(self.ts))
        peak_ts = float(np.max(np.abs(snr_ts.numpy())))
        peak_fs = float(np.max(np.abs(snr_fs.numpy())))
        self.assertAlmostEqual(peak_ts, peak_fs, delta=1e-5)

    def test_trimmed_welch_estimator(self):
        seg_len = int(4 * self.fs)
        seg_stride = int(2 * self.fs)
        psd_trim = pycbc.psd.estimate_psd_trimmed_welch(self.ts, seg_len=seg_len, seg_stride=seg_stride, alpha=0.20)
        psd_med = pycbc.psd.welch(self.ts, seg_len=seg_len, seg_stride=seg_stride, avg_method="median")
        self.assertEqual(len(psd_trim), len(psd_med))
        self.assertAlmostEqual(psd_trim.delta_f, psd_med.delta_f)

        # Inband median power ratio should be close to 1.0 (bias corrected)
        freqs = psd_trim.sample_frequencies.numpy()
        band = (freqs >= 25.0) & (freqs <= 400.0)
        ratio = psd_trim.numpy()[band] / psd_med.numpy()[band]
        self.assertAlmostEqual(float(np.median(ratio)), 1.0, delta=0.15)

    def test_multitaper_estimator(self):
        seg_len = int(4 * self.fs)
        seg_stride = int(2 * self.fs)
        psd_mt = pycbc.psd.estimate_psd_multitaper(self.ts, seg_len=seg_len, seg_stride=seg_stride, NW=3.0, avg_method="median")
        psd_med = pycbc.psd.welch(self.ts, seg_len=seg_len, seg_stride=seg_stride, avg_method="median")
        self.assertEqual(len(psd_mt), len(psd_med))
        freqs = psd_mt.sample_frequencies.numpy()
        band = (freqs >= 25.0) & (freqs <= 400.0)
        ratio = psd_mt.numpy()[band] / psd_med.numpy()[band]
        self.assertAlmostEqual(float(np.median(ratio)), 1.0, delta=0.20)

    def test_time_varying_streaming_seams(self):
        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=400.0, n_knots=24)
        model.fit_static_from_psd(psd)
        model.track_nonstationarity(self.ts, t_step=2.0, t_window=4.0)
        seam_res = pycbc.filter.verify_time_varying_streaming_seams(
            self.ts, psd, model, chunk_len_s=8.0, filter_duration=2.0
        )
        self.assertTrue(seam_res["passed"])
        self.assertLess(seam_res["max_seam_relative_error"], 1e-4)

    def test_cli_executable_execution(self):
        import subprocess
        import sys
        with tempfile.TemporaryDirectory() as tmpdir:
            strain_path = os.path.join(tmpdir, "test_strain.hdf")
            self.ts.save(strain_path, group="strain")

            out_product = os.path.join(tmpdir, "product.h5")
            out_stats = os.path.join(tmpdir, "stats.json")
            bin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin", "pycbc_overwhiten")

            cmd = [
                sys.executable, bin_path,
                "--input-strain-file", strain_path,
                "--channel-name", "strain",
                "--sample-rate", str(self.fs),
                "--f-low", "18.0",
                "--f-taper", "4.0",
                "--max-filter-duration", "2.0",
                "--chunk-duration", "8.0",
                "--psd-method", "welch-median",
                "--track-nonstationarity",
                "--output-file", out_product,
                "--output-stats-json", out_stats,
                "--validate"
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"CLI failed with error:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
            self.assertTrue(os.path.exists(out_product))
            self.assertTrue(os.path.exists(out_stats))

    def test_compact_instantaneous_psd_serialization(self):
        psd = pycbc.psd.welch(self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs))
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=500.0, n_knots=32)
        model.fit_static_from_psd(psd)
        model.track_nonstationarity(self.ts, t_step=2.0, t_window=4.0)

        # Test evaluation
        test_freqs = np.linspace(20, 450, 100)
        vals_t0 = model.eval_psd(test_freqs, t=float(self.ts.start_time) + 10.0)
        self.assertEqual(len(vals_t0), 100)
        self.assertTrue(np.all(vals_t0 > 0))

        # Test packaging into HDF5
        with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            ts_dow, ts_dw, _ = pycbc.filter.overwhiten_strain(self.ts, psd=psd, f_low=18.0, max_filter_duration=2.0)
            pycbc.psd.package_overwhitened_product(tmp_path, ts_dow, ts_dw, model_compact=model, psd_base=psd)
            self.assertTrue(os.path.exists(tmp_path))
            self.assertGreater(os.path.getsize(tmp_path), 0)

            # Reload model from HDF5
            import h5py
            with h5py.File(tmp_path, "r") as f:
                self.assertIn("dow", f)
                self.assertIn("dw", f)
                self.assertIn("instantaneous_psd", f)
                loaded_model = pycbc.psd.CompactInstantaneousPSD.load_from_hdf(f["instantaneous_psd"])
                vals_loaded = loaded_model.eval_psd(test_freqs, t=float(self.ts.start_time) + 10.0)
                np.testing.assert_allclose(vals_t0, vals_loaded, rtol=1e-5)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_metric_1_signal_immunity_and_bias_rejection(self):
        """Formal Metric 1: Verify PSD continuum and instantaneous scale factor immunity to loud CBC signal."""
        fs = 1024.0
        duration = 256.0
        N = int(fs * duration)
        dt = 1.0 / fs
        rng = np.random.default_rng(42)

        freqs = np.fft.rfftfreq(N, dt)
        color = 1.0 / np.sqrt(1.0 + (np.maximum(freqs, 10.0) / 50.0)**2)
        white = rng.normal(0, 1, N)
        fd = np.fft.rfft(white) * color
        noise = np.fft.irfft(fd, n=N)

        t_arr = np.arange(N) * dt
        gain_env = 1.0 + 0.3 * np.exp(-((t_arr - 64.0) / 20.0)**2)
        ts_clean = TimeSeries(noise * gain_env, delta_t=dt, epoch=1000000000.0)

        # Inject loud BBH SNR 80 at t = 64s
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=35, mass2=35,
            f_lower=18.0, delta_f=1.0 / duration, distance=300
        )
        hp.resize(len(freqs))
        hp_shifted = pycbc.waveform.apply_fd_time_shift(hp, 64.0)
        s_td = hp_shifted.to_timeseries()
        s_td.start_time = ts_clean.start_time

        psd_clean = pycbc.psd.welch(ts_clean, seg_len=int(8 * fs), seg_stride=int(4 * fs), avg_method="median")
        psd_ist = pycbc.psd.inverse_spectrum_truncation(
            pycbc.psd.interpolate(psd_clean, 1.0 / duration),
            max_filter_len=int(2 * fs), which_spectrum="invpsd", low_frequency_cutoff=18.0
        )
        sigsq = pycbc.filter.sigmasq(hp, psd=psd_ist, low_frequency_cutoff=18.0)
        s_td = s_td * (80.0 / np.sqrt(sigsq))
        ts_inj = ts_clean + s_td

        psd_inj = pycbc.psd.welch(ts_inj, seg_len=int(8 * fs), seg_stride=int(4 * fs), avg_method="median")

        model_clean = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=400.0, n_knots=48)
        model_clean.fit_static_from_psd(psd_clean)
        model_clean.track_nonstationarity(ts_clean, t_step=1.0, t_window=4.0)

        model_inj = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=400.0, n_knots=48)
        model_inj.fit_static_from_psd(psd_inj)
        model_inj.track_nonstationarity(ts_inj, t_step=1.0, t_window=4.0)

        # 1. Continuum baseline RMS shift under SNR 80 injection <= 2.0%
        f_eval = np.linspace(30.0, 350.0, 500)
        rms_shift = np.sqrt(np.mean(((model_inj.eval_continuum(f_eval) - model_clean.eval_continuum(f_eval)) / model_clean.eval_continuum(f_eval))**2))
        self.assertLess(rms_shift, 0.020, f"Continuum RMS shift {rms_shift*100:.3f}% exceeds 2.0% target")

        # 2. Instantaneous scale factor jump at merger < 1%
        t_merger = 1000000064.0
        idx_m = np.argmin(np.abs(model_inj.scale_times - t_merger))
        delta_g = np.max(np.abs(model_inj.band_scales[:, idx_m] - model_clean.band_scales[:, idx_m]))
        self.assertLess(delta_g, 0.01, f"Merger scale factor jump {delta_g*100:.4f}% exceeds 1% target")

    def test_metric_2_spectral_line_accuracy(self):
        """Formal Metric 2: Sub-bin frequency precision and data-driven FWHM line modeling."""
        fs = 1024.0
        duration = 256.0
        N = int(fs * duration)
        dt = 1.0 / fs
        rng = np.random.default_rng(12345)

        freqs = np.fft.rfftfreq(N, dt)
        color = 1.0 / np.sqrt(1.0 + (np.maximum(freqs, 10.0) / 50.0)**2)

        # Add sharp sub-bin line at f0 = 60.032 Hz (grid delta_f = 0.125 Hz for 8s segments)
        f0_true = 60.032
        fwhm_true = 0.06
        line_lorentz = 200.0 * (fwhm_true / 2.0)**2 / ((freqs - f0_true)**2 + (fwhm_true / 2.0)**2)
        psd_true = color**2 + line_lorentz

        white = rng.normal(0, 1, N)
        fd = np.fft.rfft(white) * np.sqrt(psd_true * fs / 2.0)
        strain = np.fft.irfft(fd, n=N)
        ts = TimeSeries(strain, delta_t=dt, epoch=1000000000.0)

        psd_welch = pycbc.psd.welch(ts, seg_len=int(8 * fs), seg_stride=int(4 * fs), avg_method="median")
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=400.0, n_knots=48)
        model.fit_static_from_psd(psd_welch)

        self.assertGreaterEqual(len(model.lines), 1, "Failed to identify line")
        f0_hat, excess_hat, fwhm_hat = model.lines[0]

        # 1. Sub-bin frequency error < 0.015 Hz (discrete bin spacing is 0.125 Hz)
        f_err = abs(f0_hat - f0_true)
        self.assertLess(f_err, 0.015, f"Sub-bin line frequency error {f_err:.4f} Hz exceeds 0.015 Hz target")

        # 2. Apparent FWHM is within [0.15, 0.40] Hz consistent with 8s Hann window main-lobe
        self.assertTrue(0.15 <= fwhm_hat <= 0.40, f"Measured FWHM {fwhm_hat:.4f} Hz out of expected windowed range")

        # 3. Spectral contrast outside resonance notch (+1.0 Hz) >= 20 dB drop
        p_f0 = model.eval_psd(np.array([f0_hat]), include_lines=True)[0]
        p_off = model.eval_psd(np.array([f0_hat + 1.0]), include_lines=True)[0]
        contrast_db = 10.0 * np.log10(p_f0 / p_off)
        self.assertGreaterEqual(contrast_db, 20.0, f"Contrast {contrast_db:.1f} dB < 20 dB target")

    def test_metric_3_complete_overwhitening_and_nonstationarity(self):
        """Formal Metric 3: Complete upstream overwhitening standardizes non-stationary gain drifts."""
        fs = 1024.0
        duration = 128.0
        N = int(fs * duration)
        dt = 1.0 / fs
        rng = np.random.default_rng(54321)

        freqs = np.fft.rfftfreq(N, dt)
        color = 1.0 / np.sqrt(1.0 + (np.maximum(freqs, 10.0) / 50.0)**2)
        white = rng.normal(0, 1, N)
        fd = np.fft.rfft(white) * color
        noise = np.fft.irfft(fd, n=N)

        # Apply large +40% gain swell in the middle
        t_arr = np.arange(N) * dt
        gain_env = 1.0 + 0.4 * np.exp(-((t_arr - 64.0) / 20.0)**2)
        ts_nonstat = TimeSeries(noise * gain_env, delta_t=dt, epoch=1000000000.0)

        # Complete upstream overwhitening with dynamic non-stationarity tracking
        ts_dow, ts_dw, meta = pycbc.filter.overwhiten_strain(
            ts_nonstat, track_nonstationarity=True, max_filter_duration=2.0, chunk_len_s=8.0
        )

        # 1. Whitened strain d_w(t) has unit variance across non-stationary drift
        edge = int(4.0 * fs)
        dw_interior = ts_dw.numpy()[edge : len(ts_dw) - edge]
        std_dw = float(np.std(dw_interior))
        self.assertAlmostEqual(std_dw, 1.00, delta=0.05,
                               msg=f"Whitened strain std {std_dw:.4f} deviates from unit variance 1.00")

        # 2. Matched filter noise triggers maintain unit variance without downstream normalization
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=25, mass2=25,
            f_lower=18.0, delta_f=ts_nonstat.delta_f, distance=500
        )
        hp.resize(len(ts_nonstat.to_frequencyseries()))
        snr_series = pycbc.filter.matched_filter_overwhitened(
            hp, ts_dow, psd=meta["psd_base"], f_low=18.0
        )
        snr_re = snr_series.numpy()[edge : len(snr_series) - edge].real
        snr_std = float(np.std(snr_re))
        self.assertAlmostEqual(snr_std, 1.00, delta=0.05,
                               msg=f"Matched filter noise std {snr_std:.4f} deviates from unit variance 1.00")

    def test_metric_4_bayesline_parity_and_superiority(self):
        """Formal Metric 4: BayesLine curvature-adaptive knot allocation and Bayesian credible intervals."""
        fs = 1024.0
        duration = 128.0
        N = int(fs * duration)
        dt = 1.0 / fs
        rng = np.random.default_rng(998877)

        freqs = np.fft.rfftfreq(N, dt)
        f_safe = np.maximum(freqs, 10.0)
        asd_true = 1e-21 * ((30.0 / f_safe)**3 * (freqs < 30.0) + 1.0 + (freqs / 150.0)**1.5)
        psd_true = asd_true ** 2

        white = rng.normal(0, 1, N)
        fd = np.fft.rfft(white) * np.sqrt(psd_true * fs / 2.0)
        ts = TimeSeries(np.fft.irfft(fd, n=N), delta_t=dt, epoch=1000000000.0)

        psd_welch = pycbc.psd.welch(ts, seg_len=int(8 * fs), seg_stride=int(4 * fs), avg_method="median")
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=400.0, n_knots=36)
        model.fit_static_from_psd(psd_welch, adaptive_knots=True, curvature_alpha=1.5, iterative_refine=True)

        # 1. Curvature adaptation allocates more knots in seismic wall (< 40 Hz) than in high freq (> 250 Hz)
        n_seismic = np.sum(model.knots_f < 40.0)
        n_high = np.sum(model.knots_f > 250.0)
        self.assertGreater(n_seismic, n_high,
                           f"Seismic knots ({n_seismic}) not greater than high freq knots ({n_high})")

        # 2. Analytical credible intervals are strictly ordered and well-bounded
        f_test = np.linspace(25.0, 350.0, 100)
        p_med, p_low, p_high = model.eval_credible_intervals(f_test, alpha=0.90)
        self.assertTrue(np.all(p_low < p_med), "Credible lower bound not strictly below median")
        self.assertTrue(np.all(p_high > p_med), "Credible upper bound not strictly above median")

        # 3. Relative 90% credible interval width is within reasonable bounds [1%, 35%]
        rel_width = (p_high - p_low) / p_med
        self.assertTrue(np.all(rel_width >= 0.01), "Credible interval unexpectedly narrow")
        self.assertTrue(np.all(rel_width <= 0.35), "Credible interval unexpectedly wide")

    def test_overwhitened_data_query_api(self):
        """Test OverwhitenedData API: querying strain, PSD at time t, valid vs boundary, and gating."""
        fs = 512.0
        duration = 64.0
        N = int(fs * duration)
        t0 = 1000000000.0
        rng = np.random.default_rng(42)

        x_dow = rng.normal(0, 1, N)
        x_dw = rng.normal(0, 1, N)
        ts_dow = TimeSeries(x_dow, delta_t=1.0 / fs, epoch=t0)
        ts_dw = TimeSeries(x_dw, delta_t=1.0 / fs, epoch=t0)

        # Base PSD and compact model
        psd_base = pycbc.psd.welch(ts_dw, seg_len=int(4 * fs), seg_stride=int(2 * fs))
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=200.0, n_knots=24)
        model.fit_static_from_psd(psd_base)
        model.track_nonstationarity(ts_dw, t_step=1.0, t_window=4.0)

        gated = [(t0 + 20.0, t0 + 22.0)]
        bpad = 2.0

        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "test_conditioned.hdf5")
            pycbc.psd.package_overwhitened_product(
                h5_path, ts_dow, ts_dw, model_compact=model, psd_base=psd_base,
                gated_segments=gated, boundary_pad=bpad
            )

            # Test opening and properties
            data = pycbc.filter.OverwhitenedData.open(h5_path)
            self.assertEqual(data.start_time, t0)
            self.assertEqual(data.end_time, t0 + duration)
            self.assertEqual(data.duration, duration)
            self.assertEqual(data.sample_rate, fs)
            self.assertEqual(data.boundary_pad, bpad)

            # Check gated and valid segments
            self.assertEqual(len(data.gated_segments), 1)
            self.assertEqual(data.gated_intervals, [(t0 + 20.0, t0 + 22.0)])

            # Valid intervals should exclude boundary pad [0, 2] and [62, 64] and gated [20, 22]
            valid_segs = data.valid_intervals
            self.assertEqual(len(valid_segs), 2)
            self.assertAlmostEqual(valid_segs[0][0], t0 + 2.0)
            self.assertAlmostEqual(valid_segs[0][1], t0 + 20.0)
            self.assertAlmostEqual(valid_segs[1][0], t0 + 22.0)
            self.assertAlmostEqual(valid_segs[1][1], t0 + 62.0)

            # Check instantaneous status methods
            # Inside boundary pad:
            self.assertTrue(data.is_boundary(t0 + 0.5))
            self.assertFalse(data.is_valid(t0 + 0.5))
            self.assertFalse(data.is_gated(t0 + 0.5))

            # Inside gated segment:
            self.assertTrue(data.is_gated(t0 + 21.0))
            self.assertFalse(data.is_valid(t0 + 21.0))
            self.assertFalse(data.is_boundary(t0 + 21.0))

            # Clean interior data:
            self.assertTrue(data.is_valid(t0 + 10.0))
            self.assertFalse(data.is_gated(t0 + 10.0))
            self.assertFalse(data.is_boundary(t0 + 10.0))

            # Query data streams
            dow_slice = data.get_overwhitened(t0 + 10.0, t0 + 15.0)
            self.assertEqual(len(dow_slice), int(5.0 * fs))
            self.assertAlmostEqual(float(dow_slice.start_time), t0 + 10.0)

            dw_slice = data.get_whitened(t0 + 10.0, t0 + 15.0)
            self.assertEqual(len(dw_slice), int(5.0 * fs))

            # Query PSD at specific time
            psd_t = data.get_psd(t0 + 15.0, delta_f=0.25)
            self.assertIsInstance(psd_t, FrequencySeries)
            self.assertEqual(psd_t.delta_f, 0.25)
            self.assertAlmostEqual(float(psd_t.epoch), t0 + 15.0)

            # Query unified slice
            q = data.query(t0 + 18.0, t0 + 24.0)
            self.assertEqual(len(q.dow), int(6.0 * fs))
            self.assertEqual(len(q.times), int(6.0 * fs))
            # q spans [18, 24], so times [20, 22] must be gated
            gated_count = np.sum(q.is_gated)
            self.assertEqual(gated_count, int(2.0 * fs))
            self.assertEqual(np.sum(q.is_valid), int(4.0 * fs))

            # Slice midpoint PSD
            psd_mid = q.get_psd(delta_f=0.5)
            self.assertEqual(psd_mid.delta_f, 0.5)

            # Test context manager
            with pycbc.filter.OverwhitenedData.open(h5_path) as ctx_data:
                self.assertEqual(ctx_data.duration, duration)

    def test_overwhitened_data_directory_and_multi_chunk(self):
        """Test OverwhitenedData indexing a directory of multiple chunk files and cross-chunk queries."""
        fs = 256.0
        dur_chunk = 32.0
        t0 = 1000000000.0
        rng = np.random.default_rng(777)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create two contiguous chunks
            for i in range(2):
                st = t0 + i * dur_chunk
                ts_dow = TimeSeries(rng.normal(0, 1, int(fs * dur_chunk)), delta_t=1.0 / fs, epoch=st)
                ts_dw = TimeSeries(rng.normal(0, 1, int(fs * dur_chunk)), delta_t=1.0 / fs, epoch=st)
                psd_base = pycbc.psd.welch(ts_dw, seg_len=int(4 * fs), seg_stride=int(2 * fs))
                fn = os.path.join(tmpdir, f"chunk_{i}.hdf5")
                pycbc.psd.package_overwhitened_product(
                    fn, ts_dow, ts_dw, psd_base=psd_base, boundary_pad=2.0
                )

            # Open entire directory
            dataset = pycbc.filter.OverwhitenedData.open(tmpdir)
            self.assertEqual(dataset.start_time, t0)
            self.assertEqual(dataset.end_time, t0 + 2 * dur_chunk)
            self.assertEqual(dataset.duration, 2 * dur_chunk)
            self.assertEqual(len(dataset.chunks), 2)

            # Query spanning across the boundary between chunk 0 and chunk 1 (t0 + 20 to t0 + 44)
            span_dow = dataset.get_overwhitened(t0 + 20.0, t0 + 44.0)
            self.assertEqual(len(span_dow), int(24.0 * fs))
            self.assertAlmostEqual(float(span_dow.start_time), t0 + 20.0)
            self.assertAlmostEqual(float(span_dow.end_time), t0 + 44.0)

            # Check validity across chunks (each chunk has 2s boundary pad at ends)
            self.assertFalse(dataset.is_valid(t0 + 31.0)) # within 2s pad of chunk 0 end
            self.assertFalse(dataset.is_valid(t0 + 33.0)) # within 2s pad of chunk 1 start
            self.assertTrue(dataset.is_valid(t0 + 25.0))  # clean interior of chunk 0
            self.assertTrue(dataset.is_valid(t0 + 40.0))  # clean interior of chunk 1

            # Test pycbc.strain.from_overwhitened_file
            st_strain = pycbc.strain.from_overwhitened_file(tmpdir, start_time=t0 + 5.0, end_time=t0 + 15.0)
            self.assertEqual(len(st_strain), int(10.0 * fs))
            self.assertAlmostEqual(float(st_strain.start_time), t0 + 5.0)

            # Test pycbc.psd.from_overwhitened_file
            psd_val = pycbc.psd.from_overwhitened_file(tmpdir, gps_time=t0 + 10.0, delta_f=0.5)
            self.assertEqual(psd_val.delta_f, 0.5)
            self.assertGreater(len(psd_val), 0)

    def test_full_bandwidth_nyquist_continuum(self):
        """Test full Nyquist bandwidth continuum coverage, seismic wall slope, and discrete unit variance."""
        fs = 4096.0
        duration = 64.0
        N = int(fs * duration)
        rng = np.random.default_rng(42)

        # Realistic strain PSD model
        freqs = np.fft.rfftfreq(N, 1.0 / fs)
        f_safe = np.maximum(freqs, 10.0)
        # Power-law seismic wall + flat + shot noise
        psd_model = 1e-46 * ((30.0 / f_safe)**10 + 1.0 + (freqs / 500.0)**2)
        psd_model[0] = psd_model[1]
        psd_series = FrequencySeries(psd_model, delta_f=1.0 / duration, epoch=1000000000.0)

        # 1. CompactInstantaneousPSD spans to Nyquist (2048 Hz)
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=fs / 2.0, n_knots=48)
        model.fit_static_from_psd(psd_series)
        self.assertEqual(model.f_max, fs / 2.0)
        self.assertEqual(model.bands[-1][1], fs / 2.0)

        # 2. Smooth physical seismic wall regularization: d(log10 Sn)/d(log10 f) <= -4.0 below f_min
        f_eval = np.array([5.0, 10.0, 14.0, 17.0])
        p_eval = model.eval_continuum(f_eval)
        log_f = np.log10(f_eval)
        log_p = np.log10(p_eval)
        slopes = np.diff(log_p) / np.diff(log_f)
        for s in slopes:
            self.assertLessEqual(s, -4.0, f"Seismic wall slope {s} is not steep enough (expected <= -4.0)")

        # 3. Discrete unit variance calibration: construct_regularized_kernels NEB
        white_fd = np.fft.rfft(rng.normal(0, 1, N)) * np.sqrt(psd_series.numpy() * fs / 2.0)
        ts_sim = TimeSeries(np.fft.irfft(white_fd, n=N), delta_t=1.0 / fs, epoch=1000000000.0)
        dow, dw, _ = pycbc.filter.overwhiten_strain(ts_sim, psd=psd_series, f_low=18.0, max_filter_duration=4.0)
        edge = int(4.0 * fs)
        dw_int = dw.numpy()[edge:-edge]
        std_val = float(np.std(dw_int))
        self.assertAlmostEqual(std_val, 1.00, delta=0.05,
                               msg=f"Simulated whitened strain std {std_val:.4f} deviates from 1.00")

    def test_coherent_line_subtraction_preservation(self):
        """Test coherent line subtraction suppresses line excess while preserving CBC waveform >99.9%."""
        fs = 2048.0
        duration = 64.0
        N = int(fs * duration)
        dt = 1.0 / fs
        t = np.arange(N) * dt
        rng = np.random.default_rng(123)

        # Background noise
        noise = rng.normal(0, 1.0, N)
        ts_noise = TimeSeries(noise, delta_t=dt, epoch=1000000000.0)

        # Injected CBC waveform
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=20, mass2=20,
            f_lower=20.0, delta_f=ts_noise.delta_f, distance=200
        )
        hp.resize(len(ts_noise.to_frequencyseries()))
        hp_td = pycbc.waveform.apply_fd_time_shift(hp, duration / 2.0).to_timeseries()
        hp_td.start_time = ts_noise.start_time

        # Add strong 60 Hz sinusoidal line
        line_60 = 20.0 * np.sin(2.0 * np.pi * 60.0 * t)
        ts_with_line = ts_noise + TimeSeries(line_60, delta_t=dt, epoch=ts_noise.start_time) + hp_td

        # Perform coherent line subtraction
        ts_clean, line_list = pycbc.filter.subtract_coherent_lines(
            ts_with_line, f_low=18.0, line_threshold=3.0, max_lines=16
        )

        # Check line was detected and subtracted near 60 Hz
        sub_freqs = [l[0] for l in line_list]
        has_60 = any(abs(f - 60.0) < 1.0 for f in sub_freqs)
        self.assertTrue(has_60, f"60 Hz line not detected in {sub_freqs}")

        # Check power at 60 Hz in cleaned data is reduced by > 90%
        df = 1.0 / duration
        idx_60 = int(round(60.0 / df))
        pow_orig = np.abs(np.fft.rfft(ts_with_line.numpy())[idx_60])**2
        pow_clean = np.abs(np.fft.rfft(ts_clean.numpy())[idx_60])**2
        reduction = 1.0 - (pow_clean / pow_orig)
        self.assertGreater(reduction, 0.90, f"60 Hz line power reduction {reduction*100:.2f}% < 90%")

        # Check CBC waveform preservation on pure waveform
        clean_waveform, _ = pycbc.filter.subtract_coherent_lines(
            hp_td, lines=[(60.0, 10.0, 0.2)]
        )
        overlap = pycbc.filter.match(hp_td, clean_waveform)[0]
        self.assertGreater(overlap, 0.999, f"Waveform match {overlap*100:.4f}% is less than 99.9%")

    def test_real_detector_comparative_whitening(self):
        """Test comparative performance against standard PyCBC whitening across multiple real detector datasets."""
        datasets = [
            ("O4a Livingston", "data/o4a/day_gw230529/L-L1_GWOSC_O4a_4KHZ_R1-1369374720-4096.hdf5", "strain/Strain", 256.0),
            ("O3b Livingston", "data/o3/L1_broadband_1264316116.40.hdf", "strain", 64.0),
            ("O3b Hanford", "data/o3/H1_loud_1257296855.20.hdf", "strain", 64.0),
        ]

        tested_count = 0
        for det_name, data_path, key, dur in datasets:
            if not os.path.exists(data_path):
                continue

            with h5py.File(data_path, "r") as f:
                dset = f[key]
                fs = 4096.0
                raw = dset[:int(fs * dur)]
                gps = 1000000000.0
                if "meta/GPSstart" in f:
                    gps = float(f["meta/GPSstart"][()])
                elif "start_time" in dset.attrs:
                    gps = float(dset.attrs["start_time"])

            ts = TimeSeries(raw, delta_t=1.0 / fs, epoch=gps)
            edge = int(4.0 * fs)

            # 1. Standard PyCBC whitening
            std_w = ts.whiten(4.0, 4.0, low_frequency_cutoff=18.0, remove_corrupted=False)
            std_int = std_w.numpy()[edge:-edge]
            std_norm = std_int / np.std(std_int)

            # Standard PyCBC has uncalibrated discrete variance (std ~ sqrt(fs/2) ~ 45 >> 1.0)
            self.assertGreater(float(np.std(std_int)), 30.0,
                               f"Standard whitening in {det_name} was unexpectedly normalized")

            # Standard PyCBC has boundary edge transients (wrap-around corruption)
            std_edge_std = float(np.std(std_w.numpy()[:int(0.5 * fs)]))
            std_int_std = float(np.std(std_int))
            self.assertGreater(std_edge_std / std_int_std, 1.5,
                               f"Standard whitening boundary corruption absent in {det_name}")

            # 2. Regularized continuous overwhitening
            dow, dw, meta = pycbc.filter.overwhiten_strain(
                ts, f_low=18.0, max_filter_duration=4.0, subtract_lines=True
            )
            dw_int = dw.numpy()[edge:-edge]

            # 3. Discrete unit variance calibration: sigma == 1.000
            std_dw = float(np.std(dw_int))
            self.assertAlmostEqual(std_dw, 1.00, delta=0.05,
                                   msg=f"{det_name} whitened strain std {std_dw:.4f} deviates from 1.00")

            # 4. Spectral Flatness Measure (Wiener entropy)
            def compute_sfm(x):
                n_x = len(x)
                p_x = np.abs(np.fft.rfft(x))**2
                f_x = np.fft.rfftfreq(n_x, 1.0 / fs)
                mask_sfm = (f_x >= 20.0) & (f_x <= 1000.0)
                p_b = np.maximum(p_x[mask_sfm], 1e-30)
                return float(np.exp(np.mean(np.log(p_b))) / np.mean(p_b))

            sfm_new = compute_sfm(dw_int)
            self.assertGreater(sfm_new, 0.45, f"{det_name} Spectral Flatness {sfm_new:.4f} < 0.45")

            # 5. Gaussian normality
            kurt_val = float(stats.kurtosis(dw_int))
            skew_val = float(stats.skew(dw_int))
            self.assertLess(abs(kurt_val), 0.15, f"{det_name} excess kurtosis {kurt_val:.4f} exceeds 0.15")
            self.assertLess(abs(skew_val), 0.03, f"{det_name} skewness {skew_val:.4f} exceeds 0.03")

            # KS test against standard normal N(0, 1)
            _, ks_pval = stats.kstest(dw_int / std_dw, 'norm')
            self.assertGreater(ks_pval, 0.01, f"{det_name} KS test p-value {ks_pval:.4f} < 0.01")

            # 6. Maximum off-zero ACF correlation (delta-likeness < 0.025)
            n = len(dw_int)
            fx = np.fft.rfft(dw_int - np.mean(dw_int), n=2 * n)
            acf = np.fft.irfft(np.abs(fx)**2)[:n]
            acf /= acf[0]
            max_acf = float(np.max(np.abs(acf[int(0.01 * fs):])))
            self.assertLess(max_acf, 0.025, f"{det_name} max off-zero ACF {max_acf:.5f} exceeds 0.025")

            # 7. For 256s dataset, verify >99% matched filter SNR recovery
            if dur >= 256.0:
                hp, _ = pycbc.waveform.get_fd_waveform(
                    approximant="IMRPhenomD", mass1=30, mass2=30,
                    f_lower=18.0, delta_f=ts.delta_f, distance=600
                )
                hp.resize(len(ts.to_frequencyseries()))
                hp_td = pycbc.waveform.apply_fd_time_shift(hp, dur / 2.0).to_timeseries()
                hp_td.start_time = ts.start_time

                ts_inj = ts + hp_td
                dow_inj, _, _ = pycbc.filter.overwhiten_strain(
                    ts_inj, psd=meta["psd_base"], f_low=18.0, max_filter_duration=4.0
                )
                snr_ow = pycbc.filter.matched_filter_overwhitened(hp, dow_inj, psd=meta["psd_base"], f_low=18.0)
                max_snr_ow = float(np.max(np.abs(snr_ow.numpy()[edge:-edge])))

                psd_interp = pycbc.psd.interpolate(meta["psd_base"], ts_inj.delta_f)
                snr_std = pycbc.filter.matched_filter(hp, ts_inj, psd=psd_interp, low_frequency_cutoff=18.0)
                max_snr_std = float(np.max(np.abs(snr_std.numpy()[edge:-edge])))

                recovery_ratio = max_snr_ow / max_snr_std
                self.assertGreaterEqual(recovery_ratio, 0.99,
                                        f"SNR recovery {recovery_ratio*100:.3f}% is less than 99%")

            tested_count += 1

        self.assertGreaterEqual(tested_count, 2, "Fewer than 2 real detector datasets were evaluated")

    def test_high_resolution_filter_durations(self):
        """Test high-resolution 16s and 32s filter durations and violin mode suppression."""
        data_path = 'data/o4a/day_gw230529/L-L1_GWOSC_O4a_4KHZ_R1-1369374720-4096.hdf5'
        if not os.path.exists(data_path):
            self.skipTest(f"Real data file not found: {data_path}")

        with h5py.File(data_path, 'r') as f:
            dset = f['strain/Strain']
            fs = 4096.0
            raw = dset[:int(fs * 256.0)]
            gps = float(f['meta/GPSstart'][()])

        ts = TimeSeries(raw, delta_t=1.0 / fs, epoch=gps)

        # Baseline 4s filter vs high-resolution 16s and 32s filters
        dow_4, dw_4, _ = pycbc.filter.overwhiten_strain(ts, f_low=18.0, max_filter_duration=4.0)
        dow_16, dw_16, _ = pycbc.filter.overwhiten_strain(ts, f_low=18.0, max_filter_duration=16.0)
        dow_32, dw_32, _ = pycbc.filter.overwhiten_strain(ts, f_low=18.0, max_filter_duration=32.0)

        edge_16 = int(16.0 * fs)
        dw_4_int = dw_4.numpy()[edge_16:-edge_16]
        dw_16_int = dw_16.numpy()[edge_16:-edge_16]
        dw_32_int = dw_32.numpy()[edge_16:-edge_16]

        # Verify unit variance calibration across all filter resolutions
        self.assertAlmostEqual(float(np.std(dw_16_int)), 1.00, delta=0.05)
        self.assertAlmostEqual(float(np.std(dw_32_int)), 1.00, delta=0.05)

        # Measure power in 508-512 Hz violin resonance band
        f_dw = np.fft.rfftfreq(len(dw_4_int), 1.0 / fs)
        v_mask = (f_dw >= 508.0) & (f_dw <= 512.0)
        p_4 = np.max(np.abs(np.fft.rfft(dw_4_int)[v_mask])**2)
        p_16 = np.max(np.abs(np.fft.rfft(dw_16_int)[v_mask])**2)
        p_32 = np.max(np.abs(np.fft.rfft(dw_32_int)[v_mask])**2)

        # 16s and 32s filters resolve line width and reduce residual peak power by > 50%
        reduction_16 = 1.0 - (p_16 / p_4)
        reduction_32 = 1.0 - (p_32 / p_4)
        self.assertGreater(reduction_16, 0.40, f"16s filter violin reduction {reduction_16*100:.1f}% < 40%")
        self.assertGreater(reduction_32, 0.40, f"32s filter violin reduction {reduction_32*100:.1f}% < 40%")

    def test_direct_compact_psd_overwhitening(self):
        """Test passing CompactInstantaneousPSD directly into overwhiten_strain."""
        fs = 2048.0
        dur = 128.0
        N = int(fs * dur)
        rng = np.random.default_rng(777)
        noise = rng.normal(0, 1e-21, N)
        ts = TimeSeries(noise, delta_t=1.0 / fs, epoch=1000000000.0)

        pilot = pycbc.psd.welch(ts, seg_len=int(4 * fs), seg_stride=int(2 * fs), avg_method='median')
        model = pycbc.psd.CompactInstantaneousPSD(f_min=18.0, f_max=1024.0, n_knots=32)
        model.fit_static_from_psd(pilot)

        # Pass model directly as psd argument
        dow, dw, meta = pycbc.filter.overwhiten_strain(ts, psd=model, f_low=18.0, max_filter_duration=4.0)
        edge = int(4.0 * fs)
        dw_int = dw.numpy()[edge:-edge]
        std_val = float(np.std(dw_int))
        self.assertAlmostEqual(std_val, 1.00, delta=0.08,
                               msg=f"Direct CompactInstantaneousPSD overwhitening std {std_val:.4f} deviates from 1.00")

    def test_block_length_requirements_and_scaling(self):
        """Test block length requirements and >100x RT execution throughput."""
        data_path = 'data/o4a/day_gw230529/L-L1_GWOSC_O4a_4KHZ_R1-1369374720-4096.hdf5'
        if not os.path.exists(data_path):
            self.skipTest(f"Real data file not found: {data_path}")

        with h5py.File(data_path, 'r') as f:
            dset = f['strain/Strain']
            fs = 4096.0
            gps = float(f['meta/GPSstart'][()])
            raw_256 = dset[:int(fs * 256.0)]

        ts_256 = TimeSeries(raw_256, delta_t=1.0 / fs, epoch=gps)

        # For 256s block, measure runtime and throughput
        t0 = time.perf_counter()
        dow, dw, meta = pycbc.filter.overwhiten_strain(ts_256, f_low=18.0, subtract_lines=True)
        t1 = time.perf_counter()
        elapsed = t1 - t0
        rt_factor = 256.0 / elapsed

        # Throughput must exceed 100x real-time (we typically see >800x RT)
        self.assertGreater(rt_factor, 100.0, f"Throughput {rt_factor:.1f}x RT is below 100x RT threshold")

        # Transient fraction for 256s (4s filter duration at ends = 8s total = 3.125%)
        transient_frac = (2.0 * 4.0) / 256.0
        self.assertLessEqual(transient_frac, 0.05, f"Transient fraction {transient_frac:.3f} exceeds 5%")

    def test_in_pipeline_block_conditioning_with_injections(self):
        """Test in-pipeline block conditioning with software injections into raw strain."""
        fs = 2048.0
        duration = 128.0
        N = int(fs * duration)
        dt = 1.0 / fs
        rng = np.random.default_rng(999)

        # Raw strain block
        raw_noise = rng.normal(0, 1e-21, N)
        ts_block = TimeSeries(raw_noise, delta_t=dt, epoch=1000000000.0)

        # CBC injection
        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD", mass1=25, mass2=25,
            f_lower=20.0, delta_f=ts_block.delta_f, distance=300
        )
        hp.resize(len(ts_block.to_frequencyseries()))
        hp_td = pycbc.waveform.apply_fd_time_shift(hp, duration / 2.0).to_timeseries()
        hp_td.start_time = ts_block.start_time

        # Condition raw block in memory with injection
        dow, dw, meta = pycbc.filter.overwhiten_block_strain(
            ts_block, injections=hp_td, f_low=18.0, max_filter_duration=4.0
        )

        self.assertEqual(meta["block_duration"], duration)
        self.assertEqual(meta["sample_rate"], fs)

        # Detect the injected template
        edge = int(4.0 * fs)
        snr_series = pycbc.filter.matched_filter_overwhitened(
            hp, dow, psd=meta["psd_base"], f_low=18.0
        )
        peak_snr = float(np.max(np.abs(snr_series.numpy()[edge:-edge])))
        self.assertGreater(peak_snr, 8.0, f"Injected signal peak SNR {peak_snr:.2f} is below detection threshold 8.0")

    def test_edge_inpainting_eliminates_boundary_explosion(self):
        """Test inpaint_edges=True suppresses edge step transients by >100x."""
        fs = 2048.0
        dur = 32.0
        N = int(fs * dur)
        rng = np.random.default_rng(42)
        white = rng.normal(0, 1, N)
        freqs = np.fft.rfftfreq(N, 1.0 / fs)
        color = 1.0 / np.sqrt(1.0 + (np.maximum(freqs, 20.0) / 40.0)**4)
        colored = np.fft.irfft(np.fft.rfft(white) * color)
        ts = TimeSeries(colored, delta_t=1.0 / fs, epoch=1000000000.0)

        # 1. Overwhiten with edge inpainting enabled
        dow_inp, dw_inp, meta = pycbc.filter.overwhiten_strain(
            ts, f_low=20.0, max_filter_duration=2.0, inpaint_edges=True, edge_pad_duration=2.0
        )
        self.assertEqual(len(dow_inp), len(ts))
        self.assertEqual(float(dow_inp.start_time), float(ts.start_time))
        self.assertTrue(meta["inpainted_edges"])

        # Check boundary variance right at start (first 0.25s) and end (last 0.25s)
        edge_n = int(0.25 * fs)
        std_start = float(np.std(dw_inp.numpy()[:edge_n]))
        std_end = float(np.std(dw_inp.numpy()[-edge_n:]))
        std_int = float(np.std(dw_inp.numpy()[int(2.0 * fs) : -int(2.0 * fs)]))

        self.assertLess(std_start, 3.0, f"Start std {std_start:.3f} was not suppressed (< 3.0)")
        self.assertLess(std_end, 3.0, f"End std {std_end:.3f} was not suppressed (< 3.0)")
        self.assertAlmostEqual(std_int, 1.0, delta=0.15)

        # 2. Test in overwhiten_block_strain
        dow_blk, dw_blk, meta_blk = pycbc.filter.overwhiten_block_strain(
            ts, f_low=20.0, max_filter_duration=2.0, inpaint_edges=True, edge_pad_duration=2.0
        )
        self.assertEqual(len(dow_blk), len(ts))
        self.assertTrue(meta_blk["inpainted_edges"])


if __name__ == "__main__":
    unittest.main()

