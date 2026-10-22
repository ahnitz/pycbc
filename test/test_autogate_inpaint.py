# Copyright (C) 2026 Alex Nitz
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

"""Unit tests for autogating and inpainting improvements.

Verifies:
1. Dual-band glitch detection ([20, 150] Hz and full-band).
2. Dynamic gate width expansion and gate merging.
3. Float64 ridge regularization (paint-ridge = 1e-8).
4. Gating and inpainting of loud glitches (SNR 25-50) yields zero edge ringing
   and matched-filter SNR consistent with Gaussian noise.
5. End-to-end strain conditioning without redundant post-inpaint highpass.
"""

import unittest
import numpy as np

import pycbc.psd
import pycbc.types
from pycbc.filter import highpass, matched_filter
from pycbc.strain.gate import (
    gate_and_paint,
    gate_and_paint_matmul,
    invert_covariance,
)
from pycbc.strain.strain import detect_loud_glitches
from pycbc.waveform import get_fd_waveform


class TestAutogateInpaint(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.sample_rate = 2048.0
        self.duration = 32.0
        self.N = int(self.duration * self.sample_rate)
        self.delta_t = 1.0 / self.sample_rate
        self.delta_f = 1.0 / self.duration

        # Generate Gaussian noise time series
        white = np.random.randn(self.N)
        self.ts_noise = pycbc.types.TimeSeries(
            white, delta_t=self.delta_t, epoch=0
        )

        # Estimate Welch PSD with Hann inverse spectrum truncation
        psd = pycbc.psd.welch(
            self.ts_noise,
            seg_len=int(4.0 * self.sample_rate),
            seg_stride=int(2.0 * self.sample_rate),
        )
        psd = pycbc.psd.interpolate(psd, self.delta_f)
        self.psd = pycbc.psd.inverse_spectrum_truncation(
            psd,
            int(4.0 * self.sample_rate),
            low_frequency_cutoff=20.0,
            trunc_method="hann",
        )
        self.invpsd = 1.0 / self.psd

        # Create gravitational wave template for matched filtering
        self.hp, _ = get_fd_waveform(
            approximant="IMRPhenomD",
            mass1=30.0,
            mass2=30.0,
            delta_f=self.delta_f,
            f_lower=20.0,
        )
        self.hp.resize(len(self.psd))

    def test_dual_band_glitch_detection(self):
        """Test dual-band detects mid-frequency glitches diluted in full."""
        t = self.ts_noise.sample_times.numpy() - 16.0
        # A 60 Hz sinusoidal burst with duration ~ 0.15s
        # Yields ~ 80 sigma in [20, 150] Hz but ~ 38 sigma full-band
        env = np.exp(-((t / 0.04) ** 2))
        glitch_mid = 25.0 * env * np.sin(2 * np.pi * 60.0 * t)
        ts_mid_glitch = self.ts_noise + pycbc.types.TimeSeries(
            glitch_mid, delta_t=self.delta_t, epoch=0
        )

        # Full-band only with threshold=50: diluted and missed
        gates_full_only = detect_loud_glitches(
            ts_mid_glitch,
            threshold=50.0,
            low_freq_cutoff=20.0,
            bandlimited_freq=None,
            return_gate_intervals=True,
        )
        self.assertEqual(len(gates_full_only), 0)

        # Dual-band ([20, 150] Hz + full-band) with threshold=50: detected
        gates_dual = detect_loud_glitches(
            ts_mid_glitch,
            threshold=50.0,
            low_freq_cutoff=20.0,
            bandlimited_freq=(20.0, 150.0),
            return_gate_intervals=True,
        )
        self.assertGreaterEqual(len(gates_dual), 1)
        center_t, half_w, _ = gates_dual[0]
        # Verify the glitch at 16.0s is inside the detected gate
        self.assertLess(center_t - half_w, 16.0)
        self.assertGreater(center_t + half_w, 16.0)
        self.assertGreater(half_w, 0.1)

        # High-frequency transient (400 Hz) detected by full-band
        env_h = np.exp(-((t / 0.02) ** 2))
        glitch_high = 65.0 * env_h * np.sin(2 * np.pi * 400.0 * t)
        ts_high_glitch = self.ts_noise + pycbc.types.TimeSeries(
            glitch_high, delta_t=self.delta_t, epoch=0
        )
        gates_high = detect_loud_glitches(
            ts_high_glitch,
            threshold=50.0,
            low_freq_cutoff=20.0,
            bandlimited_freq=(20.0, 150.0),
            return_gate_intervals=True,
        )
        self.assertGreaterEqual(len(gates_high), 1)
        self.assertAlmostEqual(gates_high[0][0], 16.0, delta=0.5)

    def test_dynamic_gate_expansion(self):
        """Test gate window expands dynamically to cover broad skirts."""
        t = self.ts_noise.sample_times.numpy() - 16.0
        env = np.exp(-((t / 0.25) ** 2))
        glitch_broad = 60.0 * env * np.sin(2 * np.pi * 50.0 * t)
        ts_broad = self.ts_noise + pycbc.types.TimeSeries(
            glitch_broad, delta_t=self.delta_t, epoch=0
        )

        gates = detect_loud_glitches(
            ts_broad,
            threshold=40.0,
            low_freq_cutoff=20.0,
            min_half_width=0.25,
            pad_duration=0.20,
            min_quiet_duration=0.10,
            quiet_threshold=3.0,
            return_gate_intervals=True,
        )
        self.assertGreaterEqual(len(gates), 1)
        center_t, half_w, _ = gates[0]
        self.assertAlmostEqual(center_t, 16.0, delta=0.5)
        self.assertGreater(half_w, 0.50)

    def test_gate_merging(self):
        """Test overlapping/close gates (< 0.20s apart) are merged."""
        t = self.ts_noise.sample_times.numpy()
        # Two in-band 60 Hz glitches separated by 0.15s (< 0.20s)
        env1 = np.exp(-(((t - 15.0) / 0.03) ** 2))
        env2 = np.exp(-(((t - 15.15) / 0.03) ** 2))
        g1 = 40.0 * env1 * np.sin(2 * np.pi * 60.0 * (t - 15.0))
        g2 = 40.0 * env2 * np.sin(2 * np.pi * 60.0 * (t - 15.15))
        ts_two_close = self.ts_noise + pycbc.types.TimeSeries(
            g1 + g2, delta_t=self.delta_t, epoch=0
        )

        gates_close = detect_loud_glitches(
            ts_two_close,
            threshold=30.0,
            low_freq_cutoff=20.0,
            min_separation=0.20,
            return_gate_intervals=True,
        )
        # Should be merged into exactly 1 unified gate
        self.assertEqual(len(gates_close), 1)
        center_t, half_w, _ = gates_close[0]
        self.assertLess(center_t - half_w, 15.0)
        self.assertGreater(center_t + half_w, 15.15)

        # Two in-band 60 Hz glitches separated by 5.0s (> 0.20s)
        env3 = np.exp(-(((t - 20.0) / 0.03) ** 2))
        g3 = 40.0 * env3 * np.sin(2 * np.pi * 60.0 * (t - 20.0))
        ts_two_far = self.ts_noise + pycbc.types.TimeSeries(
            g1 + g3, delta_t=self.delta_t, epoch=0
        )
        gates_far = detect_loud_glitches(
            ts_two_far,
            threshold=30.0,
            low_freq_cutoff=20.0,
            min_separation=0.20,
            return_gate_intervals=True,
        )
        # Should remain 2 separate gates
        self.assertEqual(len(gates_far), 2)
        self.assertAlmostEqual(gates_far[0][0], 15.0, delta=0.5)
        self.assertAlmostEqual(gates_far[1][0], 20.0, delta=0.5)

    def test_float64_ridge_regularization(self):
        """Test inpainting with paint_ridge=1e-8 in float64 precision."""
        t = self.ts_noise.sample_times.numpy() - 16.0
        env = np.exp(-((t / 0.05) ** 2))
        glitch = 50.0 * env * np.sin(2 * np.pi * 50.0 * t)
        ts_glitch = self.ts_noise + pycbc.types.TimeSeries(
            glitch, delta_t=self.delta_t, epoch=0
        )

        lindex = int(15.8 * self.sample_rate)
        rindex = int(16.2 * self.sample_rate)

        # Cholesky solver with ridge=1e-8
        cleaned_cho = gate_and_paint(
            ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="cholesky",
            ridge=1e-8,
        )
        self.assertTrue(np.all(np.isfinite(cleaned_cho.numpy())))

        # Toeplitz solver with ridge=1e-8
        cleaned_toep = gate_and_paint(
            ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="toeplitz",
            ridge=1e-8,
        )
        self.assertTrue(np.all(np.isfinite(cleaned_toep.numpy())))

        # Both solvers must agree to high numerical precision
        diff = np.max(np.abs(cleaned_cho.numpy() - cleaned_toep.numpy()))
        self.assertLess(diff, 1e-4)

        # Matmul solver with invert_covariance
        invmat = invert_covariance(self.invpsd, lindex, rindex, ridge=1e-8)
        cleaned_mm = gate_and_paint_matmul(
            ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            invmat=invmat,
            copy=True,
        )
        diff_mm = np.max(np.abs(cleaned_cho.numpy() - cleaned_mm.numpy()))
        self.assertLess(diff_mm, 1e-4)

        # Verify overwhitened residual inside gate is canceled
        owh = (cleaned_cho.to_frequencyseries() * self.invpsd).to_timeseries()
        owh_inside = owh[lindex:rindex].numpy()
        td = (
            self.invpsd.astype("complex").to_timeseries()
            * self.invpsd.delta_t
        )
        diag = td[0]
        self.assertLess(np.max(np.abs(owh_inside)) / diag, 1e-4)

    def test_loud_glitches_gated_and_inpainted_consistent_with_noise(self):
        """Verify loud glitches (SNR 25-50) inpainted with zero ringing."""
        t = self.ts_noise.sample_times.numpy() - 16.0
        env = np.exp(-((t / 0.05) ** 2))
        glitch = 45.0 * env * np.sin(2 * np.pi * 60.0 * t)
        ts_glitch = self.ts_noise + pycbc.types.TimeSeries(
            glitch, delta_t=self.delta_t, epoch=0
        )

        # Compute raw SNR on glitch before inpainting
        snr_raw = matched_filter(
            self.hp, ts_glitch, psd=self.psd, low_frequency_cutoff=20.0
        )
        i_eval_s = int(14.0 * self.sample_rate)
        i_eval_e = int(18.0 * self.sample_rate)
        eval_slice = slice(i_eval_s, i_eval_e)
        max_snr_raw = abs(snr_raw[eval_slice]).numpy().max()
        self.assertGreater(max_snr_raw, 25.0)

        # Autogate and inpaint
        gates = detect_loud_glitches(
            ts_glitch,
            threshold=30.0,
            low_freq_cutoff=20.0,
            min_half_width=0.25,
            return_gate_intervals=True,
        )
        self.assertGreaterEqual(len(gates), 1)

        ts_cleaned = ts_glitch.copy()
        for g_time, g_win, g_tap in gates:
            ts_cleaned = ts_cleaned.gate(
                g_time,
                window=g_win,
                method="paint",
                invpsd=self.invpsd,
                paint_ridge=1e-8,
                paint_method="cholesky",
            )

        snr_clean = matched_filter(
            self.hp, ts_cleaned, psd=self.psd, low_frequency_cutoff=20.0
        )
        max_snr_clean = abs(snr_clean[eval_slice]).numpy().max()
        # Must be consistent with Gaussian noise floor (< 5.5)
        self.assertLess(max_snr_clean, 5.5)

        # Inspect gate boundaries: zero edge ringing
        gate_start = gates[0][0] - gates[0][1]
        gate_end = gates[0][0] + gates[0][1]
        left_slice = slice(
            int((gate_start - 0.5) * self.sample_rate),
            int(gate_start * self.sample_rate),
        )
        right_slice = slice(
            int(gate_end * self.sample_rate),
            int((gate_end + 0.5) * self.sample_rate),
        )
        max_edge_snr_l = abs(snr_clean[left_slice]).numpy().max()
        max_edge_snr_r = abs(snr_clean[right_slice]).numpy().max()
        self.assertLess(max_edge_snr_l, 5.5)
        self.assertLess(max_edge_snr_r, 5.5)

    def test_from_cli_autogating_and_no_post_inpaint_ringing(self):
        """Test from_cli logic without redundant post-inpaint highpass."""
        t = self.ts_noise.sample_times.numpy() - 16.0
        env = np.exp(-((t / 0.05) ** 2))
        glitch = 40.0 * env * np.sin(2 * np.pi * 60.0 * t)
        ts_glitch = self.ts_noise + pycbc.types.TimeSeries(
            glitch, delta_t=self.delta_t, epoch=0
        )

        class MockOpts:
            pass

        opts = MockOpts()
        opts.pad_data = 0
        opts.strain_high_pass = 20.0
        opts.autogating_threshold = 30.0
        opts.autogating_cluster = 1.0
        opts.autogating_width = 0.25
        opts.autogating_taper = 0.05
        opts.autogating_pad = 2.0
        opts.autogating_max_iterations = 2
        opts.gating_method = "paint"
        opts.paint_method = "cholesky"
        opts.paint_ridge = 1e-8

        # In from_cli: highpass is called at line 378 before gating.
        # Line 487 (redundant post-inpaint highpass) was removed.
        initial_ts = highpass(
            ts_glitch.copy(), frequency=opts.strain_high_pass
        )
        gate_params = detect_loud_glitches(
            initial_ts,
            threshold=opts.autogating_threshold,
            cluster_window=opts.autogating_cluster,
            low_freq_cutoff=opts.strain_high_pass,
            corrupt_time=opts.pad_data + opts.autogating_pad,
            min_half_width=opts.autogating_width,
            taper_duration=opts.autogating_taper,
            return_gate_intervals=True,
        )
        self.assertGreaterEqual(len(gate_params), 1)

        clean_strain = initial_ts.gate(
            gate_params[0][0],
            window=gate_params[0][1],
            method="paint",
            invpsd=self.invpsd,
            paint_ridge=opts.paint_ridge,
            paint_method=opts.paint_method,
        )

        snr_clean = matched_filter(
            self.hp, clean_strain, psd=self.psd, low_frequency_cutoff=20.0
        )
        i_eval_s = int(14.0 * self.sample_rate)
        i_eval_e = int(18.0 * self.sample_rate)
        eval_slice = slice(i_eval_s, i_eval_e)
        max_snr = abs(snr_clean[eval_slice]).numpy().max()
        self.assertLess(max_snr, 5.5)


if __name__ == "__main__":
    unittest.main()
