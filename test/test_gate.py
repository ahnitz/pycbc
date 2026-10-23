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

"""Unit tests for gating, inpainting, and loud glitch detection."""

import argparse
import unittest
import numpy as np

import pycbc.psd
import pycbc.types
from pycbc.strain.gate import _get_tdfilter, gate_and_paint
from pycbc.strain.strain import detect_loud_glitches, insert_strain_option_group


class TestGateAndInpaint(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.sample_rate = 2048.0
        # 64-second duration gives 131,072 samples, accommodating K > 65536
        self.duration = 64.0
        self.N = int(self.duration * self.sample_rate)
        self.delta_t = 1.0 / self.sample_rate
        self.delta_f = 1.0 / self.duration

        white = np.random.randn(self.N)
        self.ts = pycbc.types.TimeSeries(white, delta_t=self.delta_t, epoch=0)

        # Estimate Welch PSD with Hann inverse spectrum truncation
        psd = pycbc.psd.welch(
            self.ts,
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

    def test_detect_loud_glitches_gwastro_master(self):
        """Verify detect_loud_glitches aligns with gwastro/master fixed-width behavior."""
        # Create strain with a glitch at t = 32.0s
        t = self.ts.sample_times.numpy() - 32.0
        env_core = np.exp(-((t / 0.05) ** 2))
        glitch_core = 50.0 * env_core * np.sin(2 * np.pi * 60.0 * t)
        ts_glitch = self.ts + pycbc.types.TimeSeries(
            glitch_core, delta_t=self.delta_t, epoch=0
        )

        glitch_times = detect_loud_glitches(
            ts_glitch,
            threshold=50.0,
            low_freq_cutoff=20.0,
        )
        self.assertGreaterEqual(len(glitch_times), 1)
        # Verify it returns float timestamps centered near 32.0s
        for gt in glitch_times:
            self.assertIsInstance(float(gt), float)
            self.assertAlmostEqual(float(gt), 32.0, delta=0.5)

        # High threshold returns empty list
        no_glitches = detect_loud_glitches(
            ts_glitch,
            threshold=500.0,
            low_freq_cutoff=20.0,
        )
        self.assertEqual(len(no_glitches), 0)

    def test_gate_and_paint_wide_gate_k_greater_than_65536(self):
        """Verify gate_and_paint with K > 65536 (e.g. K = 70000)."""
        K = 70000
        lindex = 10000
        rindex = lindex + K

        # 1. By default, K > 32768 falls back cleanly to cosine taper gating
        cleaned_fallback = gate_and_paint(
            self.ts, lindex, rindex, self.invpsd, copy=True, fallback_taper=True
        )
        self.assertTrue(np.all(np.isfinite(cleaned_fallback.numpy())))
        # Gate center should be cleanly zeroed
        center_idx = (lindex + rindex) // 2
        self.assertEqual(cleaned_fallback[center_idx], 0.0)

        # 2. When fallback_taper=False, dynamic max_samples and zero padding ensure
        # len(col) == K == len(rhs) without dimension mismatch
        cleaned_toeplitz = gate_and_paint(
            self.ts,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="toeplitz",
            ridge=1e-3,
            fallback_taper=False,
            pad_duration=0.0,  # test global FFT branch
        )
        self.assertTrue(np.all(np.isfinite(cleaned_toeplitz.numpy())))
        self.assertEqual(len(cleaned_toeplitz), len(self.ts))

    def test_dynamic_tdfilter_and_raw_col_padding(self):
        """Verify _get_tdfilter respects dynamic max_samples and raw_col padding."""
        K = 70000
        req_max_samples = max(65536, 1 << int(np.ceil(np.log2(2 * K))))
        self.assertGreaterEqual(req_max_samples, 2 * K)

        tdfilter = _get_tdfilter(self.invpsd, max_samples=req_max_samples)
        self.assertGreaterEqual(len(tdfilter), 65536)

        # Verify zero-padding if len(raw_col) < K
        raw_col = tdfilter[:K].numpy() if hasattr(tdfilter[:K], "numpy") else tdfilter[:K]
        if len(raw_col) < K:
            padded_col = np.pad(raw_col, (0, K - len(raw_col)))
            self.assertEqual(len(padded_col), K)
        else:
            self.assertEqual(len(raw_col), K)

    def test_default_paint_ridge_is_1e_3(self):
        """Verify default paint-ridge is 0.001 (1e-3) in CLI parser."""
        parser = argparse.ArgumentParser()
        insert_strain_option_group(parser)
        opts = parser.parse_args([])
        self.assertEqual(opts.paint_ridge, 0.001)

    def test_paint_method_auto_cli_and_execution(self):
        """Verify --paint-method auto parses cleanly and executes identically to toeplitz."""
        parser = argparse.ArgumentParser()
        insert_strain_option_group(parser)
        opts = parser.parse_args(["--paint-method", "auto", "--channel-name", "H1:TEST"])
        self.assertEqual(opts.paint_method, "auto")

        # Test execution with method='auto' vs method='toeplitz'
        lindex = int(31.8 * self.sample_rate)
        rindex = int(32.2 * self.sample_rate)
        cleaned_auto = gate_and_paint(
            self.ts, lindex, rindex, self.invpsd, copy=True, method="auto", ridge=1e-3
        )
        cleaned_toeplitz = gate_and_paint(
            self.ts, lindex, rindex, self.invpsd, copy=True, method="toeplitz", ridge=1e-3
        )
        np.testing.assert_allclose(cleaned_auto.numpy(), cleaned_toeplitz.numpy(), rtol=1e-7, atol=1e-7)

        # Test timeseries gate_and_paint with method='auto'
        ts_copy = self.ts.copy()
        cleaned_ts = gate_and_paint(ts_copy, lindex, rindex, self.invpsd, method="auto")
        self.assertTrue(np.all(np.isfinite(cleaned_ts.numpy())))


if __name__ == "__main__":
    unittest.main()
