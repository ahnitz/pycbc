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

"""Unit tests for regularized gating, inpainting, and strain edge conditioning."""

import argparse
import unittest
from unittest import mock

import numpy as np
from scipy import linalg

import pycbc
import pycbc.strain.strain as strain_mod
from pycbc.filter import highpass, make_frequency_series
from pycbc.psd import interpolate, inverse_spectrum_truncation, welch
from pycbc.strain import (
    StrainSegments,
    from_cli,
    insert_strain_option_group,
    insert_strain_option_group_multi_ifo,
)
from pycbc.strain.gate import (
    gate_and_paint,
    gate_and_paint_matmul,
    gate_and_paint_multi,
    invert_covariance,
)
from pycbc.types import FrequencySeries, TimeSeries


class TestGateAndPaint(unittest.TestCase):
    """Tests for regularized gating/inpainting solvers and edge conditioning."""

    def setUp(self):
        rng = np.random.default_rng(12345)
        self.sample_rate = 2048.0
        self.duration = 8.0
        self.n_samples = int(self.duration * self.sample_rate)
        self.delta_t = 1.0 / self.sample_rate
        self.delta_f = 1.0 / self.duration

        white = rng.standard_normal(self.n_samples)
        self.ts = TimeSeries(white, delta_t=self.delta_t, epoch=0)

        # Build a realistic Hann-truncated inverse PSD
        psd = welch(
            self.ts,
            seg_len=int(2 * self.sample_rate),
            seg_stride=int(self.sample_rate),
        )
        psd = interpolate(psd, self.delta_f)
        self.psd = inverse_spectrum_truncation(
            psd,
            max_filter_len=int(2 * self.sample_rate),
            low_frequency_cutoff=20.0,
            trunc_method="hann",
        )
        self.invpsd = 1.0 / self.psd

        # Superpose a loud transient glitch at t = 4.0s
        times = self.ts.sample_times.numpy() - 4.0
        glitch = (
            50.0 * np.exp(-((times / 0.05) ** 2)) * np.sin(2 * np.pi * 50.0 * times)
        )
        self.ts_glitch = self.ts + TimeSeries(glitch, delta_t=self.delta_t, epoch=0)

    def test_inpaint_solvers_agree_and_cancel_overwhitened_hole(self):
        """Cholesky, Toeplitz, auto, and matmul solvers agree and zero C^-1 x."""
        lindex = int(3.95 * self.sample_rate)
        rindex = int(4.05 * self.sample_rate)

        cleaned_toep = gate_and_paint(
            self.ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="toeplitz",
            ridge=1e-10,
        )
        cleaned_auto = gate_and_paint(
            self.ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="auto",
            ridge=1e-10,
        )
        cleaned_cho = gate_and_paint(
            self.ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            copy=True,
            method="cholesky",
            ridge=1e-10,
        )
        invmat = invert_covariance(self.invpsd, lindex, rindex, ridge=1e-10)
        cleaned_mm = gate_and_paint_matmul(
            self.ts_glitch,
            lindex,
            rindex,
            self.invpsd,
            invmat=invmat,
            copy=True,
        )

        self.assertLess(
            np.max(np.abs(cleaned_toep.numpy() - cleaned_auto.numpy())), 1e-10
        )
        self.assertLess(
            np.max(np.abs(cleaned_cho.numpy() - cleaned_toep.numpy())), 1e-4
        )
        self.assertLess(np.max(np.abs(cleaned_mm.numpy() - cleaned_cho.numpy())), 1e-4)

        owh_clean = (cleaned_cho.to_frequencyseries() * self.invpsd).to_timeseries()
        diag = (self.invpsd.astype("complex").to_timeseries() * self.invpsd.delta_t)[0]
        self.assertLess(np.max(np.abs(owh_clean[lindex:rindex].numpy())) / diag, 1e-4)

    def test_ridge_regularizes_bandlimited_psd(self):
        """Bandlimited PSDs fail without ridge=0 and succeed with ridge > 0."""
        invpsd_data = np.zeros(self.n_samples // 2 + 1)
        invpsd_data[int(30 / self.delta_f) : int(300 / self.delta_f)] = 1.0
        invpsd = FrequencySeries(invpsd_data, delta_f=self.delta_f)

        lindex = int(3.9 * self.sample_rate)
        rindex = int(4.1 * self.sample_rate)

        with self.assertRaises(linalg.LinAlgError):
            gate_and_paint(
                self.ts, lindex, rindex, invpsd, method="cholesky", ridge=0.0
            )

        cleaned = gate_and_paint(
            self.ts, lindex, rindex, invpsd, method="cholesky", ridge=1e-7
        )
        self.assertTrue(np.all(np.isfinite(cleaned.numpy())))

    def test_ridge_preserves_high_dynamic_range_inpainting(self):
        """Default ridge=1e-7 preserves low-frequency nulling on steep PSDs."""
        freqs = self.invpsd.sample_frequencies.numpy()
        # Simulate LIGO-like 6-decade inverse-PSD dynamic range
        steep_inv = np.where(
            freqs >= 15.0, (freqs / 150.0) ** 4 / (1.0 + (freqs / 150.0) ** 6), 0.0
        )
        steep_inv = FrequencySeries(
            np.maximum(steep_inv, 1e-6 * np.max(steep_inv)) * (freqs >= 15.0),
            delta_f=self.delta_f,
        )
        lindex = int(3.8 * self.sample_rate)
        rindex = int(4.2 * self.sample_rate)
        diag = (steep_inv.astype("complex").to_timeseries() * steep_inv.delta_t)[0]

        cleaned = gate_and_paint(self.ts_glitch, lindex, rindex, steep_inv)
        owh = (cleaned.to_frequencyseries() * steep_inv).to_timeseries()
        self.assertLess(np.max(np.abs(owh[lindex:rindex].numpy())) / diag, 1e-4)

    def test_joint_multi_gate_dense_burst_stability(self):
        """Dense bursts separated by 1.2s solve cleanly in float32 via PCG."""
        rng = np.random.default_rng(99)
        sr = 2048.0
        dur = 64.0
        n = int(sr * dur)
        df = 1.0 / dur
        freqs = np.fft.rfftfreq(n, 1.0 / sr)
        steep_inv = np.where(
            freqs >= 11.0, (freqs / 150.0) ** 4 / (1.0 + (freqs / 150.0) ** 6), 0.0
        )
        steep_inv = FrequencySeries(
            (np.maximum(steep_inv, 1e-6 * np.max(steep_inv)) * (freqs >= 11.0)).astype(
                np.float32
            ),
            delta_f=df,
        )
        # Generate Gaussian noise matching the steep colored PSD, scaled to ~unit in-band RMS
        inv_arr = steep_inv.numpy().astype(np.float64)
        amp = np.where(freqs >= 15.0, 1.0 / np.sqrt(np.maximum(inv_arr, 1e-30)), 0.0)
        amp /= np.sqrt(np.sum(amp[freqs >= 20.0] ** 2) * 2.0 / (n * n)) + 1e-30
        noise_fd = (
            rng.standard_normal(len(freqs)) + 1j * rng.standard_normal(len(freqs))
        ) * amp
        noise_fd[0] = 0.0
        noise_fd[-1] = np.real(noise_fd[-1])
        data = TimeSeries(
            np.fft.irfft(noise_fd, n=n).astype(np.float32),
            delta_t=1.0 / sr,
            epoch=0,
        )

        intervals = []
        for k in range(12):
            center = int((20.0 + 1.2 * k) * sr)
            li, ri = center - 512, center + 512
            data[li:ri] += 500.0
            intervals.append((li, ri))

        cleaned = gate_and_paint_multi(data, intervals, steep_inv, copy=True)
        self.assertTrue(np.all(np.isfinite(cleaned.numpy())))
        self.assertLess(np.max(np.abs(cleaned.numpy())), 1e3)
        inband = highpass(cleaned, 20.0)
        self.assertLess(np.max(np.abs(inband.numpy())), 10.0)

    def test_wide_gate_switches_and_ultrawide_falls_back_to_taper(self):
        """Wide gates (>2048) switch cholesky->toeplitz; >32768 falls back to taper."""
        lindex = int(3.25 * self.sample_rate)
        rindex = int(4.75 * self.sample_rate)
        self.assertGreater(rindex - lindex, 2048)

        cleaned_cho = gate_and_paint(
            self.ts_glitch, lindex, rindex, self.invpsd, method="cholesky", ridge=1e-6
        )
        cleaned_toep = gate_and_paint(
            self.ts_glitch, lindex, rindex, self.invpsd, method="toeplitz", ridge=1e-6
        )
        self.assertLess(
            np.max(np.abs(cleaned_cho.numpy() - cleaned_toep.numpy())), 1e-6
        )

        # Ultra-wide gate (>32768 samples) uses taper fallback cleanly
        long_sr = 2048.0
        long_n = int(24.0 * long_sr)
        long_ts = TimeSeries(np.ones(long_n), delta_t=1.0 / long_sr, epoch=0)
        long_psd = FrequencySeries(np.ones(long_n // 2 + 1), delta_f=1.0 / 24.0)
        tapered = gate_and_paint(
            long_ts, 512, 512 + 33000, long_psd, copy=True, fallback_taper=True
        )
        self.assertTrue(np.all(np.isfinite(tapered.numpy())))
        self.assertEqual(np.max(np.abs(tapered[1500:32000].numpy())), 0.0)

    def test_timeseries_gate_paint_auto(self):
        """TimeSeries.gate supports paint_method='auto' and subsample shift off."""
        gated = self.ts_glitch.gate(
            4.0001,
            window=0.1,
            method="paint",
            invpsd=self.invpsd,
            paint_method="auto",
            paint_subsample_shift=False,
            copy=True,
        )
        self.assertTrue(np.all(np.isfinite(gated.numpy())))
        # Glitch peak near 4.0s should be strongly suppressed in-band (above 20 Hz)
        self.assertLess(np.max(np.abs(highpass(gated, 25.0).numpy())), 10.0)

    def test_strain_segments_inpaint_edges(self):
        """StrainSegments inpaints zero-padded start/end boundaries when enabled."""
        rng = np.random.default_rng(7)
        sr = 2048.0
        dur = 16.0
        strain = TimeSeries(
            rng.standard_normal(int(sr * dur)), delta_t=1.0 / sr, epoch=100
        )

        # Trigger range starts/ends at strain bounds with 4s pads, forcing zero padding
        segs_plain = StrainSegments(
            strain,
            segment_length=8,
            segment_start_pad=2,
            segment_end_pad=2,
            trigger_start=100,
            trigger_end=116,
            allow_zero_padding=True,
            inpaint_edges=False,
        )
        segs_inpaint = StrainSegments(
            strain,
            segment_length=8,
            segment_start_pad=2,
            segment_end_pad=2,
            trigger_start=100,
            trigger_end=116,
            allow_zero_padding=True,
            inpaint_edges=True,
            low_frequency_cutoff=20.0,
        )

        fd_plain = segs_plain.fourier_segments()
        fd_inpaint = segs_inpaint.fourier_segments()
        self.assertEqual(len(fd_plain), len(fd_inpaint))

        first_plain_td = fd_plain[0].to_timeseries().numpy()
        first_inpaint_td = fd_inpaint[0].to_timeseries().numpy()
        # First segment starts with 2s (4096 samples) of zero-padding
        self.assertLess(np.max(np.abs(first_plain_td[:4096])), 1e-12)
        self.assertGreater(np.std(first_inpaint_td[:4096]), 0.05)

        last_plain_td = fd_plain[-1].to_timeseries().numpy()
        last_inpaint_td = fd_inpaint[-1].to_timeseries().numpy()
        # Last segment ends with 2s of zero-padding
        self.assertLess(np.max(np.abs(last_plain_td[-4096:])), 1e-12)
        self.assertGreater(np.std(last_inpaint_td[-4096:]), 0.05)

    def test_from_cli_paint_skips_second_highpass(self):
        """pycbc.strain.from_cli with gating_method='paint' omits post-gate highpass."""
        parser = argparse.ArgumentParser()
        insert_strain_option_group(parser)
        base_argv = [
            "--gps-start-time",
            "0",
            "--gps-end-time",
            "64",
            "--pad-data",
            "8",
            "--sample-rate",
            "2048",
            "--channel-name",
            "H1:TEST",
            "--fake-strain",
            "aLIGODesignSensitivityP1200087",
            "--fake-strain-seed",
            "42",
            "--strain-high-pass",
            "11.0",
            "--autogating-threshold",
            "2.0",
            "--autogating-pad",
            "1.0",
            "--autogating-width",
            "0.1",
        ]
        args_paint = parser.parse_args(
            [
                *base_argv,
                "--gating-method",
                "paint",
                "--paint-method",
                "auto",
                "--paint-ridge",
                "0.001",
            ]
        )
        args_taper = parser.parse_args([*base_argv, "--gating-method", "taper"])

        with mock.patch.object(
            strain_mod, "highpass", wraps=strain_mod.highpass
        ) as hp_paint:
            strain_paint = from_cli(args_paint, dyn_range_fac=pycbc.DYN_RANGE_FAC)
            self.assertEqual(hp_paint.call_count, 1)

        with mock.patch.object(
            strain_mod, "highpass", wraps=strain_mod.highpass
        ) as hp_taper:
            from_cli(args_taper, dyn_range_fac=pycbc.DYN_RANGE_FAC)
            self.assertEqual(hp_taper.call_count, 2)

        self.assertGreater(len(strain_paint.gating_info["auto"]), 0)
        gate_t, gate_w, _ = strain_paint.gating_info["auto"][0]
        lidx, ridx = strain_paint.get_gate_indices(gate_t, gate_w)
        # Without post-gate highpass, inpainted interior should stay smooth/finite
        self.assertTrue(np.all(np.isfinite(strain_paint[lidx:ridx].numpy())))
        self.assertGreater(len(make_frequency_series(strain_paint)), 0)

    def test_multi_ifo_paint_and_inpaint_cli_parsing(self):
        """Multi-IFO CLI parser accepts IFO:METHOD --paint-method and --inpaint-edges."""
        parser = argparse.ArgumentParser()
        insert_strain_option_group_multi_ifo(parser)
        StrainSegments.insert_segment_option_group_multi_ifo(parser)
        args = parser.parse_args(
            [
                "--paint-method",
                "H1:toeplitz",
                "L1:cholesky",
                "--paint-ridge",
                "H1:0.001",
                "L1:0.002",
                "--inpaint-edges",
            ]
        )
        self.assertEqual(args.paint_method["H1"], "toeplitz")
        self.assertEqual(args.paint_method["L1"], "cholesky")
        self.assertEqual(args.paint_ridge["H1"], 0.001)
        self.assertEqual(args.paint_ridge["L1"], 0.002)
        self.assertTrue(args.inpaint_edges)


if __name__ == "__main__":
    unittest.main()
