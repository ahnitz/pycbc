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

"""Unit tests for in-memory regularized overwhitening, overlap-save streaming
convolution, and coherent line subtraction.
"""

import unittest

import numpy as np

import pycbc.filter
import pycbc.psd
import pycbc.waveform
from pycbc.types import TimeSeries


class TestOverwhiten(unittest.TestCase):
    def setUp(self):
        self.fs = 1024.0
        self.duration = 64.0
        self.N = int(self.fs * self.duration)
        rng = np.random.default_rng(12345)
        white_noise = rng.normal(0, 1, self.N)
        freqs = np.fft.rfftfreq(self.N, 1.0 / self.fs)
        color = 1.0 / np.sqrt(1.0 + (freqs / 50.0) ** 2)
        colored_fd = np.fft.rfft(white_noise) * color
        colored_td = np.fft.irfft(colored_fd)
        self.ts = TimeSeries(colored_td, delta_t=1.0 / self.fs, epoch=1000000000.0)

    def test_construct_regularized_kernels(self):
        psd = pycbc.psd.welch(
            self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs)
        )
        invpsd, w_kernel, fir_dow, fir_w, psd_ist, df = (
            pycbc.filter.construct_regularized_kernels(
                psd,
                self.fs,
                self.duration,
                f_low=18.0,
                f_taper=4.0,
                max_filter_duration=2.0,
            )
        )
        self.assertGreater(len(invpsd), 0)
        self.assertGreater(len(w_kernel), 0)
        self.assertEqual(len(fir_dow) % 2, 1)
        self.assertEqual(len(fir_w) % 2, 1)
        k_14 = int(round(13.0 / df))
        self.assertAlmostEqual(invpsd[k_14], 0.0)
        self.assertAlmostEqual(w_kernel[k_14], 0.0)

    def test_overwhiten_and_whitening_normality(self):
        ts_dow, ts_dw, meta = pycbc.filter.overwhiten_strain(
            self.ts,
            f_low=18.0,
            f_taper=4.0,
            max_filter_duration=2.0,
            use_overlap_save=True,
            chunk_len_s=8.0,
        )
        self.assertEqual(len(ts_dow), len(self.ts))
        self.assertEqual(len(ts_dw), len(self.ts))
        edge = int(4.0 * self.fs)
        interior_dw = ts_dw.numpy()[edge : len(ts_dw) - edge]
        self.assertAlmostEqual(np.mean(interior_dw), 0.0, delta=0.1)
        self.assertAlmostEqual(np.std(interior_dw), 1.0, delta=0.2)

    def test_overlap_save_equivalence(self):
        psd = pycbc.psd.welch(
            self.ts, seg_len=int(4 * self.fs), seg_stride=int(2 * self.fs)
        )
        invpsd, w_kernel, fir_dow, fir_w, psd_ist, df = (
            pycbc.filter.construct_regularized_kernels(
                psd,
                self.fs,
                self.duration,
                f_low=18.0,
                f_taper=4.0,
                max_filter_duration=2.0,
            )
        )
        dow_stream, dw_stream = pycbc.filter.produce_conditioned_streams(
            self.ts,
            invpsd,
            w_kernel,
            fir_dow,
            fir_w,
            delta_f_full=df,
            use_overlap_save=True,
            chunk_len_s=8.0,
        )
        dow_batch, dw_batch = pycbc.filter.produce_conditioned_streams(
            self.ts,
            invpsd,
            w_kernel,
            fir_dow,
            fir_w,
            delta_f_full=df,
            use_overlap_save=False,
        )
        edge = int(4.0 * self.fs)
        diff = np.abs(dow_stream.numpy()[edge:-edge] - dow_batch.numpy()[edge:-edge])
        rel_diff = np.max(diff) / np.std(dow_batch.numpy()[edge:-edge])
        self.assertLess(rel_diff, 1e-3)

    def test_subtract_coherent_lines(self):
        fs = 1024.0
        duration = 64.0
        N = int(fs * duration)
        dt = 1.0 / fs
        t = np.arange(N) * dt
        rng = np.random.default_rng(1234)
        noise = rng.normal(0, 1.0, N)
        ts_noise = TimeSeries(noise, delta_t=dt, epoch=1000000000.0)

        hp, _ = pycbc.waveform.get_fd_waveform(
            approximant="IMRPhenomD",
            mass1=20,
            mass2=20,
            f_lower=20.0,
            delta_f=ts_noise.delta_f,
            distance=200,
        )
        hp.resize(len(ts_noise.to_frequencyseries()))
        hp_td = pycbc.waveform.apply_fd_time_shift(hp, duration / 2.0).to_timeseries()
        hp_td.start_time = ts_noise.start_time

        line_60 = 20.0 * np.sin(2.0 * np.pi * 60.0 * t)
        ts_with_line = (
            ts_noise
            + TimeSeries(line_60, delta_t=dt, epoch=ts_noise.start_time)
            + hp_td
        )

        ts_clean, line_list = pycbc.filter.subtract_coherent_lines(
            ts_with_line, f_low=18.0, line_threshold=3.0, max_lines=16
        )

        sub_freqs = [l[0] for l in line_list]
        self.assertTrue(any(abs(f - 60.0) < 1.0 for f in sub_freqs))

        df = 1.0 / duration
        idx_60 = int(round(60.0 / df))
        pow_orig = np.abs(np.fft.rfft(ts_with_line.numpy())[idx_60]) ** 2
        pow_clean = np.abs(np.fft.rfft(ts_clean.numpy())[idx_60]) ** 2
        self.assertGreater(1.0 - (pow_clean / pow_orig), 0.90)

        clean_waveform, _ = pycbc.filter.subtract_coherent_lines(
            hp_td, lines=[(60.0, 10.0, 0.2)]
        )
        overlap = pycbc.filter.match(hp_td, clean_waveform)[0]
        self.assertGreater(overlap, 0.999)

    def test_edge_inpainting_eliminates_boundary_explosion(self):
        fs = 2048.0
        dur = 32.0
        N = int(fs * dur)
        rng = np.random.default_rng(42)
        white = rng.normal(0, 1, N)
        freqs = np.fft.rfftfreq(N, 1.0 / fs)
        color = 1.0 / np.sqrt(1.0 + (np.maximum(freqs, 20.0) / 40.0) ** 4)
        colored = np.fft.irfft(np.fft.rfft(white) * color)
        ts = TimeSeries(colored, delta_t=1.0 / fs, epoch=1000000000.0)

        dow_inp, dw_inp, meta = pycbc.filter.overwhiten_strain(
            ts,
            f_low=20.0,
            max_filter_duration=2.0,
            inpaint_edges=True,
            edge_pad_duration=2.0,
        )
        self.assertEqual(len(dow_inp), len(ts))
        self.assertEqual(float(dow_inp.start_time), float(ts.start_time))
        self.assertTrue(meta["inpainted_edges"])

        edge_n = int(0.25 * fs)
        std_start = float(np.std(dw_inp.numpy()[:edge_n]))
        std_end = float(np.std(dw_inp.numpy()[-edge_n:]))
        std_int = float(np.std(dw_inp.numpy()[int(2.0 * fs) : -int(2.0 * fs)]))

        self.assertLess(std_start, 3.0)
        self.assertLess(std_end, 3.0)
        self.assertAlmostEqual(std_int, 1.0, delta=0.15)


if __name__ == "__main__":
    unittest.main()
