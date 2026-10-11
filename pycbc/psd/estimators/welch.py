# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Segmented High-Resolution Welch PSD Model (Baseline / Control)."""

import numpy as np
from scipy.interpolate import interp1d
from pycbc.types import TimeSeries, FrequencySeries
import pycbc.psd
from .base import BasePSDModel


class SegmentedWelchPSDModel(BasePSDModel):
    """High-resolution Welch median baseline model across sliding overlapping intervals."""

    def __init__(self, f_low=20.0, f_high=None, seg_len_s=16.0, seg_stride_s=8.0, avg_method='median'):
        super().__init__(f_low=f_low, f_high=f_high)
        self.seg_len_s = float(seg_len_s)
        self.seg_stride_s = float(seg_stride_s)
        self.avg_method = str(avg_method)
        self.welch_psd = None
        self.freqs = None

    def fit(self, timeseries: TimeSeries, f_low: float = None, f_high: float = None):
        if f_low is not None:
            self.f_low = float(f_low)
        self.fs = float(timeseries.sample_rate)
        if f_high is not None:
            self.f_high = float(f_high)
        elif self.f_high is None:
            self.f_high = 0.9 * (self.fs / 2.0)

        self.t_start = float(timeseries.start_time)
        self.t_end = float(timeseries.end_time)

        N_seg = int(round(self.seg_len_s * self.fs))
        N_stride = int(round(self.seg_stride_s * self.fs))
        N_seg = min(N_seg, len(timeseries) // 2)

        psd_fs = pycbc.psd.welch(
            timeseries, seg_len=N_seg, seg_stride=N_stride, avg_method=self.avg_method
        )
        self.welch_psd = psd_fs.numpy().astype(np.float64)
        self.freqs = psd_fs.sample_frequencies.numpy().astype(np.float64)

    def eval_psd(self, freqs: np.ndarray, t: float = None) -> np.ndarray:
        freqs = np.asarray(freqs, dtype=np.float64)
        f_interp = interp1d(self.freqs, self.welch_psd, kind='linear',
                            bounds_error=False, fill_value=(self.welch_psd[0], self.welch_psd[-1]))
        out_psd = f_interp(freqs)

        # Low-frequency brick wall below f_low, high-frequency cutoff above f_high
        out_psd = np.asarray(out_psd, dtype=np.float64)
        out_psd[freqs < self.f_low] = 1e30
        if self.f_high is not None:
            out_psd[freqs > self.f_high] = 1e30
        return out_psd
