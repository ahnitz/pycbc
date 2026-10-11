# Copyright (C) 2026 PyCBC Team
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""Adaptive Thomson Multitaper Spectrogram PSD Model."""

import numpy as np
from scipy.signal.windows import dpss
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.ndimage import median_filter, gaussian_filter1d
from pycbc.types import TimeSeries, FrequencySeries
from .base import BasePSDModel


class MultitaperPSDModel(BasePSDModel):
    """Adaptive Thomson Slepian Multitaper Spectrogram model.

    Computes orthogonal DPSS eigenspectra in sliding windows across the job duration
    to minimize spectral leakage at seismic walls and sharp line resonances.
    """

    def __init__(self, f_low=20.0, f_high=None, win_len_s=8.0, step_s=2.0, NW=3.0, Kmax=5):
        super().__init__(f_low=f_low, f_high=f_high)
        self.win_len_s = float(win_len_s)
        self.step_s = float(step_s)
        self.NW = float(NW)
        self.Kmax = int(Kmax)
        self.spec_times = None
        self.spec_freqs = None
        self.spectrogram = None
        self.baseline_psd = None

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

        N_win = int(round(self.win_len_s * self.fs))
        N_step = int(round(self.step_s * self.fs))

        # Generate Kmax DPSS orthogonal tapers
        tapers, ratios = dpss(N_win, self.NW, self.Kmax, return_ratios=True)
        # Normalize each taper to unit discrete energy
        for k in range(self.Kmax):
            tapers[k] /= np.sqrt(np.sum(tapers[k] ** 2))

        x = timeseries.numpy().astype(np.float64)
        N_total = len(x)
        t_centers = []
        specs = []

        f_rfft = np.fft.rfftfreq(N_win, 1.0 / self.fs)
        f_mask = (f_rfft >= max(5.0, self.f_low - 5.0)) & (f_rfft <= self.f_high)
        sub_freqs = f_rfft[f_mask]

        scale_factor = 2.0 / self.fs

        for s_idx in range(0, N_total - N_win + 1, N_step):
            e_idx = s_idx + N_win
            chunk = x[s_idx:e_idx]
            t_mid = self.t_start + (s_idx + e_idx) * 0.5 / self.fs
            t_centers.append(t_mid)

            # Thomson Multitaper average across orthogonal eigenspectra
            mt_spec = np.zeros(len(sub_freqs), dtype=np.float64)
            for k in range(self.Kmax):
                tapered_data = tapers[k] * chunk
                fft_val = np.fft.rfft(tapered_data)[f_mask]
                mt_spec += scale_factor * (np.abs(fft_val) ** 2)

            mt_spec /= self.Kmax
            specs.append(mt_spec)

        t_centers = np.array(t_centers)
        specs = np.array(specs)  # shape: (n_times, n_freqs)

        # Baseline: time-median spectrum across all multitaper windows
        base = np.median(specs, axis=0)

        # Smooth spectrogram in time to reject transient glitches
        # Median filter over ~15 seconds across time
        k_time = max(3, int(round(15.0 / self.step_s)))
        if k_time % 2 == 0:
            k_time += 1

        smooth_specs = np.zeros_like(specs)
        for fi in range(len(sub_freqs)):
            smooth_specs[:, fi] = median_filter(specs[:, fi], size=k_time)

        self.spec_times = t_centers
        self.spec_freqs = sub_freqs
        self.spectrogram = smooth_specs
        self.baseline_psd = base

    def eval_psd(self, freqs: np.ndarray, t: float = None) -> np.ndarray:
        freqs = np.asarray(freqs, dtype=np.float64)
        if t is None or self.spectrogram is None or len(self.spec_times) <= 1:
            # Baseline interpolation
            f_interp = interp1d(self.spec_freqs, self.baseline_psd, kind='linear',
                                bounds_error=False, fill_value=(self.baseline_psd[0], self.baseline_psd[-1]))
            out_psd = f_interp(freqs)
        else:
            # Find nearest time index
            t_idx = int(np.clip(np.argmin(np.abs(self.spec_times - t)), 0, len(self.spec_times) - 1))
            spec_t = self.spectrogram[t_idx, :]
            f_interp = interp1d(self.spec_freqs, spec_t, kind='linear',
                                bounds_error=False, fill_value=(spec_t[0], spec_t[-1]))
            out_psd = f_interp(freqs)

        # Apply low-frequency brick wall below f_low
        out_psd = np.asarray(out_psd, dtype=np.float64)
        out_psd[freqs < self.f_low] = 1e30
        return out_psd
