import os

import numpy as np

from pycbc.fft import FFT, IFFT
from pycbc.filter.matchedfilter import (
    correlate, get_cutoff_indices, matched_filter_core,
)
from pycbc.filter.matchedfilter_cpu import (
    fast_multiply_analytic_cython,
    find_peaks_in_block_cython,
)
from pycbc.types import Array, complex64, zeros

# Optional matchedfilter back end for the innermost product/inverse/peak step.
# Everything else -- the reference matched filter, the filter bank FFTs, the
# per-block forward FFT -- stays on pycbc's own FFT.
try:
    import matchedfilter as _mf
except ImportError:
    _mf = None


def matchedfilter_available():
    """True if the matchedfilter back end can be used."""
    return _mf is not None


class MatchedFilterRatioControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    Uses pycbc's class-based FFT interface (pycbc.fft.FFT/IFFT) for all FFT
    operations, so the actual backend (MKL, FFTW, numpy, ...) follows
    whatever --processing-scheme the caller selected, with FFT/IFFT plans
    built once per batch size and reused.

    Note: pycbc's class-based IFFT does not normalize by 1/N the way
    numpy.fft.ifft does; see _forward_block_fft for where that's corrected.
    """

    def __init__(self, snr_threshold, delta_f,
                 high_frequency_cutoff=None, fir_fft_length=4096, batch_size=64,
                 tap_sample_rate=2048, engine_sample_rate=2048,
                 engine='pycbc', false_dismissal=1e-3, coarse_band_hz=0,
                 first_stage_snr=0):
        self.delta_f = delta_f
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff

        self.tap_sr = int(tap_sample_rate)
        self.engine_sr = int(engine_sample_rate)

        self.threshold_sq = float(snr_threshold**2)

        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size

        # Taps are generated at tap_sample_rate but filtered against data at
        # engine_sample_rate; the ratio must be an exact integer so the
        # high-resolution FFT below preserves delta_f.
        exact_ratio = self.tap_sr / self.engine_sr
        self.decimation_factor = int(np.round(exact_ratio))
        if abs(exact_ratio - self.decimation_factor) > 1e-5 or self.decimation_factor < 1:
            raise ValueError(
                f"Multi-rate Error: The bank sample rate ({self.tap_sr} Hz) must be "
                f"an exact integer multiple of the engine sample "
                f"rate ({self.engine_sr} Hz).\n"
                f"Calculated ratio was {exact_ratio:.4f}. Please use an engine "
                f"sample rate that evenly divides the bank sample rate."
            )
        # Scaling by decimation_factor keeps delta_f matched between a
        # high_res_fft_len buffer at the tap rate and a fir_fft_len buffer at
        # the engine rate, so keeping only the first fir_fft_len bins below
        # is a decimation rather than a change in frequency resolution.
        self.high_res_fft_len = self.fir_fft_len * self.decimation_factor

        # Plans are keyed by (nbatch, size) and cached for the engine's
        # lifetime: n_filters is not generally a multiple of batch_size, so a
        # few smaller plans get cached alongside the dominant batch size.
        self._fft_plans = {}
        self._ifft_plans = {}
        self._ref_plan = None
        self._ref_direct = os.environ.get('PYCBC_RATIO_REFDIRECT', '1') != '0' 

        # matchedfilter back end.  "hier" gates on a low-band coarse pass and only
        # reconstructs where a detection is still possible; "flat" is the same
        # correlation with no gate, which exists to separate "matchedfilter is faster"
        # from "the gate is working".  Off by default.
        # Back end for the innermost correlate/inverse/peak step.  The
        # environment variables remain as an override, which is convenient when
        # bisecting a behaviour change without re-plumbing a workflow, but the
        # command line is what a run should set.
        mode = {'pycbc': '', 'matchedfilter': 'flat',
                'matchedfilter-hierarchical': 'hier'}.get(engine, engine)
        self._engine_mode = os.environ.get('PYCBC_RATIO_ENGINE', mode).lower()
        if self._engine_mode and _mf is None:
            raise ImportError(
                "ratio-filter-engine=%s needs matchedfilter, which is not installed"
                % engine)
        self._engine_fd = float(os.environ.get(
            'PYCBC_RATIO_ENGINE_FD', false_dismissal))
        # The band is given in Hz, which is what a reader of this search
        # thinks in; matchedfilter wants bins.  delta_f is set by the ratio
        # filter's own transform, not by the data segments.
        band_hz = float(os.environ.get('PYCBC_RATIO_BAND_HZ', coarse_band_hz))
        delta_f = self.engine_sr / float(self.fir_fft_len)
        nyquist = self.engine_sr / 2.0
        if band_hz and band_hz >= nyquist:
            raise ValueError(
                "--ratio-filter-band %g Hz is at or above the ratio filter's "
                "Nyquist frequency (%g Hz); that filters the whole band and "
                "is not a first pass" % (band_hz, nyquist))
        self._coarse_band = int(round(band_hz / delta_f)) if band_hz else 0
        self._coarse_band_hz = band_hz
        # Which device the engine's plans run on.  Probative: the ratio
        # filter was written CPU-only, so this exists to find out what the
        # GPU path actually costs END TO END -- including the host transfers
        # that the library's own benchmarks deliberately exclude, and which
        # here are paid once per filter block rather than once per batch.
        self._device = os.environ.get('PYCBC_RATIO_DEVICE') or None
        self._first_stage_snr = float(os.environ.get(
            'PYCBC_RATIO_FIRST_STAGE_SNR', first_stage_snr))
        self._ap_ref = None

    def _get_plan(self, plans, cls, nbatch, size):
        key = (nbatch, size)
        cached = plans.get(key)
        if cached is None:
            in_arr = zeros(nbatch * size, dtype=complex64)
            out_arr = zeros(nbatch * size, dtype=complex64)
            if nbatch > 1:
                plan = cls(in_arr, out_arr, nbatch=nbatch, size=size)
            else:
                plan = cls(in_arr, out_arr)
            in_view = in_arr.data.reshape(nbatch, size)
            out_view = out_arr.data.reshape(nbatch, size)
            cached = plans[key] = (plan, in_view, out_view)
        return cached

    def _get_fft_plan(self, nbatch, size):
        return self._get_plan(self._fft_plans, FFT, nbatch, size)

    def _get_ifft_plan(self, nbatch, size):
        return self._get_plan(self._ifft_plans, IFFT, nbatch, size)

    def prepare_filters(self, fir_taps, tap_counts):
        """
        Prepare frequency-domain filters for a batch of taps.
        """
        if self._engine_mode and hasattr(_mf, 'TimeDomainFilterBank'):
            self._td_bank = _mf.TimeDomainFilterBank(
                fir_taps, tap_counts,
                tap_sample_rate=self.tap_sr,
                data_sample_rate=self.engine_sr,
                engine=self._engine_mode,
                threshold=self.snr_threshold,
                false_dismissal=self._engine_fd,
                first_stage_snr=self._first_stage_snr,
                coarse_band_hz=self._coarse_band_hz,
                device=self._device,
                max_batch_size=self.batch_size,
            )
            self._tap_sizes = np.array([g['n'] for g in self._td_bank.groups])
            n_taps_max = int(np.max(tap_counts))
            self._ap_ref = None
            return self._td_bank.filters_f, n_taps_max

        n_filters, n_taps = fir_taps.shape
        if n_taps >= self.fir_fft_len:
             raise ValueError("FIR Taps (%d) exceed FFT block length (%d)" %
                              (n_taps, self.fir_fft_len))

        n_taps_max = int(np.max(tap_counts))

        # Overlap-save group boundaries.  These depend only on the bank, so
        # computing them here rather than per segment takes a quantile over
        # the whole tap-count array off the hot path.
        tap_groups = 3
        self._tap_sizes = np.quantile(
            tap_counts, np.linspace(0, 1, tap_groups + 1)[1:]).astype(int)

        filters_f = self._fft_all_filters(fir_taps, tap_counts)
        return filters_f, n_taps_max

    def make_reference_series(self, stilde, psd, ref_template):
        """Return the scaled reference SNR series and its template sigmasq."""
        h_norm = ref_template.sigmasq(psd)
        if self._ref_direct:
            snr, norm = self._reference_snr(stilde, psd, ref_template, h_norm)
        else:
            snr, _, norm = matched_filter_core(
                ref_template, stilde, psd=psd,
                low_frequency_cutoff=ref_template.f_lower,
                high_frequency_cutoff=self.f_high, h_norm=h_norm)
        scale = (norm * stilde.delta_t) / self.decimation_factor
        # The cached IFFT output is overwritten by the next reference; upper
        # series must survive while all of its middle children are processed.
        return np.asarray(snr.numpy(), dtype=np.complex64).copy() * scale, h_norm

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, indices,
                        valid_slice=None, reference_series=None,
                        profile_template=None):
        """
        Process a single data segment.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        if isinstance(ref_template, (int, float, np.floating)):
            h_norm = float(ref_template)
        elif getattr(self, '_cached_hnorm_key', None) == (id(ref_template), id(psd)):
            h_norm = self._cached_hnorm
        else:
            h_norm = self._cached_hnorm = ref_template.sigmasq(psd)
            self._cached_hnorm_key = (id(ref_template), id(psd))

        if self._engine_mode in ('hier', 'check') and self._ap_ref is None:
            profile = profile_template if profile_template is not None else ref_template
            if hasattr(self, '_td_bank'):
                self._td_bank.set_reference_from_template(stilde, psd, profile, f_high=self.f_high)
                self._ap_ref = True

        if reference_series is None:
            if self._ref_direct:
                snr, norm = self._reference_snr(stilde, psd, ref_template, h_norm)
            else:
                snr, _, norm = matched_filter_core(
                    ref_template, stilde, psd=psd,
                    low_frequency_cutoff=ref_template.f_lower,
                    high_frequency_cutoff=self.f_high, h_norm=h_norm)
            self.ref_snr = snr.numpy()
            self.ref_snr *= (norm * stilde.delta_t) / self.decimation_factor
        else:
            self.ref_snr = reference_series if isinstance(reference_series, np.ndarray) and reference_series.dtype == np.complex64 else np.asarray(reference_series, dtype=np.complex64)
            if len(self.ref_snr) != (len(stilde) - 1) * 2:
                raise ValueError('reference series length does not match data segment')

        if hasattr(self, '_td_bank'):
            res = self._td_bank.filter_series(self.ref_snr, valid_slice=valid_slice)
            local_idxs = res.template_indices
            t_idxs = res.sample_indices
            snr_vals = res.snr
            tstarts = res.block_starts
        else:
            local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
                self.ref_snr, filters_f, n_taps, valid_slice
            )

        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, tstarts, h_norm
        else:
            return [], [], [], tstarts, h_norm

    def _fft_all_filters(self, taps, counts):
        """Helper to FFT all filters, batch_size rows at a time."""
        n_filters, n_taps_alloc = taps.shape
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)

        high_res_fft_len = self.high_res_fft_len

        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start

            plan, in_view, out_view = self._get_fft_plan(batch_len, high_res_fft_len)

            in_view[:] = 0.0
            tmp_taps = taps[start:end]
            in_view[:, :n_taps_alloc] = tmp_taps

            # Roll each row so its center tap sits at index 0, with earlier
            # taps wrapping to the end -- the circular layout an FFT-based
            # FIR filter needs. get_fd_fir in waveform/bank.py undoes this
            # same roll when reconstructing a single template's time-domain
            # filter, so the two must stay in sync.
            current_counts = counts[start:end]
            roll_offsets = -(current_counts // 2)

            cols_high = np.arange(high_res_fft_len)
            rows = np.arange(batch_len)[:, None]
            shifted_cols_high = (cols_high[None, :] - roll_offsets[:, None]) % high_res_fft_len

            current_data = in_view.copy()
            in_view[:] = current_data[rows, shifted_cols_high]

            plan.execute()

            # high_res_fft_len spans the bank's full tap sample rate; only
            # the first fir_fft_len bins fall within the engine's decimated
            # frequency range, so slicing to them is the decimation step.
            fft_sliced = out_view[:, :self.fir_fft_len]

            # Conjugate so multiplying against the data spectrum and
            # inverse-transforming (_execute_blocked_kernel) yields a
            # correlation rather than a convolution.
            filters_f[start:end] = np.conj(fft_sliced)

        return filters_f

    def _reference_snr(self, stilde, psd, ref_template, h_norm):
        """The reference SNR series, driving the class-based IFFT directly.

        matched_filter_core goes through pycbc.fft's *function* API, which on
        the MKL backend builds a DFTI descriptor, uses it once and frees it on
        every call -- 2 ms of a 3 ms transform at 2^20.  The class API caches
        the descriptor, which is what it is for, and the plan cache this object
        already keeps for the block FFTs serves here too.  The buffers are
        cached with it, so the four per-call allocations go as well.

        Same arithmetic as matched_filter_core, in the same order.
        """
        N = (len(stilde) - 1) * 2
        kmin, kmax = get_cutoff_indices(
            ref_template.f_lower, self.f_high, stilde.delta_f, N)
        plan, qt, q = self._get_ref_plan(N, kmax)
        # Only the low strip is cleared per call.  correlate writes solely
        # [kmin:kmax] and execute() reads qt and writes q, so nothing ever
        # dirties the region above kmax -- and kmax derives from self.f_high,
        # an engine attribute, so it is fixed for the whole run.  Clearing it
        # once with the buffer removes 4 MB of stores per segment.  kmin comes
        # from the template's f_lower and does vary, so that strip stays.
        qt[:kmin] = 0
        # stilde arrives overwhitened (see pycbc_inspiral_fir), so there is
        # no PSD division on this path at all.
        correlate(ref_template[kmin:kmax], stilde[kmin:kmax], qt[kmin:kmax])
        plan.execute()
        norm = (4.0 * stilde.delta_f) / np.sqrt(h_norm)
        return q, norm

    def _get_ref_plan(self, size, kmax):
        """Cached IFFT plan and its buffers for the reference filter.

        Keyed on kmax as well as size: zeros() hands back a cleared buffer,
        which is what lets _reference_snr skip clearing above kmax on every
        call.  If the upper cutoff ever did change, this reallocates rather
        than leaving a stale band in the transform.
        """
        cached = self._ref_plan
        if cached is None or cached[0] != (size, kmax):
            qt = zeros(size, dtype=complex64)
            q = zeros(size, dtype=complex64)
            plan = IFFT(qt, q)
            cached = self._ref_plan = ((size, kmax), plan, qt, q)
        _, plan, qt, q = cached
        return plan, Array(qt, copy=False), q

    def _forward_block_fft(self, segment):
        """FFT one time-domain block (nbatch=1), pre-dividing by fir_fft_len.

        pycbc's class-based IFFT (unlike numpy.fft.ifft) does
        not normalize its output by 1/N. Since this block spectrum only
        ever gets multiplied by a filter spectrum and then inverse
        transformed (in _execute_blocked_kernel), pre-dividing it here by N
        once -- rather than dividing the batched product after every IFFT --
        gives an already-correctly-normalized correlation for free:
        ifft_raw(block_f/N * filter_f) == ifft_normalized(block_f * filter_f).
        """
        plan, in_view, out_view = self._get_fft_plan(1, self.fir_fft_len)
        in_view[0, :] = 0.0
        in_view[0, :len(segment)] = segment
        plan.execute()
        return out_view[0].copy() / self.fir_fft_len

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice):
        """Inner loop: Time-Blocking + Filter-Batching for stock pycbc engine."""
        nsizes = self._tap_sizes
        n_samples = len(data)
        n_filters = len(filters_f)
        N_FFT = self.fir_fft_len

        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        all_tstarts = []

        if valid_slice:
            v_start = valid_slice.start
            v_stop = valid_slice.stop
        else:
            v_start = 0
            v_stop = n_samples

        block_f_cache = {}
        for f_start in range(0, n_filters, self.batch_size):
            f_end = min(f_start + self.batch_size, n_filters)
            actual_batch_size = f_end - f_start

            ifft_plan, current_mult_view, current_corr_view = self._get_ifft_plan(
                actual_batch_size, N_FFT
            )

            # Overlap-save: valid output samples per block.
            n_taps_max = int(n_taps) if np.isscalar(n_taps) else int(n_taps[f_start:f_end].max())
            i = np.searchsorted(nsizes, n_taps_max)
            n_taps_max = nsizes[i]

            N_VALID = N_FFT - n_taps_max + 1
            STEP = N_VALID
            bad_start = n_taps_max // 2

            first_block_idx = (v_start - bad_start) // STEP
            loop_start = first_block_idx * STEP

            for t_start in range(loop_start, n_samples, STEP):
                block_valid_t0 = t_start + bad_start
                if block_valid_t0 >= v_stop:
                    break
                if block_valid_t0 + N_VALID <= v_start:
                    continue
                roi_start = max(v_start, block_valid_t0)
                roi_stop = min(v_stop, block_valid_t0 + N_VALID)
                roi_len = roi_stop - roi_start
                if roi_len <= 0:
                    continue

                buf_slice_start = roi_start - t_start
                t_end = min(t_start + N_FFT, n_samples)
                if t_start not in block_f_cache:
                    block_f_cache[t_start] = self._forward_block_fft(
                        data[t_start:t_end])

                block_f_view = block_f_cache[t_start]
                filter_batch_f = filters_f[f_start:f_end]

                fast_multiply_analytic_cython(
                    block_f_view, filter_batch_f, current_mult_view
                )
                ifft_plan.execute()
                f_list, t_list, s_list = find_peaks_in_block_cython(
                    current_corr_view, roi_start, roi_len, self.threshold_sq,
                    f_start, input_offset=buf_slice_start
                )
                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)
                    all_tstarts.extend([t_start] * len(s_list))

        return (np.array(all_f_idxs, dtype=np.int32),
                np.array(all_t_idxs, dtype=np.int64),
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
