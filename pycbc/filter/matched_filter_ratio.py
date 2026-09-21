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
from pycbc.fft import FFT, IFFT
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
                 engine='pycbc', false_dismissal=1e-3, coarse_band=0):
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
        self._coarse_band = int(os.environ.get(
            'PYCBC_RATIO_BAND', coarse_band))
        self._engine_plans = {}
        self._ap_ref = None
        self._ap_loaded = None
        self._ap_filters = None
        self._ap_ref_id = {}
        # One block per call.  matchedfilter batches D x T and larger D looks better
        # in isolation, but there the same data array is reused and stays
        # cache-warm; here every call brings a fresh 2 MB block spectrum, so the
        # ingest is the cost and batching only adds a copy on top.  Measured:
        # D=1 0.178 s, D=8 0.257, D=64 0.260.
        self._ap_ndata = int(os.environ.get('PYCBC_RATIO_ENGINE_NDATA', '1'))
        self._ap_series = os.environ.get('PYCBC_RATIO_SERIES', '1') != '0'
        self._chk_tot = 0
        self._chk_miss = 0
        self._chk_missed_snr = []
        self._chk_detail = []

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

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, indices,
                        valid_slice=None):
        """
        Process a single data segment.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        import time as _tm
        _t0 = _tm.perf_counter()
        h_norm = ref_template.sigmasq(psd)

        if self._ref_direct:
            snr, norm = self._reference_snr(stilde, psd, ref_template, h_norm)
        else:
            snr, _, norm = matched_filter_core(
                ref_template, stilde, psd=psd,
                low_frequency_cutoff=ref_template.f_lower,
                high_frequency_cutoff=self.f_high,
                h_norm=h_norm
            )

        # Scale in place.  Written as snr.numpy() * a / b this built two
        # full-length temporaries per segment -- 8 MB each at 2^20 -- for a
        # scaling by two constants.  q is the cached IFFT output buffer and is
        # fully overwritten by the next plan.execute(), so mutating it here is
        # safe.  decimation_factor is computed and validated in __init__.
        self.ref_snr = snr.numpy()
        self.ref_snr *= (norm * stilde.delta_t) / self.decimation_factor
        _t1 = _tm.perf_counter()

        local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
            self.ref_snr, filters_f, n_taps, valid_slice
        )
        _t2 = _tm.perf_counter()
        self._ph_ref = getattr(self, '_ph_ref', 0.0) + (_t1 - _t0)
        self._ph_ker = getattr(self, '_ph_ker', 0.0) + (_t2 - _t1)
        if os.environ.get('PYCBC_RATIO_PHASE'):
            import sys
            tot = self._ph_ref + self._ph_ker
            extra = ""
            if self._engine_mode == 'hier' and self._engine_plans:
                pr = tg = 0
                cfg = None
                for pl in self._engine_plans.values():
                    p_, t_ = pl.stats
                    pr += p_; tg += t_
                    cfg = pl.config
                if pr:
                    extra = ("  gate: %.1f%% of %d pair-calls triggered, band=%d U=%d K=%d"
                             % (100.0 * tg / pr, pr, cfg[0], cfg[1], cfg[2]))
            print("[ratio-phase] mode=%-5s reference-MF %.3f s (%.0f%%)  "
                  "ratio-kernel %.3f s (%.0f%%)%s"
                  % (self._engine_mode or 'stock', self._ph_ref,
                     100 * self._ph_ref / tot, self._ph_ker,
                     100 * self._ph_ker / tot, extra), file=sys.stderr)

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
        plan, qt, q = self._get_ref_plan(N)
        qt[:kmin] = 0
        qt[kmax:] = 0
        # stilde arrives overwhitened (see pycbc_inspiral_fir), so there is
        # no PSD division on this path at all.
        correlate(ref_template[kmin:kmax], stilde[kmin:kmax], qt[kmin:kmax])
        plan.execute()
        norm = (4.0 * stilde.delta_f) / np.sqrt(h_norm)
        return q, norm

    def _get_ref_plan(self, size):
        """Cached IFFT plan and its buffers for the reference filter."""
        cached = self._ref_plan
        if cached is None or cached[0] != size:
            qt = zeros(size, dtype=complex64)
            q = zeros(size, dtype=complex64)
            plan = IFFT(qt, q)
            cached = self._ref_plan = (size, plan, qt, q)
        _, plan, qt, q = cached
        return plan, Array(qt, copy=False), q

    def _set_engine_reference(self, data):
        """Reference SNR distribution for the gate.

        matchedfilter's gate needs the fraction of SNR below its band edge.  Left to
        itself it measures that from the template, which here is a short
        broadband ratio filter -- the wrong distribution entirely, since the
        filter's output reconstructs the fine template's strongly
        low-frequency SNR.  Measured on a real bank: 0.30 of the filter's own
        power sits below 256 Hz against 0.927 of the SNR it produces.

        The right distribution is the reference SNR series' own power spectrum,
        which is exactly what this data is, and it is shared by every template
        in the group -- so it is set once per segment rather than per template.
        """
        n = self.fir_fft_len
        nblk = max(1, min(8, len(data) // n))
        acc = np.zeros(n, dtype=np.float64)
        for b in range(nblk):
            seg = data[b * n:(b + 1) * n]
            if len(seg) < n:
                break
            acc += np.abs(np.fft.fft(seg.astype(np.complex64))) ** 2
        acc[n // 2 + 1:] = 0.0        # analytic: the kernel uses only these bins
        if acc.sum() <= 0:
            return None
        return (acc / acc.sum()).astype(np.float32)

    def _get_engine_plan(self, nbatch, ndata=1):
        """One plan per batch width, templates reloaded per filter batch."""
        plan = self._engine_plans.get((nbatch, ndata))
        if plan is None:
            snr = float(self.snr_threshold)
            if self._engine_mode in ('hier','check'):
                # Band override, for checking the design table's choice against
                # measurement rather than trusting it.
                bd = self._coarse_band
                if bd:
                    plan = _mf.HierarchicalFilter(
                        self.fir_fft_len, ndata=ndata, ntemplates=nbatch,
                        snr=snr, fd=self._engine_fd, band=bd,
                        oversample=2, taps=8)
                else:
                    plan = _mf.HierarchicalFilter(
                        self.fir_fft_len, ndata=ndata, ntemplates=nbatch,
                        snr=snr, fd=self._engine_fd)
            else:
                plan = _mf.MatchedFilter(
                    self.fir_fft_len, ndata=ndata, ntemplates=nbatch)
            self._engine_plans[(nbatch, ndata)] = plan
        # Only when the reference actually changes.  set_reference() re-measures
        # the recovery factors over noise realisations, which is cheap once per
        # segment and ruinous once per filter batch -- doing the latter cost
        # ~67 us per call and hid most of the gate's benefit.
        if (self._engine_mode in ('hier', 'check') and self._ap_ref is not None
                and self._ap_ref_id.get(id(plan)) is not self._ap_ref):
            plan.set_reference(self._ap_ref)
            self._ap_ref_id[id(plan)] = self._ap_ref
        return plan

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
        """
        Inner loop: Time-Blocking + Filter-Batching.
        """
        import time as _time
        self._seg_no = getattr(self, '_seg_no', 0) + 1
        if (os.environ.get('PYCBC_RATIO_CPROF')
                and self._seg_no == int(os.environ.get('PYCBC_RATIO_CPROF'))
                and not getattr(self, '_cp_done', False)):
            import cProfile, pstats, sys
            self._cp_done = True
            pr = cProfile.Profile(); pr.enable()
            try:
                self._cp_done = True
                return self._execute_blocked_kernel(
                    data, filters_f, n_taps, valid_slice)
            finally:
                pr.disable()
                st = pstats.Stats(pr, stream=sys.stderr).sort_stats('tottime')
                print("[cprof] one segment of _execute_blocked_kernel",
                      file=sys.stderr)
                st.print_stats(12)
        _t_enter = _time.perf_counter()
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

        if os.environ.get('PYCBC_RATIO_DTYPE'):
            import sys
            print("[dtype] ref_snr %s  filters_f %s  n=%d"
                  % (data.dtype, filters_f.dtype, len(data)), file=sys.stderr)
        nblocks = 0
        block_f_cache = {}
        if self._engine_mode in ('hier', 'check') and self._ap_ref is None:
            # Once, not per segment.  The reference is the SNR distribution of
            # the group's coarse template under this PSD -- a property of the
            # bank, not of the stretch of data being filtered -- and measuring
            # it runs noise realisations inside matchedfilter, which cost 21% of the
            # kernel when repeated every segment.
            self._ap_ref = self._set_engine_reference(data)
        # Do NOT reset what is loaded: the filter bank does not change between
        # segments, so re-ingesting all of it every segment is pure repetition.
        # Keyed on the array's identity as well as the batch offset, so a
        # genuinely different bank still reloads.
        for f_start in range(0, n_filters, self.batch_size):

            f_end = min(f_start + self.batch_size, n_filters)
            actual_batch_size = f_end - f_start

            if self._engine_mode:
                ap_plan = self._get_engine_plan(actual_batch_size)
                tag = (f_start, id(filters_f))
                if self._ap_loaded != tag:
                    # matchedfilter conjugates the template itself, and only the bins
                    # the analytic kernel writes may contribute -- the rest must
                    # be zero or matchedfilter would fold in the half pycbc leaves out.
                    tmpl = np.conj(filters_f[f_start:f_end]).copy()
                    tmpl[:, N_FFT // 2 + 1:] = 0
                    ap_plan.set_templates(tmpl)
                    self._ap_loaded = tag
                    self._ap_filters = filters_f
                ifft_plan = current_mult_view = current_corr_view = None
                if self._engine_mode == 'check':
                    (self._chk_ifft, self._chk_mult,
                     self._chk_corr) = self._get_ifft_plan(actual_batch_size, N_FFT)
            else:
                ifft_plan, current_mult_view, current_corr_view = self._get_ifft_plan(
                    actual_batch_size, N_FFT
                )

            # Overlap-save: valid output samples per block.
            n_taps_max = n_taps[f_start:f_end].max()
            i = np.searchsorted(nsizes, n_taps_max)
            n_taps_max = nsizes[i]

            N_VALID = N_FFT - n_taps_max + 1
            STEP = N_VALID
            bad_start = n_taps_max // 2

            first_block_idx = (v_start - bad_start) // STEP
            loop_start = first_block_idx * STEP

            # matchedfilter is a D x T engine: every block is a data segment and every
            # filter a template.  Driving it one block at a time wastes that
            # entirely -- 245 calls of D=1 where a handful of D=245 would do,
            # with the per-call cost paid 245 times and the template ingest
            # repeated.  So collect the blocks first, group them by window
            # (they are identical except at the segment edges), and hand each
            # group over in one call.
            blocks = []
            if self._engine_mode == 'hier' and self._ap_series:
                # One call per (filter batch, segment): matchedfilter walks the block
                # layout itself, doing each block's forward transform inline.
                # The layout is still computed here -- matchedfilter only executes it.
                # Vectorised: the layout is pure arithmetic on the block index,
                # so there is no reason to walk it in Python.
                ts = np.arange(loop_start, n_samples, STEP, dtype=np.int64)
                bvt0 = ts + bad_start
                keep = (bvt0 < v_stop) & (bvt0 + N_VALID > v_start)
                if keep.any():
                    last = np.flatnonzero(bvt0 < v_stop)
                    keep &= np.arange(ts.size) <= last[-1]
                ts = ts[keep]
                rs = np.maximum(v_start, bvt0[keep])
                re_ = np.minimum(v_stop, bvt0[keep] + N_VALID)
                good = re_ > rs
                bstarts, bws, bwe = ts[good], (rs - ts)[good], (re_ - ts)[good]
                if bstarts.size:
                    nblocks += bstarts.size
                    ap_plan = self._get_engine_plan(actual_batch_size, 1)
                    tag = (f_start, id(filters_f))
                    if self._ap_loaded != tag:
                        tmpl = np.conj(filters_f[f_start:f_end]).copy()
                        tmpl[:, N_FFT // 2 + 1:] = 0
                        ap_plan.set_templates(tmpl)
                        self._ap_loaded = tag
                        self._ap_filters = filters_f   # keep id() from being reused
                    _b = _time.perf_counter()
                    aidx, aval, _ = ap_plan.run_series(
                        data, bstarts, bws, bwe, binsize=N_FFT,
                        threshold=self.snr_threshold,
                        templates=(0, actual_batch_size), raw=True)
                    self._ph_ap = getattr(self, '_ph_ap', 0.0) + \
                        _time.perf_counter() - _b
                    ii = aidx[:, :, 0]
                    bi, fi = np.nonzero(ii >= 0)
                    if bi.size:
                        tsa = bstarts.astype(np.int64)
                        all_f_idxs.extend((f_start + fi).tolist())
                        all_t_idxs.extend((tsa[bi] + ii[bi, fi]).tolist())
                        all_snrs.extend(aval[:, :, 0][bi, fi].tolist())
                        all_tstarts.extend(tsa[bi].tolist())
                continue

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
                nblocks += 1
                buf_slice_start = roi_start - t_start
                t_end = min(t_start + N_FFT, n_samples)
                if t_start not in block_f_cache:
                    _a = _time.perf_counter()
                    block_f_cache[t_start] = self._forward_block_fft(
                        data[t_start:t_end])
                    self._ph_fft = getattr(self, '_ph_fft', 0.0) + \
                        _time.perf_counter() - _a
                blocks.append((t_start, buf_slice_start, roi_len))

            if self._engine_mode in ('flat', 'hier'):
                groups = {}
                for t_start, bss, rl in blocks:
                    groups.setdefault((bss, rl), []).append(t_start)
                for (bss, rl), starts in groups.items():
                    for chunk0 in range(0, len(starts), self._ap_ndata):
                        chunk = starts[chunk0:chunk0 + self._ap_ndata]
                        ap_plan = self._get_engine_plan(actual_batch_size,
                                                        len(chunk))
                        if self._ap_loaded != (f_start, id(ap_plan)):
                            tmpl = np.conj(filters_f[f_start:f_end]).copy()
                            tmpl[:, N_FFT // 2 + 1:] = 0
                            ap_plan.set_templates(tmpl)
                            self._ap_loaded = (f_start, id(ap_plan))
                        _c = _time.perf_counter()
                        if len(chunk) == 1:
                            # No copy: the cached block spectrum is already
                            # contiguous complex64, so hand the view straight in.
                            ap_plan.set_data(block_f_cache[chunk[0]][None, :])
                        else:
                            stack = np.empty((len(chunk), N_FFT), np.complex64)
                            for i, ts in enumerate(chunk):
                                stack[i] = block_f_cache[ts]
                            ap_plan.set_data(stack)
                        _b = _time.perf_counter()
                        self._ph_sd = getattr(self, '_ph_sd', 0.0) + \
                            _time.perf_counter() - _c
                        peaks = ap_plan.run(binsize=rl,
                                            threshold=self.snr_threshold,
                                            window=(bss, bss + rl))
                        self._ph_ap = getattr(self, '_ph_ap', 0.0) + \
                            _time.perf_counter() - _b
                        idx = peaks['index'][:, :, 0]
                        val = peaks['value'][:, :, 0]
                        di, fi = np.nonzero(idx >= 0)
                        if di.size:
                            ts_arr = np.asarray(chunk, dtype=np.int64)
                            all_f_idxs.extend((f_start + fi).tolist())
                            all_t_idxs.extend((ts_arr[di] + idx[di, fi]).tolist())
                            all_snrs.extend(val[di, fi].tolist())
                            all_tstarts.extend(ts_arr[di].tolist())
                continue

            for t_start, buf_slice_start, roi_len in blocks:
                roi_start = t_start + buf_slice_start
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

        if self._engine_mode == 'check' and self._chk_tot:
            import sys
            snrs = np.array(self._chk_missed_snr) if self._chk_missed_snr else np.zeros(0)
            print("[matchedfilter-check] gate dismissed %d of %d filter-block hits (%.1f%%)"
                  % (self._chk_miss, self._chk_tot,
                     100.0 * self._chk_miss / self._chk_tot), file=sys.stderr)
            if len(snrs):
                print("[matchedfilter-check] dismissed |snr|: min %.3f med %.3f max %.3f"
                      % (snrs.min(), np.median(snrs), snrs.max()), file=sys.stderr)
            for full, coarse, fb, bd, rl in self._chk_detail[:4]:
                print("[matchedfilter-check]   full=%.3f coarse(scaled)=%.3f  f_band=%.4f "
                      "band=%d roi=%d" % (full, coarse, fb, bd, rl), file=sys.stderr)
            self._chk_detail = []
            self._chk_tot = self._chk_miss = 0
            self._chk_missed_snr = []
        self._chk_detail = []
        self._ebk_ret = None
        self._kernel_seconds = getattr(self, '_kernel_seconds', 0.0) + (
            _time.perf_counter() - _t_enter)
        self._kernel_calls = getattr(self, '_kernel_calls', 0) + 1
        self._kernel_blocks = getattr(self, '_kernel_blocks', 0) + nblocks
        if os.environ.get('PYCBC_RATIO_TIMING'):
            import sys
            print("[ratio-parts2] fft=%.3f set_data=%.3f run=%.3f  rest=%.3f s"
                  % (getattr(self,'_ph_fft',0), getattr(self,'_ph_sd',0),
                     getattr(self,'_ph_ap',0),
                     self._kernel_seconds - getattr(self,'_ph_fft',0)
                     - getattr(self,'_ph_sd',0) - getattr(self,'_ph_ap',0)),
                  file=sys.stderr)
            print("[ratio-parts] fft=%.3f matchedfilter=%.3f extract=%.3f other=%.3f s"
                  % (getattr(self,'_ph_fft',0), getattr(self,'_ph_ap',0),
                     getattr(self,'_ph_out',0),
                     self._kernel_seconds - getattr(self,'_ph_fft',0)
                     - getattr(self,'_ph_ap',0) - getattr(self,'_ph_out',0)),
                  file=sys.stderr)
            print("[ratio-timing] mode=%-5s segments=%d filter-blocks=%d "
                  "kernel=%.3f s (%.3f ms/filter-block)"
                  % (self._engine_mode or 'stock', self._kernel_calls,
                     self._kernel_blocks, self._kernel_seconds,
                     1e3 * self._kernel_seconds / max(self._kernel_blocks, 1)),
                  file=sys.stderr)
        return (np.array(all_f_idxs, dtype=np.int32),
                np.array(all_t_idxs, dtype=np.int64),
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
