import logging
import time
import numpy as np
import time
import mkl_fft
from pycbc.types import zeros, complex64, TimeSeries, FrequencySeries
from pycbc.filter.matchedfilter import matched_filter_core, sigmasq
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython, find_peaks_fused_lanczos_cython, block_max_threshold

class RatioMatchedFilterControl(object):
    """
    High-performance engine for hierarchical "Ratio/FIR" matched filtering.
    Uses mkl_fft for ALL FFT operations to maximize throughput and consistency.
    """

    def __init__(self, snr_threshold, delta_f, 
                 high_frequency_cutoff=None,
                 fir_fft_length=4096,
                 batch_size=64, 
                 multiband_frequency=None):
        self.delta_f = delta_f
        self.snr_threshold = snr_threshold
        self.f_high = high_frequency_cutoff
        
        self.threshold_sq = float(snr_threshold**2)
        
        self.fir_fft_len = fir_fft_length
        self.batch_size = batch_size
        
        # 2. Intermediate Buffers (Batch x Block)
        total_batch_size = batch_size * fir_fft_length
        self.temp_freq_mult = zeros(total_batch_size, dtype=complex64).data.reshape(self.batch_size, fir_fft_length)
        self.corr_output_buffer = zeros(total_batch_size, dtype=complex64).data.reshape(self.batch_size, fir_fft_length)

        # 3. Filter Preparation Buffers
        self.filters_padded = zeros(total_batch_size, dtype=complex64)
        self.filters_f_buffer = zeros(total_batch_size, dtype=complex64)
        
        # Direct MKL handle
        self.fft_lib = mkl_fft
        
        self.multiband_frequency = multiband_frequency

    def prepare_filters(self, fir_taps, tap_counts):
        """
        Prepare frequency-domain filters for a batch of taps.
        """
        n_filters, n_taps = fir_taps.shape
        if n_taps >= self.fir_fft_len:
             raise ValueError("FIR Taps (%d) exceed FFT block length (%d)" % 
                              (n_taps, self.fir_fft_len))
        
        # Calculate max tap count for validity logic
        n_taps_max = int(np.max(tap_counts))
        
        filters_f = self._fft_all_filters(fir_taps, tap_counts)
        return filters_f, n_taps_max

    def process_segment(self, stilde, psd, ref_template, filters_f, n_taps, indices, 
                        valid_slice=None, lowband_templates=None, scale=None):
        """
        Process a single data segment.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)
        
        t1 = time.time()
        h_norm = sigmasq(ref_template, psd=psd, 
                         low_frequency_cutoff=self.multiband_frequency,
                         high_frequency_cutoff=self.f_high)    
        snr, _, norm = matched_filter_core(
            ref_template, stilde, psd=psd,
            low_frequency_cutoff=self.multiband_frequency,
            high_frequency_cutoff=self.f_high,
            h_norm=h_norm
        )
        t2 = time.time()    
                
        # Calculate lowband template SNRs
        logging.info('processing lowband templates')
        ltem = lowband_templates[0]
        hnorm_low_ref = sigmasq(ltem, psd=psd,
              low_frequency_cutoff=ltem.f_lower,
              high_frequency_cutoff=self.multiband_frequency)
                 
        # This won't work if we call this function multiple times with different templates!      
        dec = 5        
        if not hasattr(self, 'hnorms_low'):        
            hnorms_low = [sigmasq(ltem[::dec], psd=psd[::dec],
                              low_frequency_cutoff=ltem.f_lower,
                              high_frequency_cutoff=self.multiband_frequency) for ltem in lowband_templates]
            self.hnorms_low = np.array(hnorms_low)      
            
        hnorms_low = self.hnorms_low * (hnorm_low_ref / self.hnorms_low[0])
        norm_low = 4.0 * stilde.delta_f / hnorms_low ** 0.5        
        slow, shigh = scale
        hnorms_high = h_norm * shigh ** 2.0
        sigma_total = hnorms_low + hnorms_high
        rw_low = (hnorms_low / sigma_total) ** 0.5
        rw_high = (hnorms_high / sigma_total) ** 0.5

        flen = int(self.multiband_frequency / stilde.delta_f)
        dt_low = 1.0 / (flen * stilde.delta_f)
        sow = stilde[:flen] / psd[:flen]
        sow2 = sow.astype(np.complex64).numpy()
        
        kmax = min(len(lowband_templates[0]), flen)
        lowband_snrs = np.resize(sow2, len(lowband_templates) * flen).reshape((len(lowband_templates), flen)) 
        for j, ltemplate in enumerate(lowband_templates):
            lowband_snrs[j][:kmax] *= ltemplate[:kmax].conj() * (norm_low[j] * rw_low[j] * flen)
             
        t22 = time.time()
        self.fft_lib.ifft(lowband_snrs, out=lowband_snrs, axis=-1)
        t3 = time.time()        
        
        # Prenormalize the SNRs
        # Opt note: could move all this multiplication to the peak finding..
        self.ref_snr = snr.numpy() * (norm * stilde.delta_t)
        filters_f  = (filters_f * rw_high[:, None]).astype(filters_f.dtype)
           
        logging.info('.....done')                
        t4 = time.time()

        safety_sigma = 4.0
        extra_margin = 0.15
        gate_threshold = max(0.0, self.snr_threshold - safety_sigma * np.mean(rw_high) - extra_margin)

        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
            self.ref_snr, filters_f, n_taps, valid_slice, stilde,
            lowband_snrs=(lowband_snrs, dt_low, gate_threshold),
        )
        t5 = time.time()
         
        print("PRE TIMING", "RF", t5-t4, "NORM", t4 - t3, "LB", t3 - t2, t3-t22, t22-t2, "REF", t2-t1)
        
        # 4. Map indices
        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, tstarts, h_norm
        else:
            return [], [], [], tstarts, h_norm

    def _fft_all_filters(self, taps, counts):
        """Helper to FFT all filters using mkl_fft."""
        n_filters, n_taps_alloc = taps.shape
        filters_f = np.zeros((n_filters, self.fir_fft_len), dtype=np.complex64)
        
        padded_view = self.filters_padded.data
        padded_reshaped = padded_view.reshape(self.batch_size, self.fir_fft_len)

        for start in range(0, n_filters, self.batch_size):
            end = min(start + self.batch_size, n_filters)
            batch_len = end - start
            
            # 1. Zero out and Fill
            padded_reshaped[:batch_len, :] = 0
            
            tmp_taps = taps[start:end]
            padded_reshaped[:batch_len, :n_taps_alloc] = tmp_taps
            
            # 2. Variable Roll Logic
            current_counts = counts[start:end]
            roll_offsets = -(current_counts // 2)
            
            cols = np.arange(self.fir_fft_len)
            rows = np.arange(batch_len)[:, None]
            shifted_cols = (cols[None, :] - roll_offsets[:, None]) % self.fir_fft_len
            
            current_data = padded_reshaped[:batch_len].copy()
            padded_reshaped[:batch_len] = current_data[rows, shifted_cols]

            # 3. Execute FFT (Direct MKL call)
            fft_out = self.fft_lib.fft(padded_reshaped[:batch_len], axis=-1)
            
            # 4. Conjugate & Store
            kmax = self.fir_fft_len // 2
            filters_f[start:end][:kmax] = np.conj(fft_out[:kmax])
            
        return filters_f

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice, stilde, lowband_snrs=None):
        """
        Inner loop: Time-Blocking + Filter-Batching using mkl_fft.
        """
        t0 = time.time()

        N_FFT = self.fir_fft_len        
        n_taps_max = n_taps.max()
        n_samples = len(data)
        n_filters = len(filters_f)       
        N_VALID = N_FFT - n_taps_max + 1
        bad_start = n_taps_max // 2      
        
        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        all_tstarts = []

        freq_mult_view = self.temp_freq_mult
        corr_out_view = self.corr_output_buffer

        if lowband_snrs is not None:
            # Let's just use exactly matched series for now...
            lowband_snrs, dt_low, gate_threshold = lowband_snrs
            t0_low = t0_high = float(stilde.start_time)
            dt_high = stilde.delta_t

        if valid_slice:
            v_start = valid_slice.start
            v_stop = valid_slice.stop
        else:
            v_start = 0
            v_stop = n_samples

        block_f_cache = {}
        
        tspend = 0

        # Determine Loop Bounds
        first_block_idx = (v_start - bad_start) // N_VALID
        loop_start = first_block_idx * N_VALID

        # Set up the time chunk book keeping info
        t_starts = np.arange(loop_start, n_samples, N_VALID)
        valid = np.ones(len(t_starts), dtype=bool)
        block_valid_t0s = t_starts + bad_start

        # Remove parts not within the total valid data boundaries
        valid[block_valid_t0s >= v_stop] = False
        valid[block_valid_t0s + N_VALID <= v_start] = False
        
        roi_starts = np.maximum(v_start, block_valid_t0s)
        roi_stops = np.minimum(v_stop, block_valid_t0s + N_VALID)
        
        roi_len = roi_stops - roi_starts
        valid[roi_len <= 0] = False
       
        j_starts = np.round(roi_starts * dt_high / dt_low).astype(np.int32) - 1
        j_starts = np.maximum(0, j_starts)
        j_ends = np.round(roi_stops * dt_high / dt_low).astype(np.int32) + 1
        j_ends = np.minimum(lowband_snrs.shape[1], j_ends)

        valid = np.resize(valid, len(filters_f) * len(t_starts)).reshape(len(filters_f), len(t_starts))

        
        block_max_threshold(lowband_snrs, j_starts, j_ends, valid, gate_threshold)
        current_mult_view = freq_mult_view[:1]
        current_corr_view = corr_out_view[:1]

        print(valid.sum() / (valid.shape[0] * valid.shape[1]))

        tx2 = time.time()             
        for f_start in range(n_filters): 
            f_end = f_start + 1
            low_snr_block = lowband_snrs[f_start:f_end]
            low_snr_fbatch = lowband_snrs[f_start:f_end]
            filter_batch_f = filters_f[f_start:f_end]
            valid2 = valid[f_start]
            
            t_starts2 = t_starts[valid2]
            roi_starts2 = roi_starts[valid2]
            roi_stops2 = roi_stops[valid2]
            for t_start, roi_start, roi_stop in zip(t_starts2, roi_starts2, roi_stops2):
                t1 = time.time()      
                roi_len = roi_stop - roi_start

                buf_slice_start = roi_start - t_start
                t_end = min(t_start + N_FFT, n_samples)

                if t_start not in block_f_cache:
                    block_f_view = self.fft_lib.fft(data[t_start:t_end])
                    block_f_cache[t_start] = block_f_view
                block_f_view = block_f_cache[t_start]

                fast_multiply_analytic_cython(
                    block_f_view, filter_batch_f, current_mult_view
                )

                self.fft_lib.ifft(
                    current_mult_view, 
                    axis=-1, 
                    out=current_corr_view
                )             
                f_list, t_list, s_list = find_peaks_fused_lanczos_cython(
                                        current_corr_view,
                                        low_snr_block,
                                        roi_start, roi_len,
                                        self.threshold_sq,
                                        f_start,
                                        dt_high, dt_low,
                                        t0_high, t0_low,
                                        input_offset=buf_slice_start
                                    )

                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)
                    all_tstarts.extend([t_start] * len(s_list)) 
                t2 = time.time()
                tspend += t2 - t1
        tf = time.time()
        print("SKIP RATIO", "INNER", tspend, "top", tx2-t0, "total", tf-t0)        
             
        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
