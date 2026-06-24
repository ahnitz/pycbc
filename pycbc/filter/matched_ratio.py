import logging
import time
import numpy as np
import time
import mkl_fft
from pycbc.types import zeros, complex64, TimeSeries, FrequencySeries
from pycbc.filter.matchedfilter import matched_filter_core, sigmasq
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython, find_peaks_fused_lanczos_cython, block_max_threshold, fast_multiply_scale_cython, fast_multiply_indexed_cython, find_peaks_indexed_cython

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
        self.hnorms_low = None

        self.flen = int(self.multiband_frequency / self.delta_f)
        self.lowband_snrs = np.empty((self.batch_size, self.flen), dtype=np.complex64) 

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
        
        lowband_templates, lowband_templates2 = lowband_templates
        
        # Cache if we keep the same ref_template, otherwise remake
        if not hasattr(stilde, 'ref_snr') or stilde.ref_template_id != id(ref_template):
            h_norm = sigmasq(ref_template, psd=psd, 
                             low_frequency_cutoff=self.multiband_frequency,
                             high_frequency_cutoff=self.f_high)    
            snr, _, norm = matched_filter_core(
                ref_template, stilde, psd=psd,
                low_frequency_cutoff=self.multiband_frequency,
                high_frequency_cutoff=self.f_high,
                h_norm=h_norm
            )
            ref_snr = snr.numpy() * (norm * stilde.delta_t)
            stilde.ref_hnorm = h_norm    
            stilde.ref_snr = ref_snr
            stilde.ref_template_id = id(ref_template)
        
        ref_snr = stilde.ref_snr
        ref_hnorm = stilde.ref_hnorm
        
        # Calculate lowband template SNRs
        logging.info('processing lowband templates')
        ltem = lowband_templates[0]
        hnorm_low_ref = sigmasq(ltem, psd=psd,
              low_frequency_cutoff=ltem.f_lower,
              high_frequency_cutoff=self.multiband_frequency)
                 
        # This won't work if we call this function multiple times with different templates!      
        dec = 5        
        if self.hnorms_low is None:  
            print('recalculating hnorms')     
            hnorms_low = [sigmasq(ltem[::dec], psd=psd[::dec],
                              low_frequency_cutoff=ltem.f_lower,
                              high_frequency_cutoff=self.multiband_frequency) for ltem in lowband_templates]
            self.hnorms_low = np.array(hnorms_low)      
            
        hnorms_low = self.hnorms_low * (hnorm_low_ref / self.hnorms_low[0])

        #hnorms_low = [sigmasq(ltem[::dec], psd=psd[::dec],
        #                  low_frequency_cutoff=ltem.f_lower,
        #                  high_frequency_cutoff=self.multiband_frequency) for ltem in lowband_templates]
        #hnorms_low = np.array(hnorms_low)

        norm_low = 4.0 * stilde.delta_f / hnorms_low ** 0.5        
        slow, shigh = scale
        hnorms_high = ref_hnorm * shigh ** 2.0
        sigma_total = hnorms_low + hnorms_high
        rw_low = (hnorms_low / sigma_total) ** 0.5
        rw_high = (hnorms_high / sigma_total) ** 0.5

        dt_low = 1.0 / (self.flen * stilde.delta_f)
        sow = stilde[:self.flen] / psd[:self.flen]
        sow2 = sow.astype(np.complex64).numpy()
        kmax = min(len(lowband_templates[0]), self.flen)        

        t2 = time.time()   
        
        lowband_snrs = self.lowband_snrs[0:len(lowband_templates)]
        scale_factors = (norm_low * rw_low * self.flen).astype(np.float32)
        fast_multiply_scale_cython(sow2, lowband_templates2, scale_factors, lowband_snrs)
        
        t22 = time.time()
        self.fft_lib.ifft(lowband_snrs, out=lowband_snrs, axis=-1)
        t3 = time.time()        
        
        # Prenormalize the SNRs
        # Opt note: could move all this multiplication to the peak finding..
        filters_f  = (filters_f * rw_high[:, None]).astype(filters_f.dtype)
           
        logging.info('.....done')                
        t4 = time.time()

        safety_sigma = 4.0
        extra_margin = 0.05
        gate_threshold = np.maximum(0.0, self.snr_threshold * rw_low - safety_sigma * rw_high - extra_margin)

        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals, snr_low_vals, tstarts = self._execute_blocked_kernel(
            ref_snr, filters_f, n_taps, valid_slice, stilde,
            lowband_snrs=(lowband_snrs, dt_low, gate_threshold),
        )
        
        # Undo coherent addition weighting so it is back to SNR normalized
        snr_high_vals = snr_vals - snr_low_vals
        snr_high_vals = snr_high_vals / rw_high[local_idxs]
        snr_low_vals = snr_low_vals / rw_low[local_idxs]
        t5 = time.time()
         
        print("PRE TIMING", "RF", t5-t4, "NORM", t4 - t3, "LB", t3 - t2, t3-t22, t22-t2, "REF", t2-t1)
        # 4. Map indices
        if len(local_idxs) > 0:
            global_ids = indices[local_idxs]
            return global_ids, t_idxs, snr_vals, snr_low_vals, snr_high_vals, tstarts, hnorms_low[local_idxs]
        else:
            return [], [], [], [], [], [], []

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
        all_low_snrs = []
        all_tstarts = []
        indices = np.arange(n_filters)
        
        if not hasattr(stilde, 'block_f_cache'):
            stilde.block_f_cache = {}

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
        
        tspend = 0
        tspend1 = 0
        tspend2 = 0
        tspend3 = 0

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
         
        n_time_blocks = len(t_starts)
        time_chunks_with_work = np.any(valid, axis=0)
        valid_t_idxs = np.where(time_chunks_with_work)[0]

        #print(valid.sum() / (valid.shape[0] * valid.shape[1]), time_chunks_with_work.sum(), n_time_blocks)

        tx2 = time.time()    
        for t_idx in valid_t_idxs:
            t_start = t_starts[t_idx]

            # Compute or fetch FFT once per time block
            if t_start not in stilde.block_f_cache:
                t_end = min(t_start + N_FFT, n_samples)
                stilde.block_f_cache[t_start] = self.fft_lib.fft(data[t_start:t_end])
            block_f_view = stilde.block_f_cache[t_start]

            active_f_idxs = indices[valid[:, t_idx]]
            n_active = len(active_f_idxs)
            
            if n_active == 0:
                continue

            roi_start = roi_starts[t_idx]
            roi_stop = roi_stops[t_idx]
            roi_len = roi_stop - roi_start
            buf_slice_start = roi_start - t_start

            # --- BULK EXECUTION PHASE ---
            t1 = time.time()      
            
            # 1. Dynamically size the views for this specific time block
            current_mult_view = freq_mult_view[:n_active]
            current_corr_view = corr_out_view[:n_active]

            # 2. BULK C-MULTIPLY (1 Python-to-C hop for all ~10 filters)
            fast_multiply_indexed_cython(
                block_f_view, 
                filters_f, 
                active_f_idxs, 
                current_mult_view
            )
            t11 = time.time()

            # 3. BULK MKL IFFT (1 Engine Spin-up for all ~10 rows)
            self.fft_lib.ifft(
                current_mult_view,
                axis=-1, 
                out=current_corr_view
            )     
            t12 = time.time()        

            # 4. BULK C-PEAK FINDER (1 Python-to-C hop for all ~10 rows)
            f_list, t_list, s_list, s_low_list = find_peaks_indexed_cython(
                current_corr_view,
                lowband_snrs,
                active_f_idxs,
                roi_start, roi_len,
                self.threshold_sq,
                dt_high, dt_low,
                t0_high, t0_low,
                input_offset=buf_slice_start
            )
            
            # Append aggregated results for the entire batch
            if f_list:
                all_f_idxs.extend(f_list)
                all_t_idxs.extend(t_list)
                all_snrs.extend(s_list)
                all_low_snrs.extend(s_low_list)
                all_tstarts.extend([t_start] * len(s_list)) 
            
            t13 = time.time()
            
            # Aggregate the timers
            tspend += (t13 - t1)
            tspend1 += (t11 - t1)
            tspend2 += (t12 - t11)
            tspend3 += (t13 - t12)
                
        tf = time.time()
        print("SKIP RATIO", "INNER", tspend, tspend1, tspend2, tspend3, "top", tx2-t0, "total", tf-t0)        
             
        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_low_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
