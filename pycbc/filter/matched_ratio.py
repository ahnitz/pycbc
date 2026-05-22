import logging
import numpy as np
import mkl_fft
from pycbc.types import zeros, complex64, TimeSeries
from pycbc.filter.matchedfilter import matched_filter_core, sigmasq
from pycbc.filter.matchedfilter_cpu import fast_multiply_analytic_cython, find_peaks_in_block_cython, find_peaks_fused_lanczos_cython

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
                        valid_slice=None, lowband_templates=None):
        """
        Process a single data segment.
        """
        if valid_slice is None:
            valid_slice = getattr(stilde, 'analyze', None)

        import time
        
        t1 =time.time()
        # 2. Calculate Reference SNR
        flow = self.multiband_frequency if self.multiband_frequency is not None else ref_template.f_lower        
        h_norm = sigmasq(ref_template, psd=psd, 
                         low_frequency_cutoff=flow,
                         high_frequency_cutoff=self.f_high)    
        snr, _, norm = matched_filter_core(
            ref_template, stilde, psd=psd,
            low_frequency_cutoff=flow,
            high_frequency_cutoff=self.f_high,
            h_norm=h_norm
        )
        
        # Calculate lowband template SNRs
        logging.info('processing lowband templates')

        # Fix up this normalization step to do mchirp rescaling?
        h_norm2 = sigmasq(lowband_templates[0], psd=psd,
                         low_frequency_cutoff=lowband_templates[0].f_lower,
                         high_frequency_cutoff=self.multiband_frequency)
        norm2 =  4.0 * stilde.delta_f / h_norm2 ** 0.5

        # Reweighting factors so you can just add the high / low
        # snrs directly
        h_norm = 1.0 / norm**2.0
        h_norm2 = 1.0 / norm2**2.0
        sigma_total = (h_norm + h_norm2)**0.5
        rw_low = h_norm2 ** 0.5 / sigma_total
        rw_high = h_norm ** 0.5 / sigma_total
        
        flen = int(self.multiband_frequency / stilde.delta_f)
        sow = stilde[:flen] / psd[:flen]
        sow2 = sow.numpy() * rw_low * norm2 * flen

        t2 = time.time()
        lowband_snrs = np.zeros((len(lowband_templates), flen), dtype=np.complex64)
        
        for j, ltemplate in enumerate(lowband_templates):
            lowband_snrs[j] = ltemplate[:flen].conj()      
            lowband_snrs[j] *= sow2
       
        t22 = time.time()
        self.fft_lib.ifft(lowband_snrs, out=lowband_snrs, axis=-1)
        t3 = time.time()        
        
        # Prenormalize the SNRs
        # Opt note: could move all this multiplication to the peak finding..
        self.ref_snr = snr.numpy() * (norm * stilde.delta_t * rw_high)
           
        logging.info('.....done')                
        t4 = time.time()

        gate_threshold_sq = max(0.0, self.snr_threshold - 4.0 * rw_high) ** 2.0

        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
            self.ref_snr, filters_f, n_taps, valid_slice, stilde,
            lowband_snrs=(lowband_snrs, 1.0 / self.multiband_frequency, gate_threshold_sq),
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
            filters_f[start:end] = np.conj(fft_out)
            
        return filters_f

    def _execute_blocked_kernel(self, data, filters_f, n_taps, valid_slice, stilde, lowband_snrs=None):
        """
        Inner loop: Time-Blocking + Filter-Batching using mkl_fft.
        """
        tap_groups = 3
        nsizes = np.quantile(n_taps, np.linspace(0, 1, tap_groups+1)[1:]).astype(int)
        n_samples = len(data)
        n_filters = len(filters_f)
        
        N_FFT = self.fir_fft_len
        
        all_f_idxs = []
        all_t_idxs = []
        all_snrs = []
        all_tstarts = []

        freq_mult_view = self.temp_freq_mult
        corr_out_view = self.corr_output_buffer

        if lowband_snrs is not None:
            # Let's just use exactly matched series for now...
            lowband_snrs, dt_low, gate_threshold_sq = lowband_snrs
            t0_low = t0_high = float(stilde.start_time)
            dt_high = stilde.delta_t

        if valid_slice:
            v_start = valid_slice.start
            v_stop = valid_slice.stop
        else:
            v_start = 0
            v_stop = n_samples

        block_f_cache = {}
        # --- OUTER LOOP: Time Blocks ---
        for f_start in range(0, n_filters, self.batch_size):  
        
        
            f_end = min(f_start + self.batch_size, n_filters)
            actual_batch_size = f_end - f_start
 
            low_snr_block = lowband_snrs[f_start:f_end]
            current_mult_view = freq_mult_view[:actual_batch_size]
            current_corr_view = corr_out_view[:actual_batch_size]
            
            # Valid output samples per block (Overlap-Save)
            n_taps_max = n_taps[f_start:f_end].max()
            i = np.searchsorted(nsizes, n_taps_max)
            n_taps_max = nsizes[i]
            
            N_VALID = N_FFT - n_taps_max + 1
            STEP = N_VALID
            bad_start = n_taps_max // 2

            # Determine Loop Bounds
            first_block_idx = (v_start - bad_start) // STEP
            loop_start = first_block_idx * STEP

            for t_start in range(loop_start, n_samples, STEP):

                # ... inside _execute_blocked_kernel, inside the time block loop:

                # 1. Calculate the lowband index range for this time block (same logic as Cython)
                # We need the indices corresponding to [t_start, t_start + N_VALID]
                # t0_high/dt_high defines the global highband time axis.
                # dt_low is 1.0 / self.multiband_frequency
                
                c0 = (float(t_start) * dt_high) / dt_low
                c1 = (float(t_start + N_VALID) * dt_high) / dt_low
                
                j_start = int(np.floor(c0))
                j_end = int(np.ceil(c1))
                
                # Ensure bounds are safe for slicing
                j_start = max(0, j_start)
                j_end = min(lowband_snrs.shape[1], j_end)

                # 2. Extract the lowband block for this batch
                # lowband_snrs is shape (n_filters, flen)
                low_snr_batch = lowband_snrs[f_start:f_end, j_start:j_end]

                # 3. Probabilistic Gate
                # Using 3.4 as the amplitude gate, 3.4^2 = 11.56
                # If the max squared amplitude in this slice is too low, skip this batch
                if np.max(np.abs(low_snr_batch)**2) < gate_threshold_sq:
                    continue
                    
                # ... proceed to FFT, multiply, and IFFT ...

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
                    block_in_view = np.zeros(self.fir_fft_len, dtype=complex64)
                    block_in_view[0:t_end-t_start] = data[t_start:t_end]
                    
                    block_f_view = self.fft_lib.fft(block_in_view)
                    block_f_cache[t_start] = block_f_view
                
                block_f_view = block_f_cache[t_start]
                filter_batch_f = filters_f[f_start:f_end]
                low_snr_batch = lowband_snrs[f_start:f_end]
                
                fast_multiply_analytic_cython(
                    block_f_view, filter_batch_f, current_mult_view
                )

                self.fft_lib.ifft(
                    current_mult_view, 
                    axis=-1, 
                    out=current_corr_view
                )

                if lowband_snrs is not None:
                    #t1 = time.time()
                    f_list, t_list, s_list = find_peaks_fused_lanczos_cython(
                                            current_corr_view,
                                            low_snr_block,
                                            roi_start,
                                            roi_len,
                                            self.threshold_sq,
                                            f_start,
                                            dt_high,
                                            dt_low,
                                            t0_high,
                                            t0_low,
                                            input_offset=buf_slice_start
                                        )
                else:      
                    f_list2, t_list2, s_list2 = find_peaks_in_block_cython(
                        current_corr_view, 
                        roi_start,          
                        roi_len,            
                        self.threshold_sq, 
                        f_start,
                        input_offset=buf_slice_start
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
