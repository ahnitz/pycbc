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
        h_norm = sigmasq(ref_template, psd=psd, low_frequency_cutoff=flow, high_frequency_cutoff=self.f_high)
        
        snr, _, norm = matched_filter_core(
            ref_template, stilde, psd=psd,
            low_frequency_cutoff=flow,
            high_frequency_cutoff=self.f_high,
            h_norm=h_norm
        )
        
        self.ref_snr = snr.numpy() * (norm * stilde.delta_t)
        t2 = time.time()

        # Calculate lowband template SNRs
        
        logging.info('processing lowband templates')
        #lowband_snrs = []

        # Fix up this normalization step to do mchirp rescaling?
        h_norm2 = sigmasq(lowband_templates[0], psd=psd,
                         low_frequency_cutoff=lowband_templates[0].f_lower,
                         high_frequency_cutoff=self.multiband_frequency)
        norm2 =  4.0 * stilde.delta_f / h_norm2 ** 0.5
        if lowband_templates is not None:
            flen = len(lowband_templates[0])
            sow = stilde[:flen] / psd[:flen]
            
            #for ltemplate in lowband_templates:
            #    snr, _, norm2 = matched_filter_core(
            #        ltemplate, sow,
            #        low_frequency_cutoff=ltemplate.f_lower,
            #        high_frequency_cutoff=self.multiband_frequency,
            #        h_norm=h_norm2,
            #    )
            #    lowband_snrs += [snr.copy()]        

            flen2 = (flen - 1) * 2
            filter_batch_f = np.zeros((len(lowband_templates), flen2), dtype=np.complex64)
            for j, ltemplate in enumerate(lowband_templates):
                ltem = ltemplate.numpy().conj()
                kmax = int(self.multiband_frequency / ltemplate.delta_f)
                ltem[kmax:] = 0
                filter_batch_f[j, :flen] = ltem      
                
            sow2 = sow.numpy().copy()
            sow2.resize(flen2)
            
            tlen = len(lowband_templates) * flen2
            current_mult_view = zeros(tlen, dtype=np.complex64).data.reshape((len(lowband_templates), flen2))
            current_corr_view =  zeros(tlen, dtype=np.complex64).data.reshape((len(lowband_templates), flen2))

            t5 = time.time()
            
            if False:
                fast_multiply_analytic_cython(
                    sow2, filter_batch_f, current_mult_view
                )
                
                t55 = time.time()
                self.fft_lib.ifft(
                    current_mult_view, 
                    axis=-1, 
                    out=current_corr_view
                )
                lowband_snrs = [current_corr_view[i,:]
                                for i in range(len(lowband_templates))]
                                             
            else:      
                lowband_snrs = []
                for j in range(len(lowband_templates)):
                    x = sow2 * filter_batch_f[j, :]
                    lowband_snrs += [self.fft_lib.ifft(x)]
                        
            t6 = time.time()
            #breakpoint()
            
            print(t6-t5)

        lowband_snrs = [TimeSeries(x,
         delta_t=sow.delta_t, 
         epoch=stilde.start_time) for x in lowband_snrs]     
        
        # Reweight 
        h_norm = 1.0 / norm**2.0
        h_norm2 = 1.0 / norm2**2.0
        sigma_total = (h_norm + h_norm2)**0.5
        
        rw_low = h_norm2 ** 0.5 / sigma_total
        rw_high = h_norm ** 0.5 / sigma_total
        print(rw_low, rw_high)
        
        self.ref_snr *= rw_high
        lowband_snrs = [x * (rw_low * norm2 * flen2) for x in lowband_snrs]
        t3 = time.time()
       
        logging.info('.....done')                
            
        # 3. Execute Blocked Kernel
        local_idxs, t_idxs, snr_vals, tstarts = self._execute_blocked_kernel(
            self.ref_snr, filters_f, n_taps, valid_slice, stilde, lowband_snrs=lowband_snrs,
        )
        t4 = time.time()
        print("MF ref timing", t2 - t1, t3 - t2, t4 - t3)
        
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

            if lowband_snrs is not None:
                low_snr_batch = lowband_snrs[f_start:f_end]
                # Convert the batch into a contiguous 2D C-aligned array
                low_snr_block = np.ascontiguousarray(
                    np.stack([obj.numpy() if hasattr(obj, 'numpy') else obj.data for obj in low_snr_batch]),
                    dtype=np.complex64
                )
                # Extract grid references
                t0_low = float(lowband_snrs[0].start_time)
                dt_low = lowband_snrs[0].delta_t
                dt_high = stilde.delta_t
                t0_high = float(stilde.start_time)

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

                # --- NEW: Coherent Lowband SNR Addition with Bandlimited Cubic Approximation ---
# --- NEW: Coherent Lowband SNR Addition via Full FFT/IFFT Exact Interpolation ---

# --- NEW: Coherent Lowband SNR Addition via Fast Lanczos-3 Windowed Sinc Interpolation ---
# --- NEW: Coherent Lowband SNR Addition via Fast Lanczos-3 Windowed Sinc Interpolation ---
# --- NEW: Coherent Lowband SNR Addition via Arbitrary-Order Lanczos Interpolation ---
                if False and lowband_snrs is not None:
                    # =========================================================================
                    # CONTROLLABLE ACCURACY PARAMETER
                    # lanczos_a defines the kernel radius. Total stencil size = 2 * lanczos_a.
                    # Higher values (e.g., 6, 8, 12) increase accuracy toward exact FFT upsampling.
                    # =========================================================================
                    lanczos_a = 32 
                    
                    # 1. Get highband metadata
                    dt_high = stilde.delta_t
                    t0_high = float(stilde.start_time)
                    roi_high_samples = np.arange(roi_start, roi_stop)
                    
                    # 2. Get lowband metadata (using the first item as a grid reference)
                    t0_low = float(lowband_snrs[0].start_time)
                    dt_low = lowband_snrs[0].delta_t
                    
                    # Calculate continuous fractional lowband indices directly
                    exact_low_indices = ((t0_high - t0_low) + roi_high_samples * dt_high) / dt_low
                    idx0 = np.floor(exact_low_indices).astype(np.int64)
                    w = exact_low_indices - idx0
                    
                    # Construct a dynamic stencil coordinate matrix
                    stencil_offsets = np.arange(-lanczos_a + 1, lanczos_a + 1)
                    weights_list = []
                    
                    # Vectorized calculation of windowed sinc weights for the chosen radius
                    for k in stencil_offsets:
                        t = w - k
                        wk = np.sinc(t) * np.sinc(t / lanczos_a)
                        weights_list.append(wk)
                        
                    # Convert to a 2D array of shape: (2 * lanczos_a, roi_len)
                    weights = np.array(weights_list)
                    
                    # Normalize weights along the stencil axis to guarantee absolute gain stability
                    w_sum = np.sum(weights, axis=0)
                    weights /= w_sum
                    
                    # 3. Slice the buffer view corresponding to the ROI
                    buf_end = buf_slice_start + roi_len
                    
                    # 4. In-place interpolated addition for the filter batch
                    for b in range(actual_batch_size):
                        # Extract the raw numpy array to bypass PyCBC container overhead
                        low_snr_obj = low_snr_batch[b]
                        low_snr_arr = low_snr_obj.numpy() if hasattr(low_snr_obj, 'numpy') else low_snr_obj.data
                        
                        # Accumulate the stencil contributions across the entire ROI vector
                        interp_snr = np.zeros(roi_len, dtype=np.complex64)
                        for i, k in enumerate(stencil_offsets):
                            interp_snr += weights[i] * low_snr_arr[idx0 + k]
                            
                        # Add coherently into the pre-allocated highband buffer view
                        current_corr_view[b, buf_slice_start:buf_end] += interp_snr

                if False and lowband_snrs is not None:
                    # 1. Get highband metadata
                    dt_high = stilde.delta_t
                    t0_high = float(stilde.start_time)
                    
                    # 2. Get lowband metadata (using the first item as a grid reference)
                    t0_low = float(lowband_snrs[0].start_time)
                    dt_low = lowband_snrs[0].delta_t
                    
                    # Calculate the exact static sample index offset between the two grids
                    global_time_offset_idx = int(np.round((t0_high - t0_low) / dt_high))
                    
                    # Map the current highband ROI boundaries directly to the upsampled grid
                    slice_start = global_time_offset_idx + roi_start
                    slice_end = global_time_offset_idx + roi_stop
                    
                    # Slice the buffer view corresponding to the ROI
                    buf_end = buf_slice_start + roi_len
                    
                    # 3. Perform the exact Fourier upsampling per filter in the batch
                    for b in range(actual_batch_size):
                        # Extract the raw numpy array to bypass PyCBC container overhead
                        low_snr_obj = low_snr_batch[b]
                        low_snr_arr = low_snr_obj.numpy() if hasattr(low_snr_obj, 'numpy') else low_snr_obj.data
                        
                        N_low = len(low_snr_arr)
                        N_high = int(np.round(N_low * dt_low / dt_high))
                        
                        # Forward FFT to frequency domain
                        Y = self.fft_lib.fft(low_snr_arr)
                        
                        # Zero-pad the high frequencies to match the highband grid dimensions
                        Y_padded = np.zeros(N_high, dtype=np.complex64)
                        mid = N_low // 2
                        
                        if N_low % 2 == 0:
                            # Even length: split the Nyquist component symmetrically
                            Y_padded[:mid] = Y[:mid]
                            Y_padded[mid] = Y[mid] * 0.5
                            Y_padded[mid + (N_high - N_low)] = Y[mid] * 0.5
                            Y_padded[mid + (N_high - N_low) + 1:] = Y[mid + 1:]
                        else:
                            # Odd length
                            mid_odd = (N_low + 1) // 2
                            Y_padded[:mid_odd] = Y[:mid_odd]
                            Y_padded[mid_odd + (N_high - N_low):] = Y[mid_odd:]
                            
                        # Inverse FFT back to time domain (applying energy conservation scaling)
                        low_snr_upsampled = self.fft_lib.ifft(Y_padded) * (float(N_high) / N_low)
                        
                        # Extract the exact matching slice and add coherently into the highband view
                        current_corr_view[b, buf_slice_start:buf_end] += low_snr_upsampled[slice_start:slice_end]
                # ------------------------------------------

                # -----------------------------------------------------------------
                # --- NEW DROP-IN 1: Prepare 2D Lowband Array for Cython Batching ---
                # -----------------------------------------------------------------
                import time
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
                    #t2 = time.time()
                    
                    if False:         
                        f_list2, t_list2, s_list2 = find_peaks_in_block_cython(
                            current_corr_view, 
                            roi_start,          
                            roi_len,            
                            self.threshold_sq, 
                            f_start,
                            input_offset=buf_slice_start
                        )
                   # t3 = time.time()
                    #print(t3-t2, t2-t1)

                if f_list:
                    all_f_idxs.extend(f_list)
                    all_t_idxs.extend(t_list)
                    all_snrs.extend(s_list)
                    all_tstarts.extend([t_start] * len(s_list)) 
                    
        return (np.array(all_f_idxs, dtype=np.int32), 
                np.array(all_t_idxs, dtype=np.int64), 
                np.array(all_snrs, dtype=np.complex64),
                np.array(all_tstarts, dtype=np.int32))
