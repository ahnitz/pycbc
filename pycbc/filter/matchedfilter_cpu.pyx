# Copyright (C) 2018  Alex Nitz, Josh Willis
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
# cython: embedsignature=True
import numpy
from .matchedfilter import _BaseCorrelator
cimport numpy, cython
from cython.parallel import prange
from libc.math cimport sqrt

# --- Typedefs for Ratio Filter Kernels ---
ctypedef numpy.complex64_t complex64_t
ctypedef numpy.float32_t float32_t
ctypedef numpy.int32_t int32_t
ctypedef numpy.int64_t int64_t

ctypedef fused COMPLEXTYPE:
    float complex
    double complex

@cython.boundscheck(False)
@cython.wraparound(False)
def _batch_correlate(numpy.ndarray [long, ndim=1] x,
                     numpy.ndarray [float complex, ndim=1] y,
                     numpy.ndarray [long, ndim=1] z,
                     size, num_vectors):
    cdef unsigned int nvec = num_vectors
    cdef unsigned int vsize = size

    cdef float complex* xp
    cdef float complex* zp

    cdef unsigned int i, j

    for i in prange(nvec, nogil=True):
        xp = <float complex*> x[i]
        zp = <float complex*> z[i]
        for j in range(vsize):
            zp[j] = xp[j].conjugate() * y[j]

def batch_correlate_execute(self, y):
    num_vectors = self.num_vectors # pylint:disable=unused-variable
    size = self.size # pylint:disable=unused-variable
    _batch_correlate(self.x.data, y.data, self.z.data, size, num_vectors)

def correlate_numpy(x, y, z):
    z.data[:] = numpy.conjugate(x.data)[:]
    z *= y

@cython.boundscheck(False)
@cython.wraparound(False)
def _correlate(COMPLEXTYPE[:] x,
               COMPLEXTYPE[:] y,
               COMPLEXTYPE[:] z):
    cdef unsigned int xmax = x.shape[0]
    cdef unsigned int i
    for i in prange(xmax, nogil=True):
        z[i] = x[i].conjugate() * y[i]

def correlate(x, y, z):
    _correlate(x.data, y.data, z.data)

class CPUCorrelator(_BaseCorrelator):
    def __init__(self, x, y, z):
        self.x = numpy.array(x.data, copy=False)
        self.y = numpy.array(y.data, copy=False)
        self.z = numpy.array(z.data, copy=False)

    def correlate(self):
        _correlate(self.x, self.y, self.z)

def _correlate_factory(x, y, z):
    return CPUCorrelator

# -----------------------------------------------------------------------------
# Ratio Filter Optimization Kernels
# -----------------------------------------------------------------------------

@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)
@cython.initializedcheck(False)
def fast_multiply_indexed_cython(
    numpy.complex64_t[::1] data_f,           
    numpy.complex64_t[:, ::1] filters_master,
    numpy.int64_t[::1] active_idxs,          # The list of surviving filter indices
    numpy.complex64_t[:, ::1] out_batch      
):
    """
    Pulls non-contiguous rows directly from the master filter bank 
    using a list of indices. Computes the analytic half-signal.
    """
    cdef long n_active = active_idxs.shape[0]
    cdef long n_fft = filters_master.shape[1]
    cdef long n_half_plus_one = (n_fft // 2) + 1
    
    cdef long i, j, row_idx

    for i in range(n_active):
        row_idx = active_idxs[i]
        for j in range(n_half_plus_one):
            out_batch[i, j] = data_f[j] * filters_master[row_idx, j]

@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)   
def fast_multiply_analytic_cython(
    numpy.ndarray[complex64_t, ndim=1, mode="c"] data_f,
    numpy.ndarray[complex64_t, ndim=2, mode="c"] filter_batch_f,
    numpy.ndarray[complex64_t, ndim=2, mode="c"] out_batch
):
    """
    Cython version of the "half-only" analytic signal multiply.
    
    This kernel is single-threaded and relies on the C compiler's
    autovectorizer (enabled by -march=native) to use AVX.
    """
    
    # --- C-level variable declarations ---
    cdef long batch_size = filter_batch_f.shape[0]
    cdef long n_fft = filter_batch_f.shape[1]
    
    # We only compute up to the Nyquist bin
    cdef long n_half_plus_one = (n_fft // 2) + 1 
    
    cdef long i, j # Loop iterators

    # This is a pure C-loop, no Python overhead.
    for i in range(batch_size):
        for j in range(n_half_plus_one):
            # Direct C-level complex multiplication
            out_batch[i, j] = data_f[j] * filter_batch_f[i, j]

@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)
@cython.initializedcheck(False)
def fast_multiply_scale_cython(
    numpy.complex64_t[::1] data_f,           # Modern Memoryview, strictly contiguous
    numpy.complex64_t[:, ::1] filter_batch_f,# Modern Memoryview, contiguous in last dimension
    numpy.float32_t[::1] scale,              # Modern Memoryview, strictly contiguous
    numpy.complex64_t[:, ::1] out_batch      # Modern Memoryview, contiguous in last dimension
):
    """
    Cython version of the "half-only" analytic signal multiply.
    
    This kernel is single-threaded and relies on the C compiler's
    autovectorizer (enabled by -march=native) to use AVX.
    """
    
    cdef long batch_size = filter_batch_f.shape[0]
    cdef long n_fft = filter_batch_f.shape[1]
    
    cdef long i, j 
    cdef float current_scale

    for i in range(batch_size):
        current_scale = scale[i] # Explicit loop hoisting
        for j in range(n_fft):
            out_batch[i, j] = current_scale * data_f[j] * filter_batch_f[i, j].conjugate()
            
@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)   
def find_peaks_in_block_cython(
    numpy.ndarray[complex64_t, ndim=2, mode="c"] corr_output,
    long t_start,
    long n_valid,
    float threshold_sq,
    long f_start_offset,
    long input_offset=0
):
    """
    Cython version of the manually vectorized "max-reduction" kernel.
    Returns three separate, flat lists (f_idx, t_idx, snr)
    """
    
    cdef long n_filters_in_batch = corr_output.shape[0]
    cdef list f_idx_list = []
    cdef list t_idx_list = []
    cdef list snr_list = []
    cdef int VEC_WIDTH = 8
    
    cdef float32_t current_max_snr_sq_vec[8]
    cdef int64_t current_max_idx_vec[8]
    cdef complex64_t current_max_z_vec[8]
    
    cdef long f_batch_idx, i, idx, read_idx
    cdef int v_lane
    cdef int32_t f_global_idx
    cdef complex64_t z
    cdef float32_t mag_sq
    cdef float32_t final_max_snr_sq
    cdef int64_t final_max_idx
    
    for f_batch_idx in range(n_filters_in_batch):
        f_global_idx = <int32_t>(f_start_offset + f_batch_idx)

        # --- Initialize C stack arrays ---
        for v_lane in range(VEC_WIDTH):
            current_max_snr_sq_vec[v_lane] = threshold_sq
            current_max_idx_vec[v_lane] = -1

        # --- Main vectorized loop ---
        for i in range(n_valid // VEC_WIDTH):
            for v_lane in range(VEC_WIDTH):
                idx = i * VEC_WIDTH + v_lane
                
                # Apply Offset Here
                read_idx = idx + input_offset
                
                z = corr_output[f_batch_idx, read_idx]
                mag_sq = z.real * z.real + z.imag * z.imag

                if mag_sq > current_max_snr_sq_vec[v_lane]:
                    current_max_snr_sq_vec[v_lane] = mag_sq
                    # Return Global Time Index (t_start corresponds to idx=0)
                    current_max_idx_vec[v_lane] = t_start + idx
                    current_max_z_vec[v_lane] = z

        # --- Epilogue ---
        for i in range((n_valid // VEC_WIDTH) * VEC_WIDTH, n_valid):
            read_idx = i + input_offset
            z = corr_output[f_batch_idx, read_idx]
            mag_sq = z.real * z.real + z.imag * z.imag
            
            if mag_sq > current_max_snr_sq_vec[0]:
                current_max_snr_sq_vec[0] = mag_sq
                current_max_idx_vec[0] = t_start + i
                current_max_z_vec[0] = z
        
        # --- Final Reduction ---
        final_max_snr_sq = threshold_sq
        final_max_idx = -1
        final_max_z = 0 + 0j
        for v_lane in range(VEC_WIDTH):
            if current_max_snr_sq_vec[v_lane] > final_max_snr_sq:
                final_max_snr_sq = current_max_snr_sq_vec[v_lane]
                final_max_idx = current_max_idx_vec[v_lane]
                final_max_z = current_max_z_vec[v_lane]
        
        if final_max_idx != -1:
            f_idx_list.append(f_global_idx)
            t_idx_list.append(final_max_idx)
            snr_list.append(final_max_z)
            
    return (f_idx_list, t_idx_list, snr_list)


import numpy as np
cimport numpy as numpy
cimport cython
cimport libc.math
import logging

# =====================================================================
# GLOBAL MODULE LEVEL: Pregenerate Static Kaiser-128 Weights LUT
# =====================================================================
_w_grid = np.linspace(0, 1, 1024, endpoint=False)[:, None]
_k_grid = np.arange(-63, 65)[None, :]  # 128 taps total
_t_grid = _w_grid - _k_grid

# Kaiser Window setup for 128 taps (N_half = 64)
beta = 4.0
_x = _t_grid / 64.0
_x_sq = np.clip(1.0 - _x**2, 0.0, 1.0)
_kaiser_window = np.i0(beta * np.sqrt(_x_sq)) / np.i0(beta)

_ideal_sinc = np.sinc(_t_grid)
_lut_raw = _ideal_sinc * _kaiser_window
_lut_raw /= _lut_raw.sum(axis=1, keepdims=True)

KAISER_128_LUT = np.ascontiguousarray(_lut_raw, dtype=np.float32)


@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)    
@cython.initializedcheck(False)
def find_peaks_indexed_cython(
    numpy.complex64_t[:, ::1] corr_output,     # Dense 2D output from the IFFT
    numpy.complex64_t[:, ::1] lowband_master,  # The FULL 2D master lowband array
    numpy.int64_t[::1] active_idxs,            # 1D array of global filter indices
    long t_start,
    long n_valid,
    float threshold_sq,
    double dt_high,
    double dt_low,
    double t0_high,
    double t0_low,
    long input_offset=0
):
    """
    Two-Pass Fused Kernel utilizing pristine block-level max scans.
    Reads dense highband data, but sparsely looks up lowband data using active_idxs.
    """
    cdef long n_filters_in_batch = corr_output.shape[0]
    cdef list f_idx_list = []
    cdef list t_idx_list = []
    cdef list snr_list = []
    cdef list snr_low_list = []
    
    cdef long f_batch_idx, j, idx, read_idx
    cdef int k_idx
    cdef int32_t f_global_idx
    
    cdef complex64_t z_high, z_total, z_low0, z_low1, interp_low_snr
    cdef float32_t mag_sq, high_mag_bound, low_mag0_bound, low_mag1_bound, low_envelope_bound
    
    cdef float32_t current_max_snr_sq, current_max_snr_amp 
    cdef int64_t current_max_idx
    cdef complex64_t current_max_z, current_max_z_low
    
    cdef float[:, :] weights_lut = KAISER_128_LUT
    
    cdef double const_offset = (t0_high - t0_low) + <double>t_start * dt_high
    cdef double data_scale = dt_high / dt_low
    cdef double base_scale = const_offset / dt_low
    
    cdef double low_idx_start = const_offset / dt_low
    cdef double low_idx_end = (const_offset + <double>(n_valid - 1) * dt_high) / dt_low
    cdef long j_start = <long>libc.math.floor(low_idx_start)
    cdef long j_end = <long>libc.math.floor(low_idx_end)
    
    cdef long idx_min_j, idx_max_j
    cdef double exact_low_idx, w
    cdef int w_lut_idx
    
    cdef double PI = 3.14159265358979323846
    cdef double theta, c_rot, s_rot, rot_r, rot_i
    cdef float sign
    cdef float32_t max_low_envelope, l_amp, block_max_high_sq, block_max_high_amp
    
    for f_batch_idx in range(n_filters_in_batch):
        # Maps the local batch row to the true global template identity
        f_global_idx = <int32_t>active_idxs[f_batch_idx]

        current_max_snr_sq = threshold_sq
        current_max_snr_amp = <float>libc.math.sqrt(threshold_sq) 
        current_max_idx = -1
        current_max_z.real = 0.0
        current_max_z.imag = 0.0
        current_max_z_low.real = 0.0
        current_max_z_low.imag = 0.0

        # === PASS 1a: Scan Lowband ===
        max_low_envelope = 0.0
        for j in range(j_start, j_end + 2):
            # NEW: We look up using f_global_idx instead of f_batch_idx
            z_low0 = lowband_master[f_global_idx, j]
            l_amp = <float32_t>libc.math.sqrt(z_low0.real * z_low0.real + z_low0.imag * z_low0.imag)
            if l_amp > max_low_envelope:
                max_low_envelope = l_amp
        max_low_envelope *= 1.25
        
        # === PASS 1b: Unbroken Highband Max Scan ===
        block_max_high_sq = 0.0
        for idx in range(n_valid):
            read_idx = idx + input_offset
            # HIGHBAND uses f_batch_idx because it's dense
            z_high = corr_output[f_batch_idx, read_idx]
            mag_sq = z_high.real * z_high.real + z_high.imag * z_high.imag
            if mag_sq > block_max_high_sq:
                block_max_high_sq = mag_sq
                
        block_max_high_amp = <float32_t>libc.math.sqrt(block_max_high_sq)

        if (block_max_high_amp + max_low_envelope) <= current_max_snr_amp:
            continue

        # === PASS 2: Fallback Detailed Segment Loop ===
        for j in range(j_start, j_end + 1):
            idx_min_j = <long>libc.math.ceil((<double>j * dt_low - const_offset) / dt_high)
            idx_max_j = <long>libc.math.ceil((<double>(j + 1) * dt_low - const_offset) / dt_high)
            
            if idx_min_j < 0: idx_min_j = 0
            if idx_max_j > n_valid: idx_max_j = n_valid
            if idx_min_j >= idx_max_j: continue
            
            # LOWBAND uses f_global_idx
            z_low0 = lowband_master[f_global_idx, j]
            z_low1 = lowband_master[f_global_idx, j + 1]
            low_mag0_bound = libc.math.fabs(z_low0.real) + libc.math.fabs(z_low0.imag)
            low_mag1_bound = libc.math.fabs(z_low1.real) + libc.math.fabs(z_low1.imag)
            low_envelope_bound = (low_mag0_bound if low_mag0_bound > low_mag1_bound else low_mag1_bound) * 1.15
            
            for idx in range(idx_min_j, idx_max_j):
                read_idx = idx + input_offset
                
                # HIGHBAND uses f_batch_idx
                z_high = corr_output[f_batch_idx, read_idx]
                high_mag_bound = libc.math.fabs(z_high.real) + libc.math.fabs(z_high.imag)
                
                if (high_mag_bound + low_envelope_bound) <= current_max_snr_amp:
                    continue
                
                exact_low_idx = base_scale + <double>idx * data_scale
                w = exact_low_idx - <double>j
                
                w_lut_idx = <int>(w * 1024.0)
                if w_lut_idx > 1023: w_lut_idx = 1023
                elif w_lut_idx < 0: w_lut_idx = 0
                
                interp_low_snr.real = 0.0
                interp_low_snr.imag = 0.0
                
                sign = 1.0 if ((j - 63) % 2) == 0 else -1.0
                
                for k_idx in range(128):
                    # LOWBAND uses f_global_idx
                    interp_low_snr.real += sign * weights_lut[w_lut_idx, k_idx] * lowband_master[f_global_idx, j + k_idx - 63].real
                    interp_low_snr.imag += sign * weights_lut[w_lut_idx, k_idx] * lowband_master[f_global_idx, j + k_idx - 63].imag
                    sign = -sign
                
                theta = PI * exact_low_idx
                c_rot = libc.math.cos(theta)
                s_rot = libc.math.sin(theta)
                
                rot_r = interp_low_snr.real * c_rot - interp_low_snr.imag * s_rot
                rot_i = interp_low_snr.real * s_rot + interp_low_snr.imag * c_rot
                
                interp_low_snr.real = rot_r
                interp_low_snr.imag = rot_i
                
                z_total.real = z_high.real + interp_low_snr.real
                z_total.imag = z_high.imag + interp_low_snr.imag
                mag_sq = z_total.real * z_total.real + z_total.imag * z_total.imag
                
                if mag_sq > current_max_snr_sq:
                    current_max_snr_sq = mag_sq
                    current_max_snr_amp = <float>libc.math.sqrt(mag_sq)
                    current_max_idx = t_start + idx
                    current_max_z = z_total
                    current_max_z_low = interp_low_snr
                    
        if current_max_idx != -1:
            f_idx_list.append(f_global_idx)
            t_idx_list.append(current_max_idx)
            snr_list.append(current_max_z)
            snr_low_list.append(current_max_z_low)
            
    return (f_idx_list, t_idx_list, snr_list, snr_low_list)


@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.cdivision(True)    
def find_peaks_fused_lanczos_cython(
    numpy.ndarray[complex64_t, ndim=2, mode="c"] corr_output,
    numpy.ndarray[complex64_t, ndim=2, mode="c"] low_snr_block,
    long t_start,
    long n_valid,
    float threshold_sq,
    long f_start_offset,
    double dt_high,
    double dt_low,
    double t0_high,
    double t0_low,
    long input_offset=0
):
    """
    Two-Pass Fused Kernel utilizing pristine block-level max scans.
    UPDATED: 128-Tap Kaiser-Sinc for ultra-low SNR sub-sample stability.
    """
    cdef long n_filters_in_batch = corr_output.shape[0]
    cdef list f_idx_list = []
    cdef list t_idx_list = []
    cdef list snr_list = []
    cdef list snr_low_list = []
    
    cdef long f_batch_idx, j, idx, read_idx
    cdef int k_idx
    cdef int32_t f_global_idx
    
    cdef complex64_t z_high, z_total, z_low0, z_low1, interp_low_snr
    cdef float32_t mag_sq, high_mag_bound, low_mag0_bound, low_mag1_bound, low_envelope_bound
    
    cdef float32_t current_max_snr_sq
    cdef float32_t current_max_snr_amp 
    cdef int64_t current_max_idx
    cdef complex64_t current_max_z
    cdef complex64_t current_max_z_low
    
    cdef float[:, :] weights_lut = KAISER_128_LUT
    
    # Grid transformation constants
    cdef double const_offset = (t0_high - t0_low) + <double>t_start * dt_high
    cdef double data_scale = dt_high / dt_low
    cdef double base_scale = const_offset / dt_low
    
    cdef double low_idx_start = const_offset / dt_low
    cdef double low_idx_end = (const_offset + <double>(n_valid - 1) * dt_high) / dt_low
    cdef long j_start = <long>libc.math.floor(low_idx_start)
    cdef long j_end = <long>libc.math.floor(low_idx_end)
    
    cdef long idx_min_j, idx_max_j
    cdef double exact_low_idx, w
    cdef int w_lut_idx
    
    # Analytic Lanczos Variables
    cdef double PI = 3.14159265358979323846
    cdef double theta, c_rot, s_rot, rot_r, rot_i
    cdef float sign
    
    # Pass 1 Optimization Scans
    cdef float32_t max_low_envelope, l_amp, block_max_high_sq, block_max_high_amp
    
    for f_batch_idx in range(n_filters_in_batch):
        f_global_idx = <int32_t>(f_start_offset + f_batch_idx)

        current_max_snr_sq = threshold_sq
        current_max_snr_amp = <float>libc.math.sqrt(threshold_sq) 
        current_max_idx = -1
        current_max_z.real = 0.0
        current_max_z.imag = 0.0
        current_max_z_low.real = 0.0
        current_max_z_low.imag = 0.0

        # =====================================================================
        # PASS 1a: Scan Short Lowband Segment for Max Amplitude
        # =====================================================================
        max_low_envelope = 0.0
        for j in range(j_start, j_end + 2):
            z_low0 = low_snr_block[f_batch_idx, j]
            l_amp = <float32_t>libc.math.sqrt(z_low0.real * z_low0.real + z_low0.imag * z_low0.imag)
            if l_amp > max_low_envelope:
                max_low_envelope = l_amp
        max_low_envelope *= 1.25
        
        # =====================================================================
        # PASS 1b: Unbroken Highband Max Scan
        # =====================================================================
        block_max_high_sq = 0.0
        for idx in range(n_valid):
            read_idx = idx + input_offset
            z_high = corr_output[f_batch_idx, read_idx]
            mag_sq = z_high.real * z_high.real + z_high.imag * z_high.imag
            if mag_sq > block_max_high_sq:
                block_max_high_sq = mag_sq
                
        block_max_high_amp = <float32_t>libc.math.sqrt(block_max_high_sq)

        # =====================================================================
        # GLOBAL BLOCK GATE CHECK
        # =====================================================================
        if (block_max_high_amp + max_low_envelope) <= current_max_snr_amp:
            continue

        # =====================================================================
        # PASS 2: Fallback Detailed Segment Loop
        # =====================================================================
        for j in range(j_start, j_end + 1):
            idx_min_j = <long>libc.math.ceil((<double>j * dt_low - const_offset) / dt_high)
            idx_max_j = <long>libc.math.ceil((<double>(j + 1) * dt_low - const_offset) / dt_high)
            
            if idx_min_j < 0: idx_min_j = 0
            if idx_max_j > n_valid: idx_max_j = n_valid
            if idx_min_j >= idx_max_j: continue
            
            z_low0 = low_snr_block[f_batch_idx, j]
            z_low1 = low_snr_block[f_batch_idx, j + 1]
            low_mag0_bound = libc.math.fabs(z_low0.real) + libc.math.fabs(z_low0.imag)
            low_mag1_bound = libc.math.fabs(z_low1.real) + libc.math.fabs(z_low1.imag)
            low_envelope_bound = (low_mag0_bound if low_mag0_bound > low_mag1_bound else low_mag1_bound) * 1.15
            
            for idx in range(idx_min_j, idx_max_j):
                read_idx = idx + input_offset
                z_high = corr_output[f_batch_idx, read_idx]
                high_mag_bound = libc.math.fabs(z_high.real) + libc.math.fabs(z_high.imag)
                
                if (high_mag_bound + low_envelope_bound) <= current_max_snr_amp:
                    continue
                
                exact_low_idx = base_scale + <double>idx * data_scale
                w = exact_low_idx - <double>j
                
                w_lut_idx = <int>(w * 1024.0)
                if w_lut_idx > 1023: w_lut_idx = 1023
                elif w_lut_idx < 0: w_lut_idx = 0
                
                interp_low_snr.real = 0.0
                interp_low_snr.imag = 0.0
                
                # --- EXACT PHASE SIGN ALIGNMENT (128 Taps) ---
                # Based on center tap offset of -63
                sign = 1.0 if ((j - 63) % 2) == 0 else -1.0
                
                for k_idx in range(128):
                    interp_low_snr.real += sign * weights_lut[w_lut_idx, k_idx] * low_snr_block[f_batch_idx, j + k_idx - 63].real
                    interp_low_snr.imag += sign * weights_lut[w_lut_idx, k_idx] * low_snr_block[f_batch_idx, j + k_idx - 63].imag
                    sign = -sign
                
                # --- PHASE ROTATION ---
                theta = PI * exact_low_idx
                c_rot = libc.math.cos(theta)
                s_rot = libc.math.sin(theta)
                
                rot_r = interp_low_snr.real * c_rot - interp_low_snr.imag * s_rot
                rot_i = interp_low_snr.real * s_rot + interp_low_snr.imag * c_rot
                
                interp_low_snr.real = rot_r
                interp_low_snr.imag = rot_i
                
                z_total.real = z_high.real + interp_low_snr.real
                z_total.imag = z_high.imag + interp_low_snr.imag
                mag_sq = z_total.real * z_total.real + z_total.imag * z_total.imag
                
                if mag_sq > current_max_snr_sq:
                    current_max_snr_sq = mag_sq
                    current_max_snr_amp = <float>libc.math.sqrt(mag_sq)
                    current_max_idx = t_start + idx
                    current_max_z = z_total
                    current_max_z_low = interp_low_snr
                    
        if current_max_idx != -1:
            f_idx_list.append(f_global_idx)
            t_idx_list.append(current_max_idx)
            snr_list.append(current_max_z)
            snr_low_list.append(current_max_z_low)
            
    return (f_idx_list, t_idx_list, snr_list, snr_low_list)
    
    
import numpy as np
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def block_max_threshold(
    float complex[:, ::1] snrs,   
    int[:] starts, 
    int[:] ends, 
    unsigned char[:, ::1] valid,       
    double[:] threshold,
):
    cdef int num_blocks = starts.shape[0]
    cdef int num_rows = snrs.shape[0]
    cdef int j, i, k, js, je
    
    cdef float val_real, val_imag, mag_sq, local_max_sq

    cdef float threshold_sq = 0

    with nogil:
        for i in range(num_rows):
            threshold_sq = threshold[i] * threshold[i]
            for j in range(num_blocks):
                
                if valid[i, j] == 1:
                    js = starts[j]
                    je = ends[j]
                    
                    # Reset the max tracker for this block
                    local_max_sq = -1.0
                    
                    # ---------------------------------------------------
                    # PURE MATH INNER LOOP (No Branches, No Breaks)
                    # The compiler will convert this to SIMD AVX instructions
                    # ---------------------------------------------------
                    for k in range(js, je):
                        val_real = snrs[i, k].real
                        val_imag = snrs[i, k].imag
                        mag_sq = val_real * val_real + val_imag * val_imag
                        
                        # A simple max assignment is easily vectorized
                        if mag_sq > local_max_sq:
                            local_max_sq = mag_sq

                    # Check threshold once at the very end
                    if local_max_sq < threshold_sq:
                        valid[i, j] = 0
