# Copyright (C) 2016 Collin Capano
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
""" Functions for applying gates to data.
"""

import numpy as np
from scipy import linalg
from . import strain


def _gates_from_cli(opts, gate_opt):
    """Parses the given `gate_opt` into something understandable by
    `strain.gate_data`.
    """
    gates = {}
    if getattr(opts, gate_opt) is None:
        return gates
    for gate in getattr(opts, gate_opt):
        try:
            ifo, central_time, half_dur, taper_dur = gate.split(':')
            central_time = float(central_time)
            half_dur = float(half_dur)
            taper_dur = float(taper_dur)
        except ValueError:
            raise ValueError("--gate {} not formatted correctly; ".format(
                gate) + "see help")
        try:
            gates[ifo].append((central_time, half_dur, taper_dur))
        except KeyError:
            gates[ifo] = [(central_time, half_dur, taper_dur)]
    return gates


def gates_from_cli(opts):
    """Parses the --gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, 'gate')


def psd_gates_from_cli(opts):
    """Parses the --psd-gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, 'psd_gate')


def apply_gates_to_td(strain_dict, gates):
    """Applies the given dictionary of gates to the given dictionary of
    strain.

    Parameters
    ----------
    strain_dict : dict
        Dictionary of time-domain strain, keyed by the ifos.
    gates : dict
        Dictionary of gates. Keys should be the ifo to apply the data to,
        values are a tuple giving the central time of the gate, the half
        duration, and the taper duration.

    Returns
    -------
    dict
        Dictionary of time-domain strain with the gates applied.
    """
    # copy data to new dictionary
    outdict = dict(strain_dict.items())
    for ifo in gates:
        outdict[ifo] = strain.gate_data(outdict[ifo], gates[ifo])
    return outdict


def apply_gates_to_fd(stilde_dict, gates):
    """Applies the given dictionary of gates to the given dictionary of
    strain in the frequency domain.

    Gates are applied by IFFT-ing the strain data to the time domain, applying
    the gate, then FFT-ing back to the frequency domain.

    Parameters
    ----------
    stilde_dict : dict
        Dictionary of frequency-domain strain, keyed by the ifos.
    gates : dict
        Dictionary of gates. Keys should be the ifo to apply the data to,
        values are a tuple giving the central time of the gate, the half
        duration, and the taper duration.

    Returns
    -------
    dict
        Dictionary of frequency-domain strain with the gates applied.
    """
    # copy data to new dictionary
    outdict = dict(stilde_dict.items())
    # create a time-domin strain dictionary to apply the gates to
    strain_dict = dict([[ifo, outdict[ifo].to_timeseries()] for ifo in gates])
    # apply gates and fft back to the frequency domain
    for ifo,d in apply_gates_to_td(strain_dict, gates).items():
        outdict[ifo] = d.to_frequencyseries()
    return outdict


def add_gate_option_group(parser):
    """Adds the options needed to apply gates to data.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    gate_group = parser.add_argument_group("Options for gating data")

    gate_group.add_argument("--gate", nargs="+", type=str,
                            metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
                            help="Apply one or more gates to the data before "
                                 "filtering.")
    gate_group.add_argument("--gate-overwhitened", action="store_true",
                            help="Overwhiten data first, then apply the "
                                 "gates specified in --gate. Overwhitening "
                                 "allows for sharper tapers to be used, "
                                 "since lines are not blurred.")
    gate_group.add_argument("--psd-gate", nargs="+", type=str,
                            metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
                            help="Apply one or more gates to the data used "
                                 "for computing the PSD. Gates are applied "
                                 "prior to FFT-ing the data for PSD "
                                 "estimation.")
    return gate_group


def _get_tdfilter(invpsd, max_samples=65536):
    """Retrieve or compute cached time-domain impulse response of inverse PSD."""
    cached = getattr(invpsd, '_tdfilter', None)
    if cached is not None:
        return cached

    flen = len(invpsd)
    full_N = (flen - 1) * 2
    if full_N <= max_samples:
        tdfilter = (invpsd.astype('complex').to_timeseries() * invpsd.delta_t).numpy()
    else:
        from pycbc.types import FrequencySeries
        step = (flen - 1) // (max_samples // 2)
        invpsd_sub = FrequencySeries(invpsd[::step].numpy(), delta_f=invpsd.delta_f * step)
        tdfilter = (invpsd_sub.astype('complex').to_timeseries() * invpsd_sub.delta_t).numpy()

    invpsd._tdfilter = tdfilter
    return tdfilter


def gate_and_paint(data, lindex, rindex, invpsd, copy=True, method='toeplitz',
                   ridge=1e-10, pad_duration=16.0):
    """Gates and in-paints data using a hole-filling solver.

    Parameters
    ----------
    data : TimeSeries
        The data to gate.
    lindex : int
        The start index of the gate.
    rindex : int
        The end index of the gate.
    invpsd : FrequencySeries
        The inverse of the PSD.
    copy : bool, optional
        Copy the data before applying the gate. Otherwise, the gate will
        be applied in-place. Default is True.
    method : {'toeplitz', 'cholesky'}, optional
        Algorithm used to solve the linear system for the inpainting projection.
        'toeplitz' (default) uses regularized scipy.linalg.solve_toeplitz (Levinson recursion, O(K^2)).
        'cholesky' uses a regularized Cholesky factorization (O(K^3)).
    ridge : float, optional
        Diagonal Tikhonov regularization parameter relative to the diagonal
        element of the inverse covariance operator. Default is 1e-3.
    pad_duration : float, optional
        Padding duration in seconds on each side of the gate for local
        convolution inpainting. Default is 16.0 seconds. If None or if data
        duration is too short, falls back to full-array global FFT.

    Returns
    -------
    TimeSeries :
        The gated and in-painted time series.
    """
    # Uses the hole-filling method of
    # https://arxiv.org/pdf/1908.05644.pdf
    if copy:
        data = data.copy()
    data[lindex:rindex] = 0.0
    K = rindex - lindex
    if hasattr(invpsd, 'precision') and hasattr(data, 'precision') and invpsd.precision != data.precision:
        invpsd = invpsd.astype(np.float32 if data.precision == 'single' else np.float64)

    tdfilter = _get_tdfilter(invpsd)
    diag = tdfilter[0]
    delta_t = float(data.delta_t)
    M = min(len(tdfilter) // 2, int(round(pad_duration / delta_t))) if pad_duration is not None else 0

    if M > 0 and len(data) > 2 * M + K:
        from scipy import signal
        idx_start = lindex - M
        idx_end = rindex + M
        pad_left = 0
        pad_right = 0
        if idx_start < 0:
            pad_left = -idx_start
            idx_start = 0
        if idx_end > len(data):
            pad_right = idx_end - len(data)
            idx_end = len(data)

        data_chunk = data[idx_start:idx_end].numpy()
        if pad_left > 0 or pad_right > 0:
            data_chunk = np.pad(data_chunk, (pad_left, pad_right), mode='constant')

        filt_neg = tdfilter[-M:]      # lags -M to -1
        filt_pos = tdfilter[:M+1]     # lags 0 to M
        kernel = np.concatenate([filt_neg, filt_pos])
        rhs = signal.fftconvolve(data_chunk, kernel[::-1], mode='valid')
    else:
        owhgated_data = (data.to_frequencyseries() * invpsd).to_timeseries()
        rhs = owhgated_data[lindex:rindex].numpy()

    if method == 'cholesky':
        col = tdfilter[:K].numpy() if hasattr(tdfilter[:K], 'numpy') else tdfilter[:K]
        col = col / diag
        T = linalg.toeplitz(col)
        if ridge > 0:
            T += ridge * np.eye(K)
        c, lower = linalg.cho_factor(T)
        proj = linalg.cho_solve((c, lower), rhs / diag)
    elif method == 'toeplitz':
        col = tdfilter[:K].numpy().copy() if hasattr(tdfilter[:K], 'numpy') else tdfilter[:K].copy()
        if ridge > 0:
            col[0] += ridge * diag
        proj = linalg.solve_toeplitz(col, rhs)
    else:
        raise ValueError(f"Unknown inpainting method: {method}")

    # remove the projection into the null space
    data[lindex:rindex] -= proj
    return data


def invert_covariance(invpsd, lindex, rindex, ridge=1e-10):
    """Calculate the uninverted covariance matrix with optional regularization.

    Parameters
    ----------
    invpsd : FrequencySeries
        The inverse of the PSD.
    lindex : int
        The start index of the gate.
    rindex : int
        The end index of the gate.
    ridge : float, optional
        Regularization ridge parameter. Default is 1e-10.

    Returns
    -------
    array :
        The uninverted covariance matrix associated with the inverse PSD in the
        time window [lindex, rindex].
    """
    K = rindex - lindex
    tdfilter = _get_tdfilter(invpsd)
    diag = tdfilter[0]
    col = tdfilter[:K].numpy() if hasattr(tdfilter[:K], 'numpy') else tdfilter[:K]
    mat = linalg.toeplitz(col / diag)
    if ridge > 0:
        mat += ridge * np.eye(K)
    invmat = linalg.inv(mat) / diag
    return invmat


def gate_and_paint_matmul(data, lindex, rindex, invpsd, invmat=None, ridge=1e-10, copy=True, pad_duration=16.0):
    """Gates and in-paints data using explicit matrix multiplication.

    Parameters
    ----------
    data : TimeSeries
        The data to gate.
    lindex : int
        The start index of the gate.
    rindex : int
        The end index of the gate.
    invpsd : FrequencySeries
        The inverse of the PSD.
    invmat : array, optional
        The uninverted covariance matrix. If None, calculate on function call.
    ridge : float, optional
        Regularization ridge parameter if invmat is calculated. Default is 1e-10.
    copy : bool, optional
        Copy the data before applying the gate. Otherwise, the gate will
        be applied in-place. Default is True.
    pad_duration : float, optional
        Padding duration in seconds on each side of the gate for local
        convolution inpainting. Default is 16.0 seconds.

    Returns
    -------
    TimeSeries :
        The gated and in-painted time series.
    """
    if copy:
        data = data.copy()
    data[lindex:rindex] = 0.0
    K = rindex - lindex

    if hasattr(invpsd, 'precision') and hasattr(data, 'precision') and invpsd.precision != data.precision:
        invpsd = invpsd.astype(np.float32 if data.precision == 'single' else np.float64)

    tdfilter = _get_tdfilter(invpsd)
    delta_t = float(data.delta_t)
    M = min(len(tdfilter) // 2, int(round(pad_duration / delta_t))) if pad_duration is not None else 0

    if M > 0 and len(data) > 2 * M + K:
        from scipy import signal
        idx_start = lindex - M
        idx_end = rindex + M
        pad_left = 0
        pad_right = 0
        if idx_start < 0:
            pad_left = -idx_start
            idx_start = 0
        if idx_end > len(data):
            pad_right = idx_end - len(data)
            idx_end = len(data)

        data_chunk = data[idx_start:idx_end].numpy()
        if pad_left > 0 or pad_right > 0:
            data_chunk = np.pad(data_chunk, (pad_left, pad_right), mode='constant')

        filt_neg = tdfilter[-M:]
        filt_pos = tdfilter[:M+1]
        kernel = np.concatenate([filt_neg, filt_pos])
        rhs = signal.fftconvolve(data_chunk, kernel[::-1], mode='valid')
    else:
        owhgated_data = (data.to_frequencyseries() * invpsd).to_timeseries()
        rhs = owhgated_data[lindex:rindex].numpy()

    # invert the matrix if not provided
    if invmat is None:
        invmat = invert_covariance(invpsd, lindex, rindex, ridge=ridge)

    # remove the projection into the null space
    proj = invmat @ rhs
    data[lindex:rindex] -= proj
    return data
