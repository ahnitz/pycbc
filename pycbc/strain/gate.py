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
"""Functions for applying gates to data."""

import logging

import numpy as np
from scipy import linalg, signal

from pycbc.types import FrequencySeries, TimeSeries

from . import strain

logger = logging.getLogger("pycbc.strain.gate")


def _gates_from_cli(opts, gate_opt):
    """
    Parses the given `gate_opt` into something understandable by
    `strain.gate_data`.
    """
    gates = {}
    if getattr(opts, gate_opt) is None:
        return gates
    for gate in getattr(opts, gate_opt):
        try:
            ifo, central_time, half_dur, taper_dur = gate.split(":")
            central_time = float(central_time)
            half_dur = float(half_dur)
            taper_dur = float(taper_dur)
        except ValueError as err:
            raise ValueError(
                f"--gate {gate} not formatted correctly; " + "see help"
            ) from err
        try:
            gates[ifo].append((central_time, half_dur, taper_dur))
        except KeyError:
            gates[ifo] = [(central_time, half_dur, taper_dur)]
    return gates


def gates_from_cli(opts):
    """
    Parses the --gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, "gate")


def psd_gates_from_cli(opts):
    """
    Parses the --psd-gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, "psd_gate")


def apply_gates_to_td(strain_dict, gates):
    """
    Applies the given dictionary of gates to the given dictionary of
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
    """
    Applies the given dictionary of gates to the given dictionary of
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
    for ifo, d in apply_gates_to_td(strain_dict, gates).items():
        outdict[ifo] = d.to_frequencyseries()
    return outdict


def add_gate_option_group(parser):
    """
    Adds the options needed to apply gates to data.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.

    """
    gate_group = parser.add_argument_group("Options for gating data")

    gate_group.add_argument(
        "--gate",
        nargs="+",
        type=str,
        metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
        help="Apply one or more gates to the data before filtering.",
    )
    gate_group.add_argument(
        "--gate-overwhitened",
        action="store_true",
        help="Overwhiten data first, then apply the "
        "gates specified in --gate. Overwhitening "
        "allows for sharper tapers to be used, "
        "since lines are not blurred.",
    )
    gate_group.add_argument(
        "--psd-gate",
        nargs="+",
        type=str,
        metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
        help="Apply one or more gates to the data used "
        "for computing the PSD. Gates are applied "
        "prior to FFT-ing the data for PSD "
        "estimation.",
    )
    return gate_group


def _get_tdfilter(invpsd, max_samples=1048576):
    """Return a cached time-domain inverse-covariance filter for ``invpsd``."""
    flen = len(invpsd)
    full_len = (flen - 1) * 2
    target_samples = min(full_len, max_samples)
    cached = getattr(invpsd, "_tdfilter", None)
    if cached is not None and len(cached) >= target_samples:
        return cached

    if full_len <= max_samples:
        tdfilter = (invpsd.astype("complex").to_timeseries() * invpsd.delta_t).numpy()
    else:
        step = max(1, (flen - 1) // (max_samples // 2))
        invpsd_sub = FrequencySeries(
            invpsd[::step].numpy(), delta_f=invpsd.delta_f * step
        )
        tdfilter = (
            invpsd_sub.astype("complex").to_timeseries() * invpsd_sub.delta_t
        ).numpy()

    if cached is None or len(tdfilter) > len(cached):
        invpsd._tdfilter = tdfilter
    return tdfilter


def _compute_inpaint_rhs(data, lindex, rindex, invpsd, pad_duration):
    """Compute the inverse-covariance filter and RHS vector inside the hole."""
    hole_len = rindex - lindex
    if (
        hasattr(invpsd, "precision")
        and hasattr(data, "precision")
        and invpsd.precision != data.precision
    ):
        invpsd = invpsd.astype(np.float32 if data.precision == "single" else np.float64)

    req_max_samples = max(1048576, 1 << int(np.ceil(np.log2(max(2 * hole_len, 2)))))
    tdfilter = _get_tdfilter(invpsd, max_samples=req_max_samples)
    pad_samples = (
        min(len(tdfilter) // 2, round(pad_duration / float(data.delta_t)))
        if pad_duration is not None
        else 0
    )

    if pad_samples > 0 and len(data) > 2 * pad_samples + hole_len:
        idx_start = lindex - pad_samples
        idx_end = rindex + pad_samples
        pad_left = max(0, -idx_start)
        pad_right = max(0, idx_end - len(data))
        idx_start = max(0, idx_start)
        idx_end = min(len(data), idx_end)

        data_chunk = data[idx_start:idx_end].numpy()
        if pad_left > 0 or pad_right > 0:
            data_chunk = np.pad(data_chunk, (pad_left, pad_right), mode="constant")

        kernel = np.concatenate([tdfilter[-pad_samples:], tdfilter[: pad_samples + 1]])
        rhs = signal.fftconvolve(
            data_chunk.astype(np.float64),
            kernel[::-1].astype(np.float64),
            mode="valid",
        )
    else:
        owhgated_data = (data.to_frequencyseries() * invpsd).to_timeseries()
        rhs = owhgated_data[lindex:rindex].numpy().astype(np.float64)

    return tdfilter, rhs


def gate_and_paint(
    data,
    lindex,
    rindex,
    invpsd,
    copy=True,
    method="toeplitz",
    ridge=1e-7,
    pad_duration=16.0,
    fallback_taper=True,
):
    """Gate and inpaint data using a regularized hole-filling linear solver.

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
        Copy the data before applying the gate. Otherwise, apply in-place.
        Default is True.
    method : {'toeplitz', 'cholesky', 'auto'}, optional
        Solver used for the inpainting projection. ``'toeplitz'`` (default, or
        ``'auto'``) uses regularized Levinson recursion (O(K^2)).
        ``'cholesky'`` uses regularized Cholesky factorization (O(K^3)), and
        switches automatically to ``'toeplitz'`` when ``rindex - lindex > 2048``.
    ridge : float, optional
        Relative diagonal Tikhonov regularization parameter. Default is 1e-7.
    pad_duration : float, optional
        One-sided padding duration in seconds for local convolution inpainting.
        If None or the data is shorter than the padded gate window, falls back
        to global FFT evaluation. Default is 16.0.
    fallback_taper : bool, optional
        If True (default) and gate width exceeds 32768 samples, fall back to
        cosine taper gating to avoid ill-conditioned, expensive inpainting.

    Returns
    -------
    TimeSeries :
        The gated and in-painted time series.

    """
    # Uses the hole-filling method of https://arxiv.org/pdf/1908.05644.pdf
    if copy:
        data = data.copy()
    hole_len = rindex - lindex
    if hole_len <= 0:
        return data

    # Fall back to cosine tapering on ultra-wide gates (>32768 samples)
    if fallback_taper and hole_len > 32768:
        logger.warning(
            "Gate length K=%d > 32768 samples; falling back to cosine taper gating.",
            hole_len,
        )
        t_center = float(data.start_time + (lindex + rindex) / 2.0 * data.delta_t)
        half_width = float(hole_len / 2.0 * data.delta_t)
        return strain.gate_data(data, [(t_center, half_width, min(0.25, half_width))])

    data[lindex:rindex] = 0.0
    tdfilter, rhs = _compute_inpaint_rhs(data, lindex, rindex, invpsd, pad_duration)
    diag = float(tdfilter[0])
    raw_col = tdfilter[:hole_len]
    if len(raw_col) < hole_len:
        raw_col = np.pad(raw_col, (0, hole_len - len(raw_col)))

    resolved_method = "toeplitz" if method in (None, "auto") else method
    if resolved_method == "cholesky" and hole_len > 2048:
        logger.info(
            "Gate length K=%d > 2048 samples; switching from cholesky to toeplitz.",
            hole_len,
        )
        resolved_method = "toeplitz"

    # Solve the regularized Toeplitz linear system inside the hole
    if resolved_method == "cholesky":
        matrix = linalg.toeplitz(raw_col.astype(np.float64) / diag)
        if ridge > 0:
            matrix += ridge * np.eye(hole_len, dtype=np.float64)
        factor, lower = linalg.cho_factor(matrix)
        proj = linalg.cho_solve((factor, lower), rhs / diag)
    elif resolved_method == "toeplitz":
        col = raw_col.astype(np.float64).copy()
        if ridge > 0:
            col[0] += ridge * diag
        proj = linalg.solve_toeplitz(col, rhs)
    else:
        raise ValueError(f"Unknown inpainting method: {method}")

    # Remove the projection into the null space
    data[lindex:rindex] -= proj.astype(data.dtype)
    return data


def invert_covariance(invpsd, lindex, rindex, ridge=1e-7):
    """Calculate the regularized uninverted covariance matrix inside a gate.

    Parameters
    ----------
    invpsd : FrequencySeries
        The inverse of the PSD.
    lindex : int
        The start index of the gate.
    rindex : int
        The end index of the gate.
    ridge : float, optional
        Relative diagonal Tikhonov regularization parameter. Default is 1e-7.

    Returns
    -------
    array :
        The uninverted covariance matrix associated with the inverse PSD in the
        time window [lindex, rindex].

    """
    hole_len = rindex - lindex
    req_max_samples = max(1048576, 1 << int(np.ceil(np.log2(max(2 * hole_len, 2)))))
    tdfilter = _get_tdfilter(invpsd, max_samples=req_max_samples)
    diag = float(tdfilter[0])
    raw_col = tdfilter[:hole_len]
    if len(raw_col) < hole_len:
        raw_col = np.pad(raw_col, (0, hole_len - len(raw_col)))
    mat = linalg.toeplitz((raw_col / diag).astype(np.float64))
    if ridge > 0:
        mat += ridge * np.eye(hole_len, dtype=np.float64)
    return linalg.inv(mat) / diag


def gate_and_paint_matmul(
    data,
    lindex,
    rindex,
    invpsd,
    invmat=None,
    ridge=1e-7,
    copy=True,
    pad_duration=16.0,
):
    """Gate and inpaint data using explicit matrix multiplication.

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
        The uninverted covariance matrix. If None, calculated from ``invpsd``.
    ridge : float, optional
        Relative diagonal Tikhonov regularization parameter when computing
        ``invmat``. Default is 1e-7.
    copy : bool, optional
        Copy the data before applying the gate. Otherwise, apply in-place.
        Default is True.
    pad_duration : float, optional
        One-sided padding duration in seconds for local convolution inpainting.
        Default is 16.0.

    Returns
    -------
    TimeSeries :
        The gated and in-painted time series.

    """
    if copy:
        data = data.copy()
    if rindex <= lindex:
        return data
    data[lindex:rindex] = 0.0
    _, rhs = _compute_inpaint_rhs(data, lindex, rindex, invpsd, pad_duration)

    if invmat is None:
        invmat = invert_covariance(invpsd, lindex, rindex, ridge=ridge)

    # Remove the projection into the null space
    proj = invmat.astype(np.float64) @ rhs
    data[lindex:rindex] -= proj.astype(data.dtype)
    return data


def gate_and_paint_multi(
    data,
    intervals,
    invpsd,
    copy=True,
    method="toeplitz",
    ridge=1e-7,
    pad_duration=16.0,
    fallback_taper=True,
    maxiter=64,
    tol=1e-4,
):
    """Gate and jointly inpaint multiple intervals using block-Toeplitz PCG.

    Clusters nearby gates within ``pad_duration`` so dense glitches cannot
    exponentially amplify across sequential single-gate solves, while isolated
    gates dispatch directly to :func:`gate_and_paint`.

    Parameters
    ----------
    data : TimeSeries
        The data to gate.
    intervals : sequence of tuple of (int, int)
        Sample index intervals ``(lindex, rindex)`` to gate and inpaint.
    invpsd : FrequencySeries
        The inverse of the PSD.
    copy : bool, optional
        Copy the data before applying the gates. Otherwise, apply in-place.
        Default is True.
    method : {'toeplitz', 'cholesky', 'auto'}, optional
        Solver used for single-hole gates and block-diagonal preconditioning.
        Default is ``'toeplitz'``.
    ridge : float, optional
        Relative diagonal Tikhonov regularization parameter. Default is 1e-7.
    pad_duration : float, optional
        One-sided padding duration in seconds for local convolution inpainting.
        Default is 16.0.
    fallback_taper : bool, optional
        If True (default), fall back to cosine taper gating for any merged
        interval exceeding 32768 samples.
    maxiter : int, optional
        Maximum preconditioned conjugate-gradient iterations per dense cluster.
        Default is 64.
    tol : float, optional
        Relative residual tolerance for preconditioned conjugate gradient.
        Default is 1e-4.

    Returns
    -------
    TimeSeries :
        The gated and jointly inpainted time series.

    """
    if copy:
        data = data.copy()

    # Sort and merge overlapping or touching intervals inside [0, len(data)]
    cleaned = sorted(
        (max(0, int(li)), min(len(data), int(ri)))
        for li, ri in intervals
        if min(len(data), int(ri)) > max(0, int(li))
    )
    if not cleaned:
        return data

    merged = [list(cleaned[0])]
    for li, ri in cleaned[1:]:
        if li <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], ri)
        else:
            merged.append([li, ri])

    # Apply ultra-wide cosine-taper fallback before clustering remaining holes
    paint_intervals = []
    for li, ri in merged:
        if fallback_taper and (ri - li) > 32768:
            gate_and_paint(
                data,
                li,
                ri,
                invpsd,
                copy=False,
                method=method,
                ridge=ridge,
                pad_duration=pad_duration,
                fallback_taper=True,
            )
        else:
            paint_intervals.append((li, ri))
    if not paint_intervals:
        return data

    if len(paint_intervals) == 1:
        li, ri = paint_intervals[0]
        return gate_and_paint(
            data,
            li,
            ri,
            invpsd,
            copy=False,
            method=method,
            ridge=ridge,
            pad_duration=pad_duration,
            fallback_taper=fallback_taper,
        )

    max_hole = max(ri - li for li, ri in paint_intervals)
    req_max = max(1048576, 1 << int(np.ceil(np.log2(max(2 * max_hole, 2)))))
    tdfilter = _get_tdfilter(invpsd, max_samples=req_max)
    diag = float(tdfilter[0])
    pad_samples = (
        min(len(tdfilter) // 2, round(pad_duration / float(data.delta_t)))
        if pad_duration is not None
        else len(data) // 2
    )

    # Cluster intervals separated by <= pad_samples so dense bursts solve jointly
    clusters = [[paint_intervals[0]]]
    for li, ri in paint_intervals[1:]:
        if pad_samples > 0 and (li - clusters[-1][-1][1]) <= pad_samples:
            clusters[-1].append((li, ri))
        elif pad_duration is None:
            clusters[0].append((li, ri))
        else:
            clusters.append([(li, ri)])

    invpsd_f64 = invpsd.astype(np.float64)

    for cluster in clusters:
        if len(cluster) == 1:
            li, ri = cluster[0]
            gate_and_paint(
                data,
                li,
                ri,
                invpsd,
                copy=False,
                method=method,
                ridge=ridge,
                pad_duration=pad_duration,
                fallback_taper=fallback_taper,
            )
            continue

        # Local padded slice containing the cluster, or full-array periodic FFT if short
        c_start = cluster[0][0]
        c_end = cluster[-1][1]
        use_local = pad_samples > 0 and len(data) > 2 * pad_samples + (c_end - c_start)
        if use_local:
            idx_start = c_start - pad_samples
            idx_end = c_end + pad_samples
            pad_left = max(0, -idx_start)
            pad_right = max(0, idx_end - len(data))
            idx_start = max(0, idx_start)
            idx_end = min(len(data), idx_end)
            offset = pad_left - idx_start
            buf = data[idx_start:idx_end].numpy().astype(np.float64)
            if pad_left > 0 or pad_right > 0:
                buf = np.pad(buf, (pad_left, pad_right), mode="constant")
            kernel = np.concatenate(
                [tdfilter[-pad_samples:], tdfilter[: pad_samples + 1]]
            ).astype(np.float64)

            def apply_cinv(vec_full, kern=kernel):
                return signal.fftconvolve(vec_full, kern, mode="same")
        else:
            offset = 0
            buf = data.numpy().astype(np.float64)

            def apply_cinv(vec_full, inv_f=invpsd_f64, dt=float(data.delta_t)):
                ts = TimeSeries(vec_full, delta_t=dt, copy=False)
                return (
                    (ts.to_frequencyseries() * inv_f)
                    .to_timeseries()
                    .numpy()
                    .astype(np.float64)
                )

        local_slices = [(li + offset, ri + offset) for li, ri in cluster]
        mask = np.zeros(len(buf), dtype=bool)
        for li, ri in local_slices:
            buf[li:ri] = 0.0
            mask[li:ri] = True

        # Build block-diagonal Toeplitz preconditioner and RHS from zeroed data
        precond_cols = {}
        for li, ri in local_slices:
            klen = ri - li
            if klen not in precond_cols:
                col = tdfilter[:klen].astype(np.float64).copy()
                if ridge > 0:
                    col[0] += ridge * diag
                precond_cols[klen] = col

        b_vec = apply_cinv(buf)[mask]
        b_norm = float(np.linalg.norm(b_vec))
        if b_norm == 0.0:
            for li, ri in cluster:
                data[li:ri] = 0.0
            continue

        def apply_precond(vec_h, slices=local_slices, cols=precond_cols):
            out = np.empty_like(vec_h)
            pos = 0
            for li, ri in slices:
                klen = ri - li
                out[pos : pos + klen] = linalg.solve_toeplitz(
                    cols[klen], vec_h[pos : pos + klen]
                )
                pos += klen
            return out

        # Solve SPD (C^-1_HH + ridge * diag * I) x_H = b_H with PCG
        x_vec = np.zeros_like(b_vec)
        r_vec = b_vec.copy()
        z_vec = apply_precond(r_vec)
        p_vec = z_vec.copy()
        rz = float(np.dot(r_vec, z_vec))
        for _ in range(maxiter):
            v_full = np.zeros(len(buf), dtype=np.float64)
            v_full[mask] = p_vec
            ap_vec = apply_cinv(v_full)[mask]
            if ridge > 0:
                ap_vec = ap_vec + (ridge * diag) * p_vec
            denom = float(np.dot(p_vec, ap_vec))
            if denom <= 0.0:
                break
            alpha = rz / denom
            x_vec += alpha * p_vec
            r_vec -= alpha * ap_vec
            if float(np.linalg.norm(r_vec)) / b_norm < tol:
                break
            z_vec = apply_precond(r_vec)
            rz_new = float(np.dot(r_vec, z_vec))
            if rz == 0.0:
                break
            p_vec = z_vec + (rz_new / rz) * p_vec
            rz = rz_new

        pos = 0
        for li, ri in cluster:
            klen = ri - li
            data[li:ri] = TimeSeries(
                (-x_vec[pos : pos + klen]).astype(data.dtype),
                delta_t=data.delta_t,
                copy=False,
            )
            pos += klen

    return data
