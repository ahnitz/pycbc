"""
This module contains functions for calculating single-ifo ranking
statistic values
"""

import logging

import numpy

logger = logging.getLogger("pycbc.events.ranking")


def effsnr(snr, reduced_x2, fac=250.0, **kwargs):  # pylint:disable=unused-argument
    """Calculate the effective SNR statistic. See (S5y1 paper) for definition."""
    snr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    rchisq = numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64)
    esnr = snr / (1 + snr**2 / fac) ** 0.25 / rchisq**0.25

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return esnr
    return esnr[0]


def newsnr(snr, reduced_x2, q=6.0, n=2.0, chisq_thresh=1.0, **kwargs):  # pylint:disable=unused-argument
    """
    Calculate the re-weighted SNR statistic ('newSNR') from given SNR and
    reduced chi-squared values. See http://arxiv.org/abs/1208.3491 for
    definition. Previous implementation in glue/ligolw/lsctables.py
    """
    nsnr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    reduced_x2 = numpy.abs(numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64))

    # newsnr is only different from snr if reduced chisq > chisq_thresh
    ind = numpy.where(reduced_x2 > chisq_thresh)[0]
    scaled_x2 = reduced_x2[ind] / chisq_thresh
    nsnr[ind] *= (0.5 * (1.0 + scaled_x2 ** (q / n))) ** (-1.0 / q)

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_offset125(snr, reduced_x2, **kwargs):
    """Calculate NewSNR with chisq_thresh=1.25"""
    return newsnr(snr, reduced_x2, chisq_thresh=1.25, **kwargs)


def newsnr_sgveto(snr, brchisq, sgchisq, **kwargs):
    """Combined SNR derived from NewSNR and Sine-Gaussian Chisq"""
    nsnr = numpy.array(newsnr(snr, brchisq, **kwargs), ndmin=1)
    sgchisq = numpy.array(sgchisq, ndmin=1)
    t = numpy.array(sgchisq > 4, ndmin=1)
    if len(t):
        nsnr[t] = nsnr[t] / (sgchisq[t] / 4.0) ** 0.5

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto_offset125(snr, brchisq, sgchisq, **kwargs):
    """Combined SNR derived from NewSNR with chisq_thresh=1.25 and Sine-Gaussian Chisq"""
    return newsnr_sgveto(snr, brchisq, sgchisq, chisq_thresh=1.25, **kwargs)


def newsnr_sgveto_autochisq(
    snr,
    brchisq,
    sgchisq,
    autochisq,
    q=6.0,
    n=1.0,
    c_sg=0.012,
    hollow=None,
    hollow_knee=0.0,
    hollow_scale=2.5,
    **kwargs,
):
    """Multi-metric signal-consistency distance from the BBH injection envelope
    in (snr, reduced_power_chisq, reduced_auto_chisq, sg_chisq, hollow) space.
    Signals inside the 95% injection envelope experience zero downweighting;
    triggers outside the envelope are smoothly downweighted by their combined
    normalized distance D = sqrt(z_pow^2 + z_auto^2 + z_sg^2 + z_hollow^2)."""
    nsnr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    s2_pa = numpy.minimum(nsnr**2, 300.0)
    s2_sg = numpy.minimum(nsnr**2, 5000.0)

    d2 = numpy.zeros_like(nsnr)

    if brchisq is not None:
        rpow = numpy.abs(numpy.array(brchisq, ndmin=1, dtype=numpy.float64))
        tau_p = 1.50 + 0.0018 * s2_pa
        sig_p = 0.45 + 0.0008 * s2_pa
        zp = numpy.maximum(0.0, (rpow - tau_p) / sig_p)
        d2 += zp**2

    if autochisq is not None:
        rauto = numpy.abs(numpy.array(autochisq, ndmin=1, dtype=numpy.float64))
        tau_a = 1.85 + 0.0028 * s2_pa
        sig_a = 0.65 + 0.0012 * s2_pa
        za = numpy.maximum(0.0, (rauto - tau_a) / sig_a)
        d2 += za**2

    if sgchisq is not None:
        rsg = numpy.abs(numpy.array(sgchisq, ndmin=1, dtype=numpy.float64))
        tau_s = 3.20 + c_sg * s2_sg
        sig_s = 1.20 + 0.0035 * s2_sg
        zs = numpy.maximum(0.0, (rsg - tau_s) / sig_s)
        d2 += zs**2

    if hollow is not None and hollow_knee > 0.0:
        rhollow = numpy.abs(numpy.array(hollow, ndmin=1, dtype=numpy.float64))
        zh = numpy.maximum(0.0, (rhollow - hollow_knee) / hollow_scale)
        d2 += zh**2

    dist = numpy.sqrt(d2)
    ind = numpy.where(dist > 0.0)[0]
    if len(ind):
        nsnr[ind] *= (0.5 * (1.0 + (1.0 + dist[ind]) ** (q / n))) ** (-1.0 / q)

    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]



def newsnr_sgveto_psdvar(
    snr, brchisq, sgchisq, psd_var_val, min_expected_psdvar=0.65, **kwargs
):
    """
    Combined SNR derived from SNR, reduced Allen chisq, sine-Gaussian chisq and
    PSD variation statistic
    """
    # If PSD var is lower than the 'minimum usually expected value' stop this
    # being used in the statistic. This low value might arise because a
    # significant fraction of the "short" PSD period was gated (for instance).
    psd_var_val = numpy.array(psd_var_val, copy=True)
    psd_var_val[psd_var_val < min_expected_psdvar] = 1.0
    scaled_snr = snr * (psd_var_val**-0.5)
    scaled_brchisq = brchisq * (psd_var_val**-1.0)
    nsnr = newsnr_sgveto(scaled_snr, scaled_brchisq, sgchisq, **kwargs)

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto_psdvar_threshold(
    snr,
    brchisq,
    sgchisq,
    psd_var_val,
    min_expected_psdvar=0.65,
    brchisq_threshold=10.0,
    psd_var_val_threshold=10.0,
    **kwargs,
):
    """
    newsnr_sgveto_psdvar with thresholds applied.

    This is the newsnr_sgveto_psdvar statistic with additional options
    to threshold on chi-squared or PSD variation.
    """
    nsnr = newsnr_sgveto_psdvar(
        snr,
        brchisq,
        sgchisq,
        psd_var_val,
        min_expected_psdvar=min_expected_psdvar,
        **kwargs,
    )
    nsnr = numpy.array(nsnr, ndmin=1)
    nsnr[brchisq > brchisq_threshold] = 1.0
    nsnr[psd_var_val > psd_var_val_threshold] = 1.0

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto_psdvar_scaled(
    snr, brchisq, sgchisq, psd_var_val, scaling=0.33, min_expected_psdvar=0.65, **kwargs
):
    """
    Combined SNR derived from NewSNR, Sine-Gaussian Chisq and scaled PSD
    variation statistic.
    """
    nsnr = numpy.array(newsnr_sgveto(snr, brchisq, sgchisq, **kwargs), ndmin=1)
    psd_var_val = numpy.array(psd_var_val, ndmin=1, copy=True)
    psd_var_val[psd_var_val < min_expected_psdvar] = 1.0

    # Default scale is 0.33 as tuned from analysis of data from O2 chunks
    nsnr = nsnr / psd_var_val**scaling

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto_psdvar_scaled_threshold(
    snr, bchisq, sgchisq, psd_var_val, threshold=2.0, **kwargs
):
    """
    Combined SNR derived from NewSNR and Sine-Gaussian Chisq, and
    scaled psd variation.
    """
    nsnr = newsnr_sgveto_psdvar_scaled(snr, bchisq, sgchisq, psd_var_val, **kwargs)
    nsnr = numpy.array(nsnr, ndmin=1)
    nsnr[bchisq > threshold] = 1.0

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def get_snr(trigs, **kwargs):  # pylint:disable=unused-argument
    """
    Return SNR from a trigs/dictionary object

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
        Dictionary-like object holding single detector trigger information.
        'snr' is a required key

    Returns
    -------
    numpy.ndarray
        Array of snr values

    """
    return numpy.array(trigs["snr"][:], ndmin=1, dtype=numpy.float32)


def get_newsnr(trigs, **kwargs):
    """
    Calculate newsnr ('reweighted SNR') for a trigs/dictionary object

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
        Dictionary-like object holding single detector trigger information.
        'chisq_dof', 'snr', and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr = newsnr(trigs["snr"][:], trigs["chisq"][:] / dof, **kwargs)
    return numpy.array(nsnr, ndmin=1, dtype=numpy.float32)


def get_newsnr_offset125(trigs, **kwargs):
    """Calculate newsnr with chisq_thresh=1.25"""
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr = newsnr_offset125(trigs["snr"][:], trigs["chisq"][:] / dof, **kwargs)
    return numpy.array(nsnr, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto(trigs, **kwargs):
    """
    Calculate newsnr re-weigthed by the sine-gaussian veto

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
        Dictionary-like object holding single detector trigger information.
        'chisq_dof', 'snr', 'sg_chisq' and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg = newsnr_sgveto(
        trigs["snr"][:], trigs["chisq"][:] / dof, trigs["sg_chisq"][:], **kwargs
    )
    return numpy.array(nsnr_sg, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_offset125(trigs, **kwargs):
    """Calculate newsnr_sgveto with chisq_thresh=1.25"""
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg = newsnr_sgveto_offset125(
        trigs["snr"][:], trigs["chisq"][:] / dof, trigs["sg_chisq"][:], **kwargs
    )
    return numpy.array(nsnr_sg, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_psdvar(trigs, **kwargs):
    """
    Calculate snr re-weighted by Allen chisq, sine-gaussian veto and
    psd variation statistic

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', 'chisq' and 'psd_var_val' are required keys

    Returns
    -------
     numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg_psd = newsnr_sgveto_psdvar(
        trigs["snr"][:],
        trigs["chisq"][:] / dof,
        trigs["sg_chisq"][:],
        trigs["psd_var_val"][:],
        **kwargs,
    )
    return numpy.array(nsnr_sg_psd, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_psdvar_threshold(trigs, **kwargs):
    """
    Calculate newsnr re-weighted by the sine-gaussian veto and scaled
    psd variation statistic

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', 'chisq' and 'psd_var_val' are required keys

    Returns
    -------
     numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg_psdt = newsnr_sgveto_psdvar_threshold(
        trigs["snr"][:],
        trigs["chisq"][:] / dof,
        trigs["sg_chisq"][:],
        trigs["psd_var_val"][:],
        **kwargs,
    )
    return numpy.array(nsnr_sg_psdt, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_psdvar_scaled(trigs, **kwargs):
    """
    Calculate newsnr re-weighted by the sine-gaussian veto and scaled
    psd variation statistic

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', 'chisq' and 'psd_var_val' are required keys

    Returns
    -------
     numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg_psdscale = newsnr_sgveto_psdvar_scaled(
        trigs["snr"][:],
        trigs["chisq"][:] / dof,
        trigs["sg_chisq"][:],
        trigs["psd_var_val"][:],
        **kwargs,
    )
    return numpy.array(nsnr_sg_psdscale, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_psdvar_scaled_threshold(trigs, **kwargs):
    """
    Calculate newsnr re-weighted by the sine-gaussian veto and scaled
    psd variation statistic. A further threshold is applied to the
    reduced chisq.

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', 'chisq' and 'psd_var_val' are required keys

    Returns
    -------
     numpy.ndarray
        Array of newsnr values

    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    nsnr_sg_psdt = newsnr_sgveto_psdvar_scaled_threshold(
        trigs["snr"][:],
        trigs["chisq"][:] / dof,
        trigs["sg_chisq"][:],
        trigs["psd_var_val"][:],
        **kwargs,
    )
    return numpy.array(nsnr_sg_psdt, ndmin=1, dtype=numpy.float32)


def _extract_clean_hollow_times(
    trigs,
    min_snr=6.0,
    max_rchisq=2.5,
):
    """Extract and 0.1s-cluster clean triggers from trigs or its underlying dataset."""
    source = None
    if isinstance(trigs, dict) and "_full_trigs" in trigs:
        source = trigs["_full_trigs"]
    elif hasattr(trigs, "trigs") and hasattr(trigs.trigs, "__getitem__"):
        # SingleDetTriggers
        source = trigs.trigs
    elif (
        hasattr(trigs, "file")
        and hasattr(trigs, "ifo")
        and hasattr(trigs.file, "__getitem__")
        and trigs.ifo in trigs.file
    ):
        # ReadByTemplate
        source = trigs.file[trigs.ifo]
    elif hasattr(trigs, "__getitem__") and "snr" in trigs and "chisq" in trigs:
        # HDF5 Group or trigger dict
        source = trigs

    if source is None:
        return None

    cache_key = ("_hollow_times", float(min_snr), float(max_rchisq))
    if hasattr(source, "_hollow_cache") and cache_key in source._hollow_cache:
        return source._hollow_cache[cache_key]

    try:
        all_snr = source["snr"][:]
        all_chisq = source["chisq"][:]
        all_dof = 2.0 * source["chisq_dof"][:] - 2.0
        all_rchisq = all_chisq / numpy.maximum(1.0, all_dof)
        clean = (all_snr >= min_snr) & (all_rchisq <= max_rchisq)

        if "end_time" in source:
            clean_times = numpy.array(source["end_time"][:][clean], dtype=numpy.float64)
            if "end_time_ns" in source:
                clean_times += (
                    numpy.array(source["end_time_ns"][:][clean], dtype=numpy.float64)
                    * 1e-9
                )
        elif "time" in source:
            clean_times = numpy.array(source["time"][:][clean], dtype=numpy.float64)
        else:
            return None

        if len(clean_times) == 0:
            clustered = numpy.empty(0, dtype=numpy.float64)
        else:
            bin_idx = numpy.floor(clean_times / 0.1).astype(numpy.int64)
            unique_bins = numpy.unique(bin_idx)
            clustered = unique_bins * 0.1 + 0.05

        try:
            if not hasattr(source, "_hollow_cache"):
                source._hollow_cache = {}
            source._hollow_cache[cache_key] = clustered
        except (AttributeError, TypeError):
            pass

        return clustered
    except (KeyError, ValueError, TypeError):
        return None


def compute_hollow_counts(
    t,
    clustered_clean_times,
    inner=0.25,
    outer=12.0,
):
    """Compute count of clustered clean triggers in hollow annulus [inner, outer]."""
    if clustered_clean_times is None or len(clustered_clean_times) == 0:
        return numpy.zeros_like(t, dtype=numpy.float64)

    t_arr = numpy.array(t, ndmin=1, dtype=numpy.float64)
    left_start = numpy.searchsorted(clustered_clean_times, t_arr - outer, side="left")
    left_end = numpy.searchsorted(clustered_clean_times, t_arr - inner, side="right")
    count_left = left_end - left_start

    right_start = numpy.searchsorted(clustered_clean_times, t_arr + inner, side="left")
    right_end = numpy.searchsorted(clustered_clean_times, t_arr + outer, side="right")
    count_right = right_end - right_start

    H = (count_left + count_right).astype(numpy.float64)
    if hasattr(t, "__len__"):
        return H
    return H[0]


def get_newsnr_sgveto_autochisq(
    trigs,
    chisq_thresh=1.50,
    sg_thresh=4.0,
    auto_thresh=4.0,
    n_trad=2.0,
    n_sg=2.0,
    n_auto=2.0,
    c_sg=0.012,
    q=6.0,
    **kwargs,
):
    """Calculate newsnr re-weighted by traditional chisq, sine-gaussian chisq,
    auto-correlation chisq, and hollow-annulus non-stationarity"""
    if "sngl_ranking_c_sg" in kwargs:
        c_sg = float(kwargs.pop("sngl_ranking_c_sg"))
    dof = numpy.maximum(1.0, 2.0 * trigs["chisq_dof"][:] - 2.0)
    rchisq = trigs["chisq"][:] / dof
    try:
        sgchisq = trigs["sg_chisq"][:]
    except (KeyError, ValueError):
        sgchisq = None

    if hasattr(trigs, "auto_rchisq") and trigs.auto_rchisq is not None:
        autochisq = trigs.auto_rchisq
    else:
        try:
            autochisq = trigs["cont_chisq"][:]
            try:
                adof = numpy.maximum(1.0, trigs["cont_chisq_dof"][:])
                autochisq = autochisq / adof
            except (KeyError, ValueError):
                pass
        except (KeyError, ValueError):
            try:
                autochisq = trigs["auto_chisq"][:]
                try:
                    adof = numpy.maximum(1.0, trigs["auto_chisq_dof"][:])
                    autochisq = autochisq / adof
                except (KeyError, ValueError):
                    pass
            except (KeyError, ValueError):
                autochisq = None

    # Parse hollow-annulus non-stationarity parameters
    hollow_knee = float(
        kwargs.pop("sngl_ranking_hollow_knee", kwargs.pop("hollow_knee", 0.0))
    )
    hollow_scale = float(
        kwargs.pop("sngl_ranking_hollow_scale", kwargs.pop("hollow_scale", 2.5))
    )
    hollow_inner = float(
        kwargs.pop("sngl_ranking_hollow_inner", kwargs.pop("hollow_inner", 0.25))
    )
    hollow_outer = float(
        kwargs.pop("sngl_ranking_hollow_outer", kwargs.pop("hollow_outer", 12.0))
    )
    hollow_min_snr = float(
        kwargs.pop("sngl_ranking_hollow_min_snr", kwargs.pop("hollow_min_snr", 6.0))
    )
    hollow_max_rchisq = float(
        kwargs.pop(
            "sngl_ranking_hollow_max_rchisq", kwargs.pop("hollow_max_rchisq", 2.5)
        )
    )

    hollow_counts = None
    if hollow_knee > 0.0:
        t = None
        if "end_time" in trigs:
            t = numpy.array(trigs["end_time"][:], dtype=numpy.float64)
            if "end_time_ns" in trigs:
                t += (
                    numpy.array(trigs["end_time_ns"][:], dtype=numpy.float64)
                    * 1e-9
                )
        elif "time" in trigs:
            t = numpy.array(trigs["time"][:], dtype=numpy.float64)

        if t is not None:
            clustered = _extract_clean_hollow_times(
                trigs,
                min_snr=hollow_min_snr,
                max_rchisq=hollow_max_rchisq,
            )
            hollow_counts = compute_hollow_counts(
                t,
                clustered,
                inner=hollow_inner,
                outer=hollow_outer,
            )

    nsnr = newsnr_sgveto_autochisq(
        trigs["snr"][:],
        rchisq,
        sgchisq,
        autochisq,
        chisq_thresh=chisq_thresh,
        sg_thresh=sg_thresh,
        auto_thresh=auto_thresh,
        n_trad=n_trad,
        n_sg=n_sg,
        n_auto=n_auto,
        c_sg=c_sg,
        q=q,
        hollow=hollow_counts,
        hollow_knee=hollow_knee,
        hollow_scale=hollow_scale,
        **kwargs,
    )
    return numpy.array(nsnr, ndmin=1, dtype=numpy.float32)


sngls_ranking_function_dict = {
    "snr": get_snr,
    "newsnr": get_newsnr,
    "new_snr": get_newsnr,
    "newsnr_offset125": get_newsnr_offset125,
    "newsnr_sgveto": get_newsnr_sgveto,
    "newsnr_sgveto_offset125": get_newsnr_sgveto_offset125,
    "newsnr_sgveto_autochisq": get_newsnr_sgveto_autochisq,
    "newsnr_sgveto_psdvar": get_newsnr_sgveto_psdvar,
    "newsnr_sgveto_psdvar_threshold": get_newsnr_sgveto_psdvar_threshold,
    "newsnr_sgveto_psdvar_scaled": get_newsnr_sgveto_psdvar_scaled,
    "newsnr_sgveto_psdvar_scaled_threshold": get_newsnr_sgveto_psdvar_scaled_threshold,
}

# Lists of datasets required in the trigs object for each function
reqd_datasets = {}
reqd_datasets["snr"] = ["snr"]
reqd_datasets["newsnr"] = reqd_datasets["snr"] + ["chisq", "chisq_dof"]
reqd_datasets["new_snr"] = reqd_datasets["newsnr"]
reqd_datasets["newsnr_offset125"] = reqd_datasets["newsnr"]
reqd_datasets["newsnr_sgveto"] = reqd_datasets["newsnr"] + ["sg_chisq"]
reqd_datasets["newsnr_sgveto_offset125"] = reqd_datasets["newsnr_sgveto"]
reqd_datasets["newsnr_sgveto_autochisq"] = reqd_datasets["newsnr_sgveto"] + [
    "cont_chisq",
    "cont_chisq_dof",
    "end_time",
]
reqd_datasets["newsnr_sgveto_psdvar"] = reqd_datasets["newsnr_sgveto"] + ["psd_var_val"]
reqd_datasets["newsnr_sgveto_psdvar_threshold"] = reqd_datasets["newsnr_sgveto_psdvar"]
reqd_datasets["newsnr_sgveto_psdvar_scaled"] = reqd_datasets["newsnr_sgveto_psdvar"]
reqd_datasets["newsnr_sgveto_psdvar_scaled_threshold"] = reqd_datasets[
    "newsnr_sgveto_psdvar"
]


def get_sngls_ranking_from_trigs(trigs, statname, **kwargs):
    """
    Return ranking for all trigs given a statname.

    Compute the single-detector ranking for a list of input triggers for a
    specific statname.

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, SingleDetTriggers or ReadByTemplate
        Dictionary holding single detector trigger information.
    statname:
        The statistic to use.

    """
    # Identify correct function
    try:
        sngl_func = sngls_ranking_function_dict[statname]
    except KeyError as exc:
        err_msg = f"Single-detector ranking {statname} not recognized"
        raise ValueError(err_msg) from exc

    # NOTE: In the sngl_funcs all the kwargs are explicitly stated, so any
    #       kwargs sent here must be known to the function.
    return sngl_func(trigs, **kwargs)
