"""
This module contains functions for calculating single-ifo ranking
statistic values
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def effsnr(snr, reduced_x2, fac=250.0, **kwargs):  # pylint:disable=unused-argument
    """Calculate the effective SNR statistic. See (S5y1 paper) for definition."""
    snr = np.array(snr, ndmin=1, dtype=np.float64)
    rchisq = np.array(reduced_x2, ndmin=1, dtype=np.float64)
    esnr = snr / (1 + snr**2 / fac) ** 0.25 / rchisq**0.25

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return esnr
    return esnr[0]


def newsnr(
    snr,
    reduced_x2,
    q=6.0,
    n=2.0,
    chisq_thresh=1.0,
    chisq_scale=0.0,
    chisq_max_snr2=300.0,
    **kwargs,  # pylint:disable=unused-argument
):
    """Calculate the re-weighted SNR statistic ('newSNR').

    Parameters
    ----------
    snr : float or array_like
        Matched-filter signal-to-noise ratio values.
    reduced_x2 : float or array_like
        Reduced power chi-squared values.
    q : float, optional
        Exponent controlling asymptotic penalty steepness. Default is 6.0.
    n : float, optional
        Chi-squared normalization power. Default is 2.0.
    chisq_thresh : float, optional
        Baseline reduced chi-squared plateau above which reweighting begins.
        Default is 1.0.
    chisq_scale : float, optional
        Quadratic SNR-squared tolerance coefficient added to ``chisq_thresh``
        as ``chisq_scale * min(snr**2, chisq_max_snr2)``. Default is 0.0.
    chisq_max_snr2 : float, optional
        Upper cap on ``snr**2`` for quadratic tolerance scaling. Default is
        300.0.

    Returns
    -------
    float or numpy.ndarray
        Re-weighted SNR values matching scalar/array type of ``snr``.
    """
    snr_arr = np.array(snr, ndmin=1, dtype=np.float64)
    nsnr = snr_arr.copy()
    reduced_x2 = np.array(reduced_x2, ndmin=1, dtype=np.float64)

    # Apply SNR^2-scaled threshold tolerance above chisq_thresh
    thresh = chisq_thresh + chisq_scale * np.minimum(snr_arr**2, chisq_max_snr2)
    ind = np.where(reduced_x2 > thresh)[0]
    if len(ind):
        scaled_x2 = reduced_x2[ind] / thresh[ind]
        nsnr[ind] *= (0.5 * (1.0 + scaled_x2 ** (q / n))) ** (-1.0 / q)

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto(
    snr,
    brchisq,
    sgchisq,
    sgchisq_thresh=4.0,
    sgchisq_scale=0.0,
    sgchisq_max_snr2=5000.0,
    **kwargs,
):
    """Calculate combined SNR derived from NewSNR and sine-Gaussian chisq.

    Parameters
    ----------
    snr : float or array_like
        Matched-filter signal-to-noise ratio values.
    brchisq : float or array_like
        Reduced power chi-squared values.
    sgchisq : float or array_like
        Sine-Gaussian chi-squared values.
    sgchisq_thresh : float, optional
        Baseline sine-Gaussian chi-squared threshold. Default is 4.0.
    sgchisq_scale : float, optional
        Quadratic higher-order-mode tolerance coefficient added to
        ``sgchisq_thresh`` as ``sgchisq_scale * min(snr**2, sgchisq_max_snr2)``.
        Default is 0.0.
    sgchisq_max_snr2 : float, optional
        Upper cap on ``snr**2`` for sine-Gaussian tolerance scaling. Default is
        5000.0.

    Returns
    -------
    float or numpy.ndarray
        Re-weighted SNR values matching scalar/array type of ``snr``.
    """
    snr_arr = np.array(snr, ndmin=1, dtype=np.float64)
    nsnr = np.array(newsnr(snr_arr, brchisq, **kwargs), ndmin=1)
    sgchisq = np.array(sgchisq, ndmin=1, dtype=np.float64)

    # Downweight triggers exceeding the SNR^2-scaled sine-Gaussian threshold
    thresh = sgchisq_thresh + sgchisq_scale * np.minimum(snr_arr**2, sgchisq_max_snr2)
    mask = sgchisq > thresh
    if np.any(mask):
        nsnr[mask] /= (sgchisq[mask] / thresh[mask]) ** 0.5

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, "__len__"):
        return nsnr
    return nsnr[0]


def newsnr_sgveto_autochisq(
    snr,
    brchisq,
    sgchisq,
    autochisq,
    q=6.0,
    n=2.0,
    autochisq_thresh=1.0,
    autochisq_scale=0.0,
    autochisq_max_snr2=300.0,
    **kwargs,
):
    """Calculate combined SNR from NewSNR, sine-Gaussian chisq, and auto-chisq.

    Parameters
    ----------
    snr : float or array_like
        Matched-filter signal-to-noise ratio values.
    brchisq : float or array_like
        Reduced power chi-squared values.
    sgchisq : float or array_like
        Sine-Gaussian chi-squared values.
    autochisq : float or array_like
        Reduced autocorrelation chi-squared values.
    q : float, optional
        Exponent controlling asymptotic penalty steepness. Default is 6.0.
    n : float, optional
        Chi-squared normalization power. Default is 2.0.
    autochisq_thresh : float, optional
        Baseline reduced auto-chisq threshold. Default is 1.0.
    autochisq_scale : float, optional
        Quadratic SNR-squared tolerance coefficient for auto-chisq. Default is
        0.0.
    autochisq_max_snr2 : float, optional
        Upper cap on ``snr**2`` for auto-chisq scaling. Default is 300.0.

    Returns
    -------
    float or numpy.ndarray
        Re-weighted SNR values matching scalar/array type of ``snr``.
    """
    snr_arr = np.array(snr, ndmin=1, dtype=np.float64)
    nsnr = np.array(
        newsnr_sgveto(snr_arr, brchisq, sgchisq, q=q, n=n, **kwargs),
        ndmin=1,
    )
    rauto = np.array(autochisq, ndmin=1, dtype=np.float64)

    # Apply continuous NewSNR-style reweighting above auto-chisq threshold
    thresh = autochisq_thresh + autochisq_scale * np.minimum(
        snr_arr**2, autochisq_max_snr2
    )
    mask = rauto > thresh
    if np.any(mask):
        scaled = rauto[mask] / thresh[mask]
        nsnr[mask] *= (0.5 * (1.0 + scaled ** (q / n))) ** (-1.0 / q)

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
    psd_var_val = np.array(psd_var_val, copy=True)
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
    nsnr = np.array(nsnr, ndmin=1)
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
    nsnr = np.array(newsnr_sgveto(snr, brchisq, sgchisq, **kwargs), ndmin=1)
    psd_var_val = np.array(psd_var_val, ndmin=1, copy=True)
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
    nsnr = np.array(nsnr, ndmin=1)
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
    return np.array(trigs["snr"][:], ndmin=1, dtype=np.float32)


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
    return np.array(nsnr, ndmin=1, dtype=np.float32)


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
    return np.array(nsnr_sg, ndmin=1, dtype=np.float32)


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
    return np.array(nsnr_sg_psd, ndmin=1, dtype=np.float32)


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
    return np.array(nsnr_sg_psdt, ndmin=1, dtype=np.float32)


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
    return np.array(nsnr_sg_psdscale, ndmin=1, dtype=np.float32)


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
    return np.array(nsnr_sg_psdt, ndmin=1, dtype=np.float32)


def get_newsnr_sgveto_autochisq(trigs, **kwargs):
    """Calculate newsnr re-weighted by sine-Gaussian and auto-chisq vetoes.

    Parameters
    ----------
    trigs : dict of numpy.ndarray or h5py.Group
        Mapping holding single-detector trigger arrays, including ``'snr'``,
        ``'chisq'``, ``'chisq_dof'``, ``'sg_chisq'``, ``'cont_chisq'``, and
        ``'cont_chisq_dof'``.

    Returns
    -------
    numpy.ndarray
        Array of re-weighted SNR values.
    """
    dof = 2.0 * trigs["chisq_dof"][:] - 2.0
    red_brchisq = trigs["chisq"][:] / dof
    auto_dof = trigs["cont_chisq_dof"][:]
    red_autochisq = np.where(
        auto_dof > 0, trigs["cont_chisq"][:] / np.maximum(auto_dof, 1), 1.0
    )
    nsnr_sga = newsnr_sgveto_autochisq(
        trigs["snr"][:],
        red_brchisq,
        trigs["sg_chisq"][:],
        red_autochisq,
        **kwargs,
    )
    return np.array(nsnr_sga, ndmin=1, dtype=np.float32)


sngls_ranking_function_dict = {
    "snr": get_snr,
    "newsnr": get_newsnr,
    "new_snr": get_newsnr,
    "newsnr_sgveto": get_newsnr_sgveto,
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
reqd_datasets["newsnr_sgveto"] = reqd_datasets["newsnr"] + ["sg_chisq"]
reqd_datasets["newsnr_sgveto_autochisq"] = reqd_datasets["newsnr_sgveto"] + [
    "cont_chisq",
    "cont_chisq_dof",
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
