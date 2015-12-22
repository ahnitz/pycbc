from numpy import *
from scipy import special
from scipy.integrate import quad
import fstat


def loglike_approx(a_hat, f_plus, f_cross, d_max=1000., method="coh"):
    """
    Calculate the approximate likelihood. This works for three cases:
    left and right circularly polarized and the standard coherent analysis.
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param d_max: maximum distance for marginalization
    :param method: approximation for calculating likelihood, one of "coh", "left", "right"
    """
    if (method == "left") or (method == "right"):
        a_hat = fstat.circ_project(a_hat, f_plus, f_cross, method)
    d_hat, cosi_hat, _, _ = fstat.a_to_params(a_hat)
    snr = fstat.expected_snr(a_hat, f_plus, f_cross)

    if snr == 0:
        loglike = 0
    elif method == "coh":
        loglike = log(32 * (d_hat / d_max) ** 3 * d_hat ** 4 /
                   (f_plus ** 2 * f_cross ** 2) / (1 - cosi_hat ** 2) ** 3)
    else:
        # the width in cos iota:
        cos_fac = sqrt((f_cross ** 2 + f_plus ** 2) / (f_plus * f_cross))
        cos_width = minimum(cos_fac / snr ** 0.5, 0.5)
        loglike = log((d_hat / d_max) ** 3 / snr ** 2 * cos_width)

    return loglike, snr


def like_approx(a_hat, f_plus, f_cross, d_max=1000.):
    """
    Calculate the approximate likelihood summed over left, right and coherent.
    """
    loglike = {}
    snr = {}
    like = {}
    for method in ["left", "right", "coh"]:
        loglike[method], snr[method] = loglike_approx(a_hat, f_plus, f_cross, d_max, method= method)
        like[method] = snr[method]**2 / 2 + loglike[method]
        if snr[method] < 6:
            like[method] = 0

    if ((snr["coh"] ** 2 - snr["right"] ** 2) < 1) or ((snr["coh"] ** 2 - snr["left"] ** 2) < 1):
        like["coh"] = 0

    like_approx = logaddexp(logaddexp(like["left"], like["right"]), like["coh"])


    return like_approx, like
