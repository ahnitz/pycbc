from numpy import *
from scipy import special
from scipy.integrate import quad
import pylab
import fstat



def like_equal_d_cosi(a_hat, f, d, cosi):
    """
    For a network equally sensitive to plus and cross, calculate the
    likelihood marginalized over the 2 phases.  This a function of d and cosi.
    Note: we use uniform prior on phi, psi (=1/2pi)
    :param a_hat: the f-stat parameters of the signal
    :param f: the detector response (f = F+ = Fx)
    :param d: distance of template
    :param cosi: cos(inclination) of template
    """
    # calculate the two dimensional likelihood in circular polarization
    # marginalized over the two phases
    ar_hat, al_hat = fstat.a_to_circ_amp(a_hat)
    d0 = a_hat[0]
    al = (d0 / d) * (1 - cosi) ** 2 / 4
    ar = (d0 / d) * (1 + cosi) ** 2 / 4
    like = exp(- f ** 2 * (al - al_hat) ** 2) * \
           exp(- f ** 2 * (ar - ar_hat) ** 2) * \
           special.i0e(2 * f ** 2 * al * al_hat) * special.i0e(2 * f ** 2 * ar * ar_hat)
    return like


def like_equal_cosi(a_hat, f, x, d_max=1000., make_plot=False):
    """
    For a network equally sensitive to plus and cross, calculate the
    likelihood marginalized over the distance and 2 phases.
    This a function of cosi.   We use uniform priors, specifically 1/2pi for
    psi,phi, 1/2 for cosi, 3 d^2dd/d_max^3 for dist.
    :param a_hat: the f-stat parameters of the signal
    :param f: the detector response (f = F+ = Fx)
    :param x: cos(inclination) of template
    :param d_max: maximum distance for marginalization
    :param make_plot: plot the likelihood vs distance
    """
    ar_hat, al_hat = fstat.a_to_circ_amp(a_hat)
    d0 = a_hat[0]
    # we want to choose sensible ranges of integration:
    x4 = (1 + 6 * x ** 2 + x ** 4)
    a_peak = 2 * (al_hat * (1 - x) ** 2 + ar_hat * (1 + x) ** 2) / x4
    a_width = 2. * sqrt(2) / (f * sqrt(x4))
    l_circ = quad(lambda a: (3 * d0 ** 3) / d_max ** 3 * a ** -4 *
                            like_equal_d_cosi(a_hat, f, d0 / a, x),
                  max(a_peak - 5 * a_width, 0), a_peak + 5 * a_width,
                  epsabs=1e-2, epsrel=1e-4)[0]
    if make_plot:
        a = linspace(max(a_peak - 5 * a_width, 0), a_peak + 5 * a_width)
        lc = zeros_like(a)
        for (i, a) in enumerate(a):
            lc[i] = a ** -4 * like_equal_d_cosi(a_hat, f, d0 / a, x)
        pylab.figure()
        pylab.plot(a, lc)
    return l_circ


def like_equal(a_hat, f, d_max=1000):
    """calculate the likelihood for network equally sensitive to plus and cross.
    Marginalized over the two phases, inclination and distance.  We use uniform
    priors, specifically 1/2pi for psi,phi, 1/2 for cosi, 3 d^2dd/d_max^3 for
    distance d_max.  """
    lc = 0.5 * quad(lambda x: like_equal_cosi(a_hat, f, x, d_max),
                    -1., 1., epsabs=1e-2, epsrel=1e-4)[0]
    return lc




def like_parts_d_cosi_psi(a_hat, f_plus, f_cross, x, psi):
    """
    calculate the two dimensional likelihood, marginalized over phi
    log-likelihood can be written as:
    1/2(ahat^2 - 2*d0/d * f(x, psi) * cos(2phi - phi0) + (d0/d)^2 g(x,psi))
    return: ahat2, f, g
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param x: cos(inclination)
    :param psi: polarization
    """
    c2p = cos(2 * psi)
    s2p = sin(2 * psi)
    x2 = (1 + x ** 2) / 2
    c2phi_fac = a_hat[1] * f_plus ** 2 * x2 * c2p \
                + a_hat[2] * f_cross ** 2 * x2 * s2p \
                - a_hat[3] * f_plus ** 2 * x * s2p \
                + a_hat[4] * f_cross ** 2 * x * c2p
    s2phi_fac = - a_hat[1] * f_plus ** 2 * x * s2p \
                + a_hat[2] * f_cross ** 2 * x * c2p \
                - a_hat[3] * f_plus ** 2 * x2 * c2p \
                - a_hat[4] * f_cross ** 2 * x2 * s2p
    f = sqrt(c2phi_fac ** 2 + s2phi_fac ** 2)
    g = f_plus ** 2 * x2 ** 2 * c2p ** 2 + f_cross ** 2 * x2 ** 2 * s2p ** 2 \
        + f_plus ** 2 * x ** 2 * s2p ** 2 + f_cross ** 2 * x ** 2 * c2p ** 2
    f_resp = array([0, f_plus, f_cross, f_plus, f_cross])
    ahat2 = sum(f_resp ** 2 * a_hat ** 2)
    return ahat2, f, g


def like_d_cosi_psi(a_hat, f_plus, f_cross, d, x, psi, marg=True):
    """
    Return the likelihood marginalized over phi, using flat (1/2pi) prior
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param d: distance
    :param x: cos(inclination)
    :param psi: polarization
    :param marg: do or don't do the marginalization
    :returns: the marginalization factor.  I don't think it includes the exp(rho^2/2) term.
    """
    ahat2, f, g = like_parts_d_cosi_psi(a_hat, f_plus, f_cross, x, psi)
    d0 = a_hat[0]
    # Marginalizing over phi (from zero to 2 pi) gives:
    # 2 pi  i0e(a f) exp(-1/2(ahat^2 - 2 f a + g a^2))
    a = d0 / d
    like = exp(f * a - 0.5 * (ahat2 + g * a ** 2))
    if marg: like *= special.i0e(f * a)
    return like


def like_d_cosi(a_hat, f_plus, f_cross, d, x):
    """
    Return the likelihood marginalized over phi and psi, with a uniform
    (1/2pi) prior on both
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param d: distance
    :param x: cos(inclination)
    """
    l = 2 / pi * quad(lambda p: like_d_cosi_psi(a_hat, f_plus, f_cross, d, x, p),
                      0, pi / 2)[0]
    return l


def like_cosi_psi(a_hat, f_plus, f_cross, x, psi, d_max=1000.):
    """
    Return the likelihood marginalized over d and phi, using flat (1/2pi
    prior on phi; uniform volume on d up to d_max.
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param x: cos(inclination)
    :param psi: polarization
    :param d_max: maximum distance for marginalization
    """
    ahat2, f, g = like_parts_d_cosi_psi(a_hat, f_plus, f_cross, x, psi)
    d0 = a_hat[0]
    # Marginalizing over phi gives:
    # 2 pi  i0e(a f) exp(-1/2(ahat^2 - 2 f a + g a^2))
    like = lambda a: 3 * d0 ** 2 / d_max ** 3 * a ** (-4) * special.i0e(f * a) * \
                     exp(f * a - 0.5 * (ahat2 + g * a ** 2))
    l_psix = quad(like, max(0, f / g - 5 / sqrt(g)), f / g + 5 / sqrt(g))
    return l_psix[0]


def like_cosi(a_hat, f_plus, f_cross, x, d_max=1000.):
    """
    Return the likelihood marginalized over d, phi and psi.
    Use uniform (1/2pi) prior on phi, psi, uniform in volume prior over d out to
    d_max.
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param x: cos(inclination)
    :param d_max: maximum distance for marginalization
    """
    l_x = 2 / pi * quad(lambda p: like_cosi_psi(a_hat, f_plus, f_cross, x, p, d_max),
                        0, math.pi / 2, epsabs=1e-2, epsrel=1e-4)[0]
    return l_x


def like(a_hat, f_plus, f_cross, d_max=1000.):
    """
    Return the likelihood marginalized over all 4 f-stat parameteras.
    Use uniform (1/2pi) prior on phi, psi, uniform (1/2)on cosi,
    uniform in volume prior over d out to d_max.
    :param a_hat: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param d_max: maximum distance for marginalization
    """
    l = 1. / 2 * quad(lambda x: like_cosi(a_hat, f_plus, f_cross, x, d_max),
                      -1, 1, epsabs=1e-2, epsrel=1e-4)[0]
    return l


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
