# Copyright (C) 2020  Daniel Finstad
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
"""This module provides model classes and functions for implementing
a relative binning likelihood for parameter estimation.
"""


import logging
import numpy
import itertools
from scipy.interpolate import interp1d
from scipy.optimize import minimize

from pycbc.waveform import (get_fd_waveform_sequence,
                            get_fd_det_waveform_sequence, fd_det_sequence)
from pycbc.detector import Detector
from pycbc.types import Array, TimeSeries

from .base import ModelStats
from .gaussian_noise import (BaseGaussianNoise, catch_waveform_error)
from pycbc.waveform import FailedWaveformError
from .relbin_cpu import (likelihood_parts, likelihood_parts_v,
                         likelihood_parts_multi, likelihood_parts_multi_v,
                         likelihood_parts_det, likelihood_parts_det_multi,
                         likelihood_parts_vector,
                         likelihood_parts_v_pol,
                         likelihood_parts_v_time,
                         likelihood_parts_v_pol_time,
                         likelihood_parts_vectorp, snr_predictor,
                         likelihood_parts_vectort,
                         snr_predictor_dom)
from .tools import DistMarg


def setup_bins(f_full, f_lo, f_hi, chi=1.0,
               eps=0.1, gammas=None,
               ):
    """Construct frequency bins for use in a relative likelihood
    model. For details, see [Barak, Dai & Venumadhav 2018].

    Parameters
    ----------
    f_full : array
        The full resolution array of frequencies being used in the analysis.
    f_lo : float
        The starting frequency used in matched filtering. This will be the
        left edge of the first frequency bin.
    f_hi : float
        The ending frequency used in matched filtering. This will be the right
        edge of the last frequency bin.
    chi : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]
    eps : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]. Lower values
        result in larger number of bins.
    gammas : array, optional
        Frequency powerlaw indices to be used in computing bins.

    Returns
    -------
    nbin : int
        Number of bins.
    fbin : numpy.array of floats
        Bin edge frequencies.
    fbin_ind : numpy.array of ints
        Indices of bin edges in full frequency array.
    """
    f = numpy.linspace(f_lo, f_hi, 10000)
    # f^ga power law index
    ga = (
        gammas
        if gammas is not None
        else numpy.array([-5.0 / 3, -2.0 / 3, 1.0, 5.0 / 3, 7.0 / 3])
    )
    logging.info("Using powerlaw indices: %s", ga)
    dalp = chi * 2.0 * numpy.pi / numpy.absolute((f_lo ** ga) - (f_hi ** ga))
    dphi = numpy.sum(
        numpy.array([numpy.sign(g) * d * (f ** g) for g, d in zip(ga, dalp)]),
        axis=0,
    )
    dphi_diff = dphi - dphi[0]
    # now construct frequency bins
    nbin = int(dphi_diff[-1] / eps)
    dphi2f = interp1d(
        dphi_diff, f, kind="slinear", bounds_error=False, fill_value=0.0
    )
    dphi_grid = numpy.linspace(dphi_diff[0], dphi_diff[-1], nbin + 1)
    # frequency grid points
    fbin = dphi2f(dphi_grid)
    # indices of frequency grid points in the FFT array
    fbin_ind = numpy.searchsorted(f_full, fbin)
    for idx_fbin, idx_f_full in enumerate(fbin_ind):
        if idx_f_full == 0:
            curr_idx = 0
        elif idx_f_full == len(f_full):
            curr_idx = len(f_full) - 1
        else:
            abs1 = abs(f_full[idx_f_full] - fbin[idx_fbin])
            abs2 = abs(f_full[idx_f_full-1] - fbin[idx_fbin])
            if abs1 > abs2:
                curr_idx = idx_f_full - 1
            else:
                curr_idx = idx_f_full
        fbin_ind[idx_fbin] = curr_idx
    fbin_ind = numpy.unique(fbin_ind)
    return fbin_ind


def is_rescaled(base, other, tolerance=1e-6):
    """Whether two waveforms differ only by an overall constant.

    That is what tells a parameter which merely scales the waveform from
    one which changes how it varies with frequency.
    """
    keep = base != 0
    if not keep.any():
        return True
    ratio = other[keep] / base[keep]
    scale = numpy.median(abs(ratio))
    if scale == 0:
        return not numpy.any(other)
    return bool(abs(ratio - numpy.mean(ratio)).max() <= tolerance * scale)


class Relative(DistMarg, BaseGaussianNoise):
    r"""Model that assumes the likelihood in a region around the peak
    is slowly varying such that a linear approximation can be made, and
    likelihoods can be calculated at a coarser frequency resolution. For
    more details on the implementation, see https://arxiv.org/abs/1806.08792.

    This model requires the use of a fiducial waveform whose parameters are
    near the peak of the likelihood. The fiducial waveform and all template
    waveforms used in likelihood calculation are currently generated using
    the SPAtmplt approximant.

    For more details on initialization parameters and definition of terms, see
    :py:class:`BaseGaussianNoise`.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). All data must have the
        same frequency resolution.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    figucial_params : dict
        A dictionary of waveform parameters to be used for generating the
        fiducial waveform. Keys must be parameter names in the form
        'PARAM_ref' where PARAM is a recognized extrinsic parameter or
        an intrinsic parameter compatible with the chosen approximant.
    gammas : array of floats, optional
        Frequency powerlaw indices to be used in computing frequency bins.
    epsilon : float or 'auto', optional
        Tuning parameter used in calculating the frequency bins. Lower values
        will result in higher resolution and more bins. If 'auto', the bins
        are laid out again while running whenever they turn out not to be
        good enough for where the sampler has got to.
    accuracy : float, optional
        How much error in the log likelihood ratio to accept when 'auto' is
        deciding whether to lay the bins out again.
    earth_rotation: boolean, optional
        Default is False. If True, then vary the fp/fc polarization values
        as a function of frequency bin, using a predetermined PN approximation
        for the time offsets.
    optimize_fiducial: boolean, optional
        Default is False. If True, treat the fiducial parameters given as a
        rough starting point, maximize from there, and build the fiducial
        waveform again at what that finds. Saves having to know where the
        peak is beforehand.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise`.
    """
    name = "relative"

    def __init__(
        self,
        variable_params,
        data,
        low_frequency_cutoff,
        fiducial_params=None,
        gammas=None,
        epsilon=0.5,
        earth_rotation=False,
        earth_rotation_mode=2,
        marginalize_phase=True,
        optimize_fiducial=False,
        accuracy=0.05,
        **kwargs
    ):

        variable_params, kwargs = self.setup_marginalization(
                               variable_params,
                               marginalize_phase=marginalize_phase,
                               **kwargs)

        super(Relative, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs
        )

        # If the waveform needs us to apply the detector response,
        # set flag to true (most cases for ground-based observatories).
        self.still_needs_det_response = False
        if self.static_params['approximant'] in fd_det_sequence:
            self.still_needs_det_response = True

        # reference waveform and bin edges
        self.f, self.df, self.end_time, self.det = {}, {}, {}, {}
        self.h00, self.h00_sparse = {}, {}
        self.fedges, self.edges = {}, {}
        self.ta, self.antenna_time = {}, {}

        # filtered summary data for linear approximation
        self.sdat = {}

        # store fiducial waveform params
        self.fid_params = self.static_params.copy()
        self.fid_params.update(fiducial_params)

        # the flag used in `_loglr`
        self.return_sh_hh = False

        for k in self.static_params:
            if self.fid_params[k] == 'REPLACE':
               self.fid_params.pop(k)

        # the data with the fiducial waveform's timing taken out, and the
        # frequency range it covers. Both are what the bins are laid out
        # from, and are kept so that can be done again.
        self.shifted, self.flimits = {}, {}

        self.adapt = str(epsilon).lower() == 'auto'
        self.epsilon = 0.5 if self.adapt else float(epsilon)
        self.accuracy = float(accuracy)
        self.layout = (gammas, earth_rotation, int(earth_rotation_mode))
        self.best_loglr = -numpy.inf
        self.since_check = 0
        self.interval = 500
        self.rebins = 0

        self.setup_fiducial()
        self.setup_bin_layout(self.epsilon, *self.layout)

        if optimize_fiducial:
            for _ in range(5):
                found = self.optimize_fiducial_params()
                if not found:
                    break
                moved = max((abs(v - self.fid_params[p]) / max(abs(v), 1e-8)
                             for p, v in found.items()
                             if p in self.fid_params), default=numpy.inf)
                self.fid_params.update(found)
                self.setup_fiducial()
                self.setup_bin_layout(self.epsilon, *self.layout)
                if moved < 1e-3:
                    break
            else:
                logging.warning("The fiducial waveform was still moving "
                                "after 5 rounds of maximizing; it may not "
                                "be near the peak")

    def setup_fiducial(self):
        """Generate the fiducial waveform and the data to compare against it.

        Redone when the fiducial parameters change. The bins are laid out
        from what this leaves behind, so they have to be laid out again
        after it.
        """
        for ifo in self.data:
            # store data and frequencies
            d0 = self.data[ifo]
            self.f[ifo] = d0.sample_frequencies.numpy()
            self.df[ifo] = d0.delta_f
            self.end_time[ifo] = float(d0.end_time)

            # generate fiducial waveform
            f_lo = self.kmin[ifo] * self.df[ifo]
            f_hi = self.kmax[ifo] * self.df[ifo]
            logging.info(
                "%s: Generating fiducial waveform from %s to %s Hz",
                ifo, f_lo, f_hi,
            )

            # prune low frequency samples to avoid waveform errors
            fpoints = Array(self.f[ifo].astype(numpy.float64))
            fpoints = fpoints[self.kmin[ifo]:self.kmax[ifo]+1]

            if self.still_needs_det_response:
                wave = get_fd_det_waveform_sequence(ifos=ifo,
                                                    sample_points=fpoints,
                                                    **self.fid_params)
                curr_wav = wave[ifo]
                self.ta[ifo] = 0.
            else:
                fid_hp, fid_hc = get_fd_waveform_sequence(sample_points=fpoints,
                                                          **self.fid_params)
                # Apply detector response if not handled by
                # the waveform generator. Reference the detector at the
                # time being analyzed: Detector estimates sidereal time by
                # advancing it at a constant rate from its reference, which
                # defaults to the time of GW150914, so a detector left at
                # that default drifts further from the truth the further the
                # data sits from 2015. The sky draws already reference the
                # event time, so leaving this at the default also made the
                # two disagree.
                self.det[ifo] = Detector(
                    ifo, reference_time=self.fid_params["tc"])
                # the combined call, from detector-response-split, gets the
                # pattern and the delay from one sidereal time calculation
                fp, fc, dt = self.det[ifo].antenna_pattern_and_delay(
                    self.fid_params["ra"], self.fid_params["dec"],
                    self.fid_params["polarization"], self.fid_params["tc"])
                self.ta[ifo] = self.fid_params["tc"] + dt
                curr_wav = (fid_hp * fp + fid_hc * fc)

            # check for zeros at low and high frequencies
            # make sure only nonzero samples are included in bins
            numzeros_lo = list(curr_wav != 0j).index(True)
            if numzeros_lo > 0:
                new_kmin = self.kmin[ifo] + numzeros_lo
                f_lo = new_kmin * self.df[ifo]
                logging.info(
                    "WARNING! Fiducial waveform starts above "
                    "low-frequency-cutoff, initial bin frequency "
                    "will be %s Hz", f_lo)
            numzeros_hi = list(curr_wav[::-1] != 0j).index(True)
            if numzeros_hi > 0:
                new_kmax = self.kmax[ifo] - numzeros_hi
                f_hi = new_kmax * self.df[ifo]
                logging.info(
                    "WARNING! Fiducial waveform terminates below "
                    "high-frequency-cutoff, final bin frequency "
                    "will be %s Hz", f_hi)

            self.ta[ifo] -= self.end_time[ifo]
            curr_wav.resize(len(self.f[ifo]))
            curr_wav = numpy.roll(curr_wav, self.kmin[ifo])

            # We'll apply this to the data, in lieu of the ref waveform
            # This makes it easier to compare target signal to reference later
            tshift = numpy.exp(-2.0j * numpy.pi * self.f[ifo] * self.ta[ifo])
            self.h00[ifo] = numpy.array(curr_wav) # * tshift
            self.shifted[ifo] = self.data[ifo] * numpy.conjugate(tshift)
            self.flimits[ifo] = (f_lo, f_hi)

    def setup_bin_layout(self, epsilon, gammas, earth_rotation,
                         earth_rotation_mode):
        """Place the frequency bins and compute the summary data for them.

        Everything this needs is left behind by
        :py:meth:`setup_fiducial`, so it can be redone on its own to
        change the resolution without generating the fiducial waveform
        again.
        """
        for ifo in self.data:
            f_lo, f_hi = self.flimits[ifo]

            logging.info("Computing frequency bins")
            fbin_ind = setup_bins(
                f_full=self.f[ifo], f_lo=f_lo, f_hi=f_hi,
                gammas=gammas, eps=float(epsilon),
            )
            logging.info("Using %s bins for this model", len(fbin_ind))

            self.fedges[ifo] = self.f[ifo][fbin_ind]
            self.edges[ifo] = fbin_ind
            self.init_from_frequencies(self.shifted[ifo], self.h00,
                                       fbin_ind, ifo)
            self.antenna_time[ifo] = self.setup_antenna(
                                        earth_rotation,
                                        earth_rotation_mode,
                                        self.fedges[ifo])
        self.combine_layout()
        self.check_bin_resolution()

    def optimize_fiducial_params(self, ndraw=200, seed=0):
        """Look for fiducial parameters near the peak of the likelihood.

        Relative binning is accurate where the waveform stays close in
        shape to the fiducial one, so the fiducial wants to sit near the
        posterior. This starts from a rough set of parameters drawn from
        the prior and maximizes from there, which reuses the waveform
        machinery already set up rather than adding any.

        What is maximized is the signal to noise ratio rather than the
        likelihood. It is normalized by the amplitude of the waveform, so
        parameters that only scale the waveform, such as the distance, do
        not enter it. Maximizing the likelihood instead walks along the
        degeneracy between those parameters and stops at whatever bound it
        reaches, which says nothing about the shape the bins are placed
        for. Parameters that leave the normalized waveform alone are then
        left where they started rather than searched over.

        The model used is the one built from the current fiducial
        waveform, so it is only a guide; the caller rebuilds with the
        result, and the answer improves rather than degrades because the
        new fiducial is closer to the signal than the one it came from.

        Parameters
        ----------
        ndraw : int, optional
            How many prior draws to start the maximization from.
        seed : int, optional
            Seed for the draws, so the choice is reproducible.

        Returns
        -------
        dict
            The parameters to use for the fiducial waveform.
        """
        if self.prior_distribution is None:
            return {}
        params = [p for p in self.variable_params
                  if p in self.prior_distribution.variable_args]
        if not params:
            return {}

        state = numpy.random.get_state()
        numpy.random.seed(seed)
        try:
            draws = self.prior_distribution.rvs(size=int(ndraw))
        finally:
            numpy.random.set_state(state)

        saved = (self._current_params, self._current_stats,
                 self.return_sh_hh)
        self.return_sh_hh = True
        try:
            def snr(point):
                """Minus the network signal to noise ratio squared."""
                if self.prior_distribution(**point) == -numpy.inf:
                    return numpy.inf
                try:
                    self.update(**point)
                    sh, hh = self.loglr
                    return -abs(sh) ** 2 / hh
                except Exception:  # a point we cannot evaluate is no peak
                    return numpy.inf

            start = min(({p: draw[p] for p in params} for draw in draws),
                        key=snr)
            shape = self.shape_params(params, draws)
            if not shape:
                logging.warning("No parameter changes the shape of the "
                                "waveform, keeping the fiducial given")
                return {}

            def objective(values):
                return snr(dict(start, **dict(zip(shape, values))))

            result = minimize(objective, [start[p] for p in shape],
                              method='Nelder-Mead')
        finally:
            (self._current_params, self._current_stats,
             self.return_sh_hh) = saved

        if not numpy.isfinite(result.fun):
            logging.warning("Could not find a fiducial waveform by "
                            "maximizing the signal to noise ratio, keeping "
                            "the one given")
            return {}

        found = dict(start, **dict(zip(shape, result.x)))
        logging.info("Chose a fiducial waveform with a signal to noise "
                     "ratio of %.4g at %s", (-result.fun) ** 0.5,
                     ", ".join("%s=%.4g" % kv for kv in found.items()))
        return found

    def shape_params(self, params, draws):
        """Which of these parameters change how the waveform varies with
        frequency.

        Bins are placed to follow that variation, so it is all the fiducial
        waveform has to get right. A parameter that only multiplies the
        waveform by a constant, such as the distance, leaves it alone
        whatever value it takes, and searching over one only walks along a
        degeneracy until it reaches a bound.

        Found by moving one parameter at a time between prior draws and
        comparing the waveforms, so nothing here needs to know which
        parameters a given approximant treats as an overall scale.
        """
        shape = []
        # a handful of pairs is plenty; two draws agreeing by chance is
        # what the rest are there to catch
        for a, b in zip(draws[:5], draws[1:6], strict=False):
            point = {p: a[p] for p in params}
            base = self.waveform_profile(point)
            if base is None:
                continue
            for p in params:
                if p in shape:
                    continue
                other = self.waveform_profile(dict(point, **{p: b[p]}))
                if other is None or any(
                        not is_rescaled(x, y)
                        for x, y in zip(base, other, strict=False)):
                    shape.append(p)
            if len(shape) == len(params):
                break
        return shape

    def waveform_profile(self, point):
        """The waveform at the bin edges, or None if it cannot be made."""
        try:
            self.update(**point)
            wfs = self.get_waveforms(self.current_params)
        except Exception:
            return None
        profile = []
        for ifo in sorted(wfs):
            value = wfs[ifo]
            # with the detector response applied there is one array per
            # detector, without it a polarization pair
            for part in (value if isinstance(value, tuple) else (value,)):
                profile.append(numpy.asarray(part))
        return profile

    def update(self, **params):
        if self.adapt:
            self.keep_bins_good()
        super().update(**params)

    def loglr_at(self, point):
        """The log likelihood ratio at a point, leaving the model be."""
        saved = (self._current_params, self._current_stats)
        try:
            self._current_params = point
            self._current_stats = ModelStats()
            return float(self.loglr)
        except Exception:
            return numpy.nan
        finally:
            self._current_params, self._current_stats = saved

    def keep_bins_good(self, drop=20.0, check_every=500, smallest=0.005):
        """Lay the bins out again if they are not good enough for where
        the sampler has got to.

        The resolution the bins need depends on where the likelihood is
        being asked about, which is not known until it is being asked.
        This looks at the point just evaluated: if it is within ``drop`` of
        the best seen, so somewhere the answer matters, the resolution is
        checked and the bins are refined if that costs more than the
        accuracy wanted.

        The check costs about thirty evaluations, so it is not done on
        every one; laying the bins out again costs a few hundred, against
        the millions in a run. Refining only ever adds bins, so this
        settles rather than going back and forth.
        """
        previous = getattr(self._current_stats, 'loglr', None)
        if previous is None or not numpy.isfinite(previous):
            return
        previous = float(previous)
        self.best_loglr = max(self.best_loglr, previous)

        self.since_check += 1
        if (previous < self.best_loglr - drop
                or self.since_check < self.interval
                or self.epsilon <= smallest):
            return
        self.since_check = 0

        # the cheap screen: predict what the interpolation error is worth
        # in the likelihood, which costs no new bin layout. The error in
        # the log likelihood ratio is about the error in the waveform
        # ratio times the signal to noise ratio, and that ratio is about
        # the square root of twice the log likelihood ratio, so a small
        # departure at a loud signal is not small in the likelihood. A
        # fixed threshold on the ratio alone would pass a loud signal that
        # is badly under resolved, and waste checks on a quiet one where
        # the cost is genuinely small. Skipping only when the predicted
        # cost is below the accuracy wanted handles both.
        snr = (2.0 * max(previous, 0.0)) ** 0.5
        if snr * self.interpolation_error_from_reference() < self.accuracy:
            self.interval = min(self.interval * 2, 32000)
            return

        # it looks bad, so find out what it is worth in the likelihood
        point = dict(self._current_params)
        self.setup_bin_layout(self.epsilon / 2., *self.layout)
        finer = self.loglr_at(point)
        if numpy.isfinite(finer) and abs(finer - previous) < self.accuracy:
            # the coarser bins were good enough after all
            self.setup_bin_layout(self.epsilon, *self.layout)
            self.interval = min(self.interval * 2, 32000)
            return

        self.epsilon /= 2.
        self.interval = check_every
        self.rebins += 1
        logging.info("Refined epsilon to %.4g, %s bins: the likelihood "
                     "moved by %.3g against an accuracy of %.3g",
                     self.epsilon, len(self.fedges[list(self.data)[0]]),
                     abs(finer - previous), self.accuracy)

    def init_from_frequencies(self, data, h00, fbin_ind, ifo):
        bins = numpy.array(
            [
                (fbin_ind[i], fbin_ind[i + 1])
                for i in range(len(fbin_ind) - 1)
            ]
        )

        # store low res copy of fiducial waveform
        self.h00_sparse[ifo] = h00[ifo].copy().take(fbin_ind)

        # compute summary data
        logging.info(
            "Calculating summary data at frequency resolution %s Hz",
            self.df[ifo],
        )

        a0, a1 = self.summary_product(data, h00[ifo], bins, ifo)
        b0, b1 = self.summary_product(h00[ifo], h00[ifo], bins, ifo)
        self.sdat[ifo] = {"a0": a0, "a1": a1, "b0": abs(b0), "b1": abs(b1)}

    def combine_layout(self):
        # determine the unique ifo layouts
        self.edge_unique = []
        self.ifo_map = {}
        for ifo in self.fedges:
            if len(self.edge_unique) == 0:
                self.ifo_map[ifo] = 0
                self.edge_unique.append(Array(self.fedges[ifo]))
            else:
                for i, edge in enumerate(self.edge_unique):
                    if numpy.array_equal(edge, self.fedges[ifo]):
                        self.ifo_map[ifo] = i
                        break
                else:
                    self.ifo_map[ifo] = len(self.edge_unique)
                    self.edge_unique.append(Array(self.fedges[ifo]))
        logging.info("%s unique ifo layouts", len(self.edge_unique))

    def setup_antenna(self, earth_rotation, mode, fedges):
        # Calculate the times to evaluate fp/fc
        self.earth_rotation = earth_rotation
        if earth_rotation is not False:
            logging.info("Enabling frequency-dependent earth rotation")
            from pycbc.waveform.spa_tmplt import spa_length_in_time

            times = spa_length_in_time(
                phase_order=-1,
                mass1=self.fid_params["mass1"],
                mass2=self.fid_params["mass2"],
                f_lower=numpy.array(fedges) / mode * 2.0,
            )
            atimes = self.fid_params["tc"] - times
            self.lik = likelihood_parts_v
            self.mlik = likelihood_parts_multi_v
        else:
            atimes = self.fid_params["tc"]
            if self.still_needs_det_response:
                self.lik = likelihood_parts_det
                self.mlik = likelihood_parts_det_multi
            else:
                self.lik = likelihood_parts
                self.mlik = likelihood_parts_multi
        return atimes

    @property
    def likelihood_function(self):
        self.lformat = None
        if self.marginalize_vector_params:
            p = self.current_params

            vmarg = set(k for k in self.marginalize_vector_params
                        if not numpy.isscalar(p[k]))

            if self.earth_rotation:
                if set(['tc', 'polarization']).issubset(vmarg):
                    self.lformat = 'earth_time_pol'
                    return likelihood_parts_v_pol_time
                elif set(['polarization']).issubset(vmarg):
                    self.lformat = 'earth_pol'
                    return likelihood_parts_v_pol
                elif set(['tc']).issubset(vmarg):
                    self.lformat = 'earth_time'
                    return likelihood_parts_v_time
            else:
                if set(['ra', 'dec', 'tc']).issubset(vmarg):
                    return likelihood_parts_vector
                elif set(['tc', 'polarization']).issubset(vmarg):
                    return likelihood_parts_vector
                elif set(['tc']).issubset(vmarg):
                    return likelihood_parts_vectort
                elif set(['polarization']).issubset(vmarg):
                    return likelihood_parts_vectorp

        return self.lik

    def summary_product(self, h1, h2, bins, ifo):
        """ Calculate the summary values for the inner product <h1|h2>
        """
        # calculate coefficients
        h12 = numpy.conjugate(h1) * h2 / self.psds[ifo]

        # constant terms
        a0 = numpy.array([
                4.0 * self.df[ifo] * h12[l:h].sum()
                for l, h in bins
            ])

        # linear terms
        a1 = numpy.array([
                4.0 / (h - l) *
                (h12[l:h] * (self.f[ifo][l:h] - self.f[ifo][l])).sum()
                for l, h in bins])

        return a0, a1

    def get_waveforms(self, params):
        """ Get the waveform polarizations for each ifo
        """
        if self.still_needs_det_response:
            wfs = {}
            for ifo in self.data:
                wfs.update(get_fd_det_waveform_sequence(
                        ifos=ifo, sample_points=self.fedges[ifo], **params))
                if self.recalibration:
                    wfs[ifo] = wfs[ifo] * self.calibration_factor(ifo, params)
            return wfs

        wfs = []
        for edge in self.edge_unique:
            hp, hc = get_fd_waveform_sequence(sample_points=edge, **params)
            hp = hp.numpy()
            hc = hc.numpy()
            wfs.append((hp, hc))
        wf_ret = {ifo: wfs[self.ifo_map[ifo]] for ifo in self.data}

        if self.recalibration:
            # the correction differs by detector, so it cannot be shared
            # between detectors the way the waveform itself is
            wf_ret = {ifo: tuple(h * self.calibration_factor(ifo, params)
                                 for h in pols)
                      for ifo, pols in wf_ret.items()}

        self.wf_ret = wf_ret
        return wf_ret

    def calibration_factor(self, ifo, params):
        """The calibration correction, at this detector's bin edges only.

        The correction is a smooth function of frequency, so it is as well
        described across a bin as the waveform ratio is, and evaluating it
        at the edges costs one spline per bin rather than one per frequency
        sample of the data.
        """
        recalib = self.recalibration[ifo]
        recalib.set_params(**params)
        return recalib.calibration_factor(self.fedges[ifo])

    @property
    def multi_signal_support(self):
        """ The list of classes that this model supports in a multi-signal
        likelihood
        """
        # Check if this model *can* be included in a multi-signal model.
        # All marginalizations must currently be disabled to work: the
        # marginalized likelihood is not the sum of the components', so
        # combining the two gives a wrong answer rather than a slow one.
        # This used to only log the fact and report support anyway, which
        # left the wrong answer reachable without any sign of it.
        if (self.marginalize_vector_params or
                self.marginalize_distance or
                self.marginalize_phase):
            raise ValueError(
                "Cannot use %s inside of multi_signal while any "
                "marginalization is enabled; the combination silently gives "
                "a wrong likelihood. Disable marginalize_phase, "
                "marginalize_distance and marginalize_vector_params, or use "
                "the joint_primary_marginalized model instead."
                % type(self).__name__)
        return [type(self)]

    def calculate_hihjs(self, models):
        """ Pre-calculate the hihj inner products on a grid
        """
        self.hihj = {}
        for m1, m2 in itertools.combinations(models, 2):
            self.hihj[(m1, m2)] = {}
            for ifo in self.data:
                h1 = m1.h00[ifo]
                h2 = m2.h00[ifo]

                # Combine the grids
                edge = numpy.unique([m1.edges[ifo], m2.edges[ifo]])

                # Remove any points where either reference is zero
                keep = numpy.where((h1[edge] != 0) | (h2[edge] != 0))[0]
                edge = edge[keep]
                fedge = m1.f[ifo][edge]

                bins = numpy.array([
                        (edge[i], edge[i + 1])
                        for i in range(len(edge) - 1)
                    ])
                a0, a1 = self.summary_product(h1, h2, bins, ifo)
                self.hihj[(m1, m2)][ifo] = a0, a1, fedge

    def multi_loglikelihood(self, models):
        """ Calculate a multi-model (signal) likelihood
        """
        models = [self] + models
        loglr = 0
        # handle sum[<d|h_i> - 0.5 <h_i|h_i>]
        for m in models:
            loglr += m.loglr

        if not hasattr(self, 'hihj'):
            self.calculate_hihjs(models)

        if self.still_needs_det_response:
            for m1, m2 in itertools.combinations(models, 2):
                for det in self.data:
                    a0, a1, fedge = self.hihj[(m1, m2)][det]

                    dtc, channel, h00 = m1._current_wf_parts[det]
                    dtc2, channel2, h002 = m2._current_wf_parts[det]

                    c1c2 = self.mlik(fedge,
                                     dtc, channel, h00,
                                     dtc2, channel2, h002,
                                     a0, a1)
                    loglr += - c1c2.real # This is -0.5 * re(<h1|h2> + <h2|h1>)
        else:
            # finally add in the lognl term from this model
            for m1, m2 in itertools.combinations(models, 2):
                for det in self.data:
                    a0, a1, fedge = self.hihj[(m1, m2)][det]

                    fp, fc, dtc, hp, hc, h00 = m1._current_wf_parts[det]
                    fp2, fc2, dtc2, hp2, hc2, h002 = m2._current_wf_parts[det]

                    h1h2 = self.mlik(fedge,
                                     fp, fc, dtc, hp, hc, h00,
                                     fp2, fc2, dtc2, hp2, hc2, h002,
                                     a0, a1)
                    loglr += - h1h2.real # This is -0.5 * re(<h1|h2> + <h2|h1>)
        return loglr + self.lognl

    @catch_waveform_error
    def _loglr(self):
        r"""Computes the log likelihood ratio,
        or inner product <s|h> and <h|h> if `self.return_sh_hh` is True.

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        or
        tuple
            The inner product (<s|h>, <h|h>).
        """
        # get model params
        p = self.current_params
        wfs = self.get_waveforms(p)
        lik = self.likelihood_function
        norm = 0.0
        filt = 0j
        self._current_wf_parts = {}
        pol_phase = numpy.exp(-2.0j * p['polarization'])

        for ifo in self.data:
            freqs = self.fedges[ifo]
            sdat = self.sdat[ifo]
            h00 = self.h00_sparse[ifo]
            end_time = self.end_time[ifo]
            times = self.antenna_time[ifo]

            # project waveform to detector frame if waveform does not deal
            # with detector response. Otherwise, skip detector response.

            if self.still_needs_det_response:
                dtc = 0.

                channel = wfs[ifo].numpy()
                filter_i, norm_i = lik(freqs, dtc, channel, h00,
                                       sdat['a0'], sdat['a1'],
                                       sdat['b0'], sdat['b1'])
                self._current_wf_parts[ifo] = (dtc, channel, h00)
            else:
                hp, hc = wfs[ifo]
                det = self.det[ifo]
                fp, fc, dt = det.antenna_pattern_and_delay(
                    p["ra"], p["dec"], 0.0, times)
                dtc = p["tc"] + dt - end_time - self.ta[ifo]

                if self.lformat == 'earth_pol':
                    filter_i, norm_i = lik(freqs, fp, fc, dtc, pol_phase,
                                           hp, hc, h00,
                                           sdat['a0'], sdat['a1'],
                                           sdat['b0'], sdat['b1'])
                else:
                    f = (fp + 1.0j * fc) * pol_phase
                    fp = f.real.copy()
                    fc = f.imag.copy()
                    filter_i, norm_i = lik(freqs, fp, fc, dtc,
                                           hp, hc, h00,
                                           sdat['a0'], sdat['a1'],
                                           sdat['b0'], sdat['b1'])
                    self._current_wf_parts[ifo] = (fp, fc, dtc, hp, hc, h00)

            filt += filter_i
            norm += norm_i

        loglr = self.marginalize_loglr(filt, norm)
        if self.return_sh_hh:
            results = (filt, norm)
        else:
            results = loglr
        return results

    def _nowaveform_handler(self):
        """Returns -inf for loglr if no waveform generated.

        If `return_sh_hh` is set to True, a FailedWaveformError will be raised.
        """
        if self.return_sh_hh:
            raise FailedWaveformError("Waveform failed to generate and "
                                      "return_sh_hh set to True! I don't know "
                                      "what to return in this case.")
        return -numpy.inf

    def write_metadata(self, fp, group=None):
        """Adds writing the fiducial parameters and epsilon to file's attrs.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        super().write_metadata(fp, group=group)
        if group is None:
            attrs = fp.attrs
        else:
            attrs = fp[group].attrs
        for p, v in self.fid_params.items():
            attrs["{}_ref".format(p)] = v

    def interpolation_error_from_reference(self):
        """ Return the largest error made by interpolating the waveform
        ratio across a bin, relative to the size of the ratio.

        Relative binning assumes the ratio to the fiducial waveform is
        linear within a bin. This evaluates the ratio at the bin midpoints
        and compares it to that linear model, so it measures the
        approximation directly rather than estimating it from the spacing
        of the edges.
        """
        worst = 0.
        for ifo in self.data:
            idx = numpy.asarray(self.edges[ifo])
            mid = (idx[:-1] + idx[1:]) // 2
            # bins only one sample wide are exact, and have no midpoint
            keep = (mid > idx[:-1]) & (self.h00[ifo][mid] != 0)
            if not keep.any():
                continue

            fmid = self.f[ifo][mid]
            hp, _ = get_fd_waveform_sequence(sample_points=Array(fmid),
                                             **self.current_params)
            ratio = hp.numpy() / self.h00[ifo][mid]

            edge = self.wf_ret[ifo][0] / self.h00_sparse[ifo]
            fedges = self.fedges[ifo]
            w = (fmid - fedges[:-1]) / numpy.diff(fedges)
            linear = edge[:-1] * (1 - w) + edge[1:] * w

            err = abs(ratio - linear)[keep] / abs(ratio)[keep]
            worst = max(worst, err.max())
        return worst

    def check_bin_resolution(self, ndraw=10, threshold=1e-3, seed=0):
        """Check that the frequency bins resolve the waveform ratio.

        Relative binning assumes the ratio of the waveform to the fiducial
        one is linear across a bin. This draws a few points from the prior
        and asks how far the ratio at the bin midpoints departs from that
        assumption; where it departs a lot, a bin was needed in between.

        It costs two waveforms per draw at the resolution of the bins
        rather than of the data, so it is cheap compared with the analysis
        it precedes.

        Parameters
        ----------
        ndraw : int, optional
            How many points to draw from the prior.
        threshold : float, optional
            Warn if the interpolation error exceeds this. On GW170817,
            errors of 9.6e-3, 2.3e-3, 5.9e-4 and 2.3e-5 came with errors in
            the log likelihood ratio of 0.31, 0.056, 0.015 and 0.0013, so
            the cost is roughly 25 times the error reported here. That
            factor grows with the signal to noise ratio, so this indicates
            trouble rather than bounding the error.
        seed : int, optional
            Seed for the draws, so the check is reproducible.

        Returns
        -------
        float
            The largest error found, or 0 if it could not be checked.
        """
        if self.prior_distribution is None or self.still_needs_det_response:
            # nothing to draw from, or the ratio is not formed here
            return 0.
        state = numpy.random.get_state()
        numpy.random.seed(seed)
        try:
            draws = self.prior_distribution.rvs(size=int(ndraw))
        finally:
            numpy.random.set_state(state)

        # a diagnostic should not move the model out from under the caller
        saved = (self._current_params, self._current_stats)
        worst = 0.
        try:
            for draw in draws:
                params = {p: draw[p] for p in self.variable_params}
                try:
                    self.update(**params)
                    self.get_waveforms(self.current_params)
                    worst = max(worst,
                                self.interpolation_error_from_reference())
                except Exception:  # a bad draw should not stop the analysis
                    continue
        finally:
            self._current_params, self._current_stats = saved

        if worst > threshold:
            logging.warning(
                "The frequency bins may be too coarse: interpolating the "
                "waveform ratio across a bin is off by %.3g, against a "
                "threshold of %.3g. Consider a smaller epsilon, or a "
                "fiducial waveform closer to where the posterior is.",
                worst, threshold)
        else:
            logging.info("Bin resolution check: interpolation error %.3g "
                         "(threshold %.3g)", worst, threshold)
        return worst

    @staticmethod
    def extra_args_from_config(cp, section, skip_args=None, dtypes=None):
        """Adds reading fiducial waveform parameters from config file."""
        # add fiducial params to skip list
        skip_args += [
            option for option in cp.options(section) if option.endswith("_ref")
        ]

        # get frequency power-law indices if specified
        # NOTE these should be supplied in units of 1/3
        gammas = None
        if cp.has_option(section, "gammas"):
            skip_args.append("gammas")
            gammas = numpy.array(
                [float(g) / 3.0 for g in cp.get(section, "gammas").split()]
            )
        args = super(Relative, Relative).extra_args_from_config(
            cp, section, skip_args=skip_args, dtypes=dtypes
        )

        # get fiducial params from config
        fid_params = {
            p.replace("_ref", ""): float(cp.get("model", p))
            for p in cp.options("model")
            if p.endswith("_ref")
        }

        # add optional params with default values if not specified
        opt_params = {
            "ra": numpy.pi,
            "dec": 0.0,
            "inclination": 0.0,
            "polarization": numpy.pi,
        }
        fid_params.update(
            {p: opt_params[p] for p in opt_params if p not in fid_params}
        )
        args.update({"fiducial_params": fid_params, "gammas": gammas})
        return args


class RelativeTime(Relative):
    """ Heterodyne likelihood optimized for time marginalization. In addition
    it supports phase (dominant-mode), sky location, and polarization
    marginalization.
    """
    name = "relative_time"

    def __init__(self, *args,
                 sample_rate=4096,
                 **kwargs):
        super(RelativeTime, self).__init__(*args, **kwargs)
        auto = str(sample_rate).lower() == 'auto'
        if auto:
            # start from the data, which is as coarse as it can sensibly be
            sample_rate = self.data[tuple(self.data)[0]].sample_rate
        self.set_sample_rate(sample_rate, **kwargs)
        if auto:
            self.choose_sample_rate(**kwargs)
        self.draw_ifos(self.ref_snr, **kwargs)

    def set_sample_rate(self, rate, **kwargs):
        """Set the time resolution and redo what was built for the old one.
        """
        self.sample_rate = float(rate)
        if hasattr(self, '_ref_snr'):
            del self._ref_snr
        self.setup_peak_lock(sample_rate=self.sample_rate, **kwargs)

    def resolved_samples(self, snrs):
        """How many time samples share the weight of the SNR peak.

        The times are drawn in proportion to the likelihood along these
        series, so this is the number of samples that drawing actually has
        to choose between. One sample means the peak falls between grid
        points as often as not, and the answer depends on which.
        """
        worst = numpy.inf
        for ifo in snrs:
            logw = snrs[ifo].squared_norm().numpy() / 2.0
            weight = numpy.exp(logw - logw.max())
            worst = min(worst, weight.sum() ** 2 / (weight ** 2).sum())
        return worst

    def choose_sample_rate(self, target=2.5, most=65536.0, **kwargs):
        """Pick a time resolution that resolves the peak of the likelihood.

        The peak is narrow, and how narrow depends on the signal, so a
        fixed rate is either wasteful or wrong. Where the series already
        resolves the peak its width can be measured, and the rate that
        follows from it reached in one step. Where it does not there is
        nothing to measure, and the rate is doubled until there is.

        The rate is worth no more than it buys: the signal to noise series
        is rebuilt at every likelihood evaluation, so its length is a
        running cost rather than a setup one.

        Parameters
        ----------
        target : float, optional
            How many samples the peak should be spread over. Measured on a
            binary neutron star injection: at one sample the answer is as
            likely to be off by 0.17 as to be right, while from two
            samples upwards it was within 0.007 of direct integration at
            every rate tried.
        most : float, optional
            Never go above this rate.
        """
        resolved = self.resolved_samples(self.ref_snr)
        while resolved < target and self.sample_rate < most:
            if resolved >= 2.0:
                # resolved well enough to say how much more is needed
                steps = max(1, int(numpy.ceil(numpy.log2(target / resolved))))
            else:
                steps = 1
            self.set_sample_rate(min(self.sample_rate * 2 ** steps, most),
                                 **kwargs)
            resolved = self.resolved_samples(self.ref_snr)

        rate = self.sample_rate
        if resolved < target:
            logging.warning("Could not resolve the peak of the likelihood in "
                            "time: %.1f samples across it at %s Hz, against "
                            "the %s wanted", resolved, rate, target)
        else:
            logging.info("Chose a sample rate of %s Hz, spreading the peak "
                         "over %.1f samples", rate, resolved)
        return rate

    @property
    def ref_snr(self):
        if not hasattr(self, '_ref_snr'):
            wfs = {ifo: (self.h00_sparse[ifo],
                         self.h00_sparse[ifo]) for ifo in self.h00_sparse}
            self._ref_snr = self.get_snr(wfs)
        return self._ref_snr

    def get_snr(self, wfs):
        """ Return hp/hc maximized SNR time series
        """
        delta_t = 1.0 / self.sample_rate
        snrs = {}
        for ifo in wfs:
            sdat = self.sdat[ifo]
            dtc = self.tstart[ifo] - self.end_time[ifo] - self.ta[ifo]

            snr = snr_predictor(self.fedges[ifo],
                                dtc - delta_t * 2.0, delta_t,
                                self.num_samples[ifo] + 4,
                                wfs[ifo][0], wfs[ifo][1],
                                self.h00_sparse[ifo],
                                sdat['a0'], sdat['a1'],
                                sdat['b0'], sdat['b1'])
            snrs[ifo] = TimeSeries(snr, delta_t=delta_t,
                                   epoch=self.tstart[ifo] - delta_t * 2.0)
        return snrs

    @catch_waveform_error
    def _loglr(self):
        r"""Computes the log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # get model params
        p = self.current_params
        wfs = self.get_waveforms(p)
        lik = self.likelihood_function
        norm = 0.0
        filt = 0j
        pol_phase = numpy.exp(-2.0j * p['polarization'])

        self.snr_draw(wfs)
        p = self.current_params

        for ifo in self.data:
            freqs = self.fedges[ifo]
            sdat = self.sdat[ifo]
            h00 = self.h00_sparse[ifo]
            end_time = self.end_time[ifo]
            times = self.antenna_time[ifo]

            hp, hc = wfs[ifo]
            det = self.det[ifo]
            if self.precalc_antenna_factors and not self.earth_rotation:
                # The sky draw already evaluated the response and the delay
                # at every one of its pointings, so re-evaluating them here
                # for the sampled subset is pure repeated work. Only usable
                # when a sky draw actually happened (the factors are False
                # before the first one and None whenever the marginalization
                # was skipped) and when the response is a single number per
                # detector; with earth rotation on it is a function of
                # frequency and the stored single-time values do not apply.
                fp, fc, times = self.get_precalc_antenna_factors(ifo)
            else:
                fp, fc = det.antenna_pattern(p["ra"], p["dec"],
                                             0, times)
                times = det.time_delay_from_earth_center(p["ra"], p["dec"],
                                                         times)
            dtc = p["tc"] - end_time - self.ta[ifo]

            if self.lformat == 'earth_time_pol':
                filter_i, norm_i = lik(
                                       freqs, fp, fc, times, dtc, pol_phase,
                                       hp, hc, h00,
                                       sdat['a0'], sdat['a1'],
                                       sdat['b0'], sdat['b1'])
            else:
                f = (fp + 1.0j * fc) * pol_phase
                fp = f.real.copy()
                fc = f.imag.copy()
                if self.lformat == 'earth_time':
                    filter_i, norm_i = lik(
                                           freqs, fp, fc, times, dtc,
                                           hp, hc, h00,
                                           sdat['a0'], sdat['a1'],
                                           sdat['b0'], sdat['b1'])
                else:
                    filter_i, norm_i = lik(freqs, fp, fc, times + dtc,
                                           hp, hc, h00,
                                           sdat['a0'], sdat['a1'],
                                           sdat['b0'], sdat['b1'])
            filt += filter_i
            norm += norm_i
        loglr = self.marginalize_loglr(filt, norm)
        return loglr

    def _nowaveform_handler(self):
        """Sets loglr values if no waveform generated.
        """
        return -numpy.inf


class RelativeTimeDom(RelativeTime):
    """ Heterodyne likelihood optimized for time marginalization and only
    dominant-mode waveforms. This enables the ability to do inclination
    marginalization in addition to the other forms supportedy by RelativeTime.
    """
    name = "relative_time_dom"

    def get_snr(self, wfs):
        """ Return hp/hc maximized SNR time series
        """
        delta_t = 1.0 / self.sample_rate
        snrs = {}
        self.sh = {}
        self.hh = {}
        # once the marginalization points are precalculated, the
        # per-evaluation draw takes them from that pool and never looks at
        # this signal-to-noise series, so building it every call is
        # wasted; self.sh and self.hh are still needed for the likelihood.
        # Before the pool exists (setup, and the non-precalculated path)
        # the series is drawn from and must be built.
        need_series = not hasattr(self, 'premarg')
        for ifo in wfs:
            sdat = self.sdat[ifo]
            dtc = self.tstart[ifo] - self.end_time[ifo] - self.ta[ifo]

            sh, hh = snr_predictor_dom(self.fedges[ifo],
                                       dtc - delta_t * 2.0, delta_t,
                                       self.num_samples[ifo] + 4,
                                       wfs[ifo][0],
                                       self.h00_sparse[ifo],
                                       sdat['a0'], sdat['a1'],
                                       sdat['b0'], sdat['b1'])
            self.sh[ifo] = TimeSeries(sh, delta_t=delta_t,
                                      epoch=self.tstart[ifo] - delta_t * 2.0)
            self.hh[ifo] = hh
            if need_series:
                snrs[ifo] = TimeSeries(abs(sh[2:-2]) / hh ** 0.5,
                                       delta_t=delta_t,
                                       epoch=self.tstart[ifo])

        return snrs

    @catch_waveform_error
    def _loglr(self):
        r"""Computes the log likelihood ratio,
        or inner product <s|h> and <h|h> if `self.return_sh_hh` is True.

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        or
        tuple
            The inner product (<s|h>, <h|h>).
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params

        p2 = p.copy()
        p2.pop('inclination')
        wfs = self.get_waveforms(p2)

        sh_total = hh_total = 0
        ic = numpy.cos(p['inclination'])
        ip = 0.5 * (1.0 + ic * ic)
        pol_phase = numpy.exp(-2.0j * p['polarization'])

        snrs = self.get_snr(wfs)
        self.snr_draw(snrs=snrs)

        for ifo in self.sh:
            if self.precalc_antenna_factors:
                fp, fc, dt = self.get_precalc_antenna_factors(ifo)
            else:
                dt = self.det[ifo].time_delay_from_earth_center(p['ra'],
                                                                p['dec'],
                                                                p['tc'])
                fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                       0, p['tc'])
            dts = p['tc'] + dt
            f = (fp + 1.0j * fc) * pol_phase
            # Note, this includes complex conjugation already
            # as our stored inner products were hp* x data
            htf = (f.real * ip + 1.0j * f.imag * ic)
            sh = self.sh[ifo].at_time(dts, 
                                      interpolate='quadratic',
                                      extrapolate=0.0j)
            sh_total += sh * htf
            hh_total += self.hh[ifo] * abs(htf) ** 2.0

        loglr = self.marginalize_loglr(sh_total, hh_total)
        if self.return_sh_hh:
            results = (sh_total, hh_total)
        else:
            results = loglr
        return results

    def _nowaveform_handler(self):
        """Sets loglr values if no waveform generated.
        """
        loglr = sh_total = hh_total = -numpy.inf
        if self.return_sh_hh:
            results = (sh_total, hh_total)
        else:
            results = loglr
        return results
