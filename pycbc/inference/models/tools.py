""" Common utility functions for calculation of likelihoods
"""

import logging
import warnings
from distutils.util import strtobool

import numpy
import numpy.random
import tqdm

from scipy.special import logsumexp, i0e
from scipy.interpolate import RectBivariateSpline, interp1d
from pycbc.distributions import JointDistribution

from pycbc.constants import C_SI
from pycbc.detector import Detector


# Earth radius in seconds
EARTH_RADIUS = 0.031

# Metres per second, for the analytic delay -> sky inversion


def str_to_tuple(sval, ftype):
    """ Convenience parsing to convert str to tuple"""
    if sval is None:
        return ()
    return tuple(ftype(x.strip(' ')) for x in sval.split(','))


def str_to_bool(sval):
    """ Ensure value is a bool if it can be converted """
    if isinstance(sval, str):
        return strtobool(sval)
    return sval


def draw_sample(loglr, size=None):
    """ Draw a random index from a 1-d vector with loglr weights
    """
    if size:
        x = numpy.random.uniform(size=size)
    else:
        x = numpy.random.uniform()
    loglr = loglr - loglr.max()
    cdf = numpy.exp(loglr).cumsum()
    cdf /= cdf[-1]
    xl = numpy.searchsorted(cdf, x)
    return xl


class DistMarg():
    """Help class to add bookkeeping for likelihood marginalization"""

    def setup_marginalization(self,
                              variable_params,
                              marginalize_phase=False,
                              marginalize_distance=False,
                              marginalize_distance_param='distance',
                              marginalize_distance_samples=int(1e4),
                              marginalize_distance_interpolator=False,
                              marginalize_distance_snr_range=None,
                              marginalize_distance_density=None,
                              marginalize_vector_params=None,
                              marginalize_vector_samples=1e3,
                              marginalize_sky_initial_samples=1e6,
                              marginalize_sky_analytic=False,
                              **kwargs):
        """ Setup the model for use with distance marginalization

        This function sets up precalculations for distance / phase
        marginalization. For distance margininalization it modifies the
        model to internally remove distance as a parameter.

        Parameters
        ----------
        variable_params: list of strings
            The set of variable parameters
        marginalize_phase: bool, False
            Do analytic marginalization (appopriate only for 22 mode waveforms)
        marginalize_distance: bool, False
            Marginalize over distance
        marginalize_distance_param: str
            Name of the parameter that is used to determine the distance.
            This might be 'distance' or a parameter which can be converted
            to distance by a provided univariate transformation.
        marginalize_distance_interpolator: bool
            Use a pre-calculated interpolating function for the distance
            marginalized likelihood.
        marginalize_distance_snr_range: tuple of floats, (1, 50)
            The SNR range for the interpolating function to be defined in.
            If a sampler goes outside this range, the logl will be returned
            as -numpy.inf.
        marginalize_distance_density: tuple of intes, (1000, 1000)
            The dimensions of the interpolation grid over (sh, hh).

        Returns
        -------
        variable_params: list of strings
            Set of variable params (missing distance-related parameter).
        kwags: dict
            The keyword arguments to the model initialization, may be modified
            from the original set by this function.
        """
        def pop_prior(param):
            variable_params.remove(param)
            old_prior = kwargs['prior']

            dists = [d for d in old_prior.distributions
                     if param not in d.params]
            dprior = [d for d in old_prior.distributions
                      if param in d.params][0]
            prior = JointDistribution(variable_params,
                                      *dists, **old_prior.kwargs)
            kwargs['prior'] = prior
            return dprior

        self.reconstruct_phase = False
        self.reconstruct_distance = False
        self.reconstruct_vector = False
        self.precalc_antenna_factors = False

        # Handle any requested parameter vector / brute force marginalizations
        self.marginalize_vector_params = {}
        self.marginalized_vector_priors = {}
        self.vsamples = int(marginalize_vector_samples)

        self.marginalize_sky_initial_samples = \
            int(float(marginalize_sky_initial_samples))
        # PCG64 rather than the legacy global state: it is faster per
        # draw, and seeding it from that state leaves numpy.random.seed
        # pinning a run as before
        self._sky_rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 63))
        self.marginalize_sky_analytic = \
            str_to_bool(marginalize_sky_analytic)
        self._sky_detectors = None

        for param in str_to_tuple(marginalize_vector_params, str):
            logging.info('Marginalizing over %s, %s points from prior',
                         param, self.vsamples)
            self.marginalized_vector_priors[param] = pop_prior(param)

        # Remove in the future, backwards compatibility
        if 'polarization_samples' in kwargs:
            warnings.warn("use marginalize_vector_samples rather "
                          "than 'polarization_samples'", DeprecationWarning)
            pol_uniform = numpy.linspace(0, numpy.pi * 2.0, self.vsamples)
            self.marginalize_vector_params['polarization'] = pol_uniform
            self.vsamples = int(kwargs['polarization_samples'])
            kwargs.pop('polarization_samples')

        self.reset_vector_params()

        self.marginalize_phase = str_to_bool(marginalize_phase)

        self.distance_marginalization = False
        self.distance_interpolator = None

        marginalize_distance = str_to_bool(marginalize_distance)
        self.marginalize_distance = marginalize_distance
        if not marginalize_distance:
            return variable_params, kwargs

        if isinstance(marginalize_distance_snr_range, str):
            marginalize_distance_snr_range = \
                str_to_tuple(marginalize_distance_snr_range, float)

        if isinstance(marginalize_distance_density, str):
            marginalize_distance_density = \
                str_to_tuple(marginalize_distance_density, int)

        logging.info('Marginalizing over distance')

        # Take distance out of the variable params since we'll handle it
        # manually now
        dprior = pop_prior(marginalize_distance_param)

        if len(dprior.params) != 1 or not hasattr(dprior, 'bounds'):
            raise ValueError('Distance Marginalization requires a '
                             'univariate and bounded prior')

        # Set up distance prior vector and samples

        # (1) prior is using distance
        if dprior.params[0] == 'distance':
            logging.info("Prior is directly on distance, setting up "
                         "%s grid weights", marginalize_distance_samples)
            dmin, dmax = dprior.bounds['distance']
            dist_locs = numpy.linspace(dmin, dmax,
                                       int(marginalize_distance_samples))
            dist_weights = [dprior.pdf(distance=l) for l in dist_locs]
            dist_weights = numpy.array(dist_weights)

        # (2) prior is univariate and can be converted to distance
        elif marginalize_distance_param != 'distance':
            waveform_transforms = kwargs['waveform_transforms']
            pname = dprior.params[0]
            logging.info("Settings up transform,  prior is in terms of"
                         " %s", pname)
            wtrans = [d for d in waveform_transforms
                      if 'distance' not in d.outputs]
            if len(wtrans) == 0:
                wtrans = None
            kwargs['waveform_transforms'] = wtrans
            dtrans = [d for d in waveform_transforms
                      if 'distance' in d.outputs][0]
            v = dprior.rvs(int(1e8))
            d = dtrans.transform({pname: v[pname]})['distance']
            d.sort()
            cdf = numpy.arange(1, len(d)+1) / len(d)
            i = interp1d(d, cdf)
            dmin, dmax = d.min(), d.max()
            logging.info('Distance range %s-%s', dmin, dmax)
            x = numpy.linspace(dmin, dmax,
                               int(marginalize_distance_samples) + 1)
            xl, xr = x[:-1], x[1:]
            dist_locs = 0.5 * (xr + xl)
            dist_weights = i(xr) - i(xl)
        else:
            raise ValueError("No prior seems to determine the distance")

        dist_weights /= dist_weights.sum()
        dist_ref = 0.5 * (dmax + dmin)
        self.dist_locs = dist_locs
        self.distance_marginalization = dist_ref / dist_locs, dist_weights
        self.distance_interpolator = None

        if str_to_bool(marginalize_distance_interpolator):
            setup_args = {}
            if marginalize_distance_snr_range:
                setup_args['snr_range'] = marginalize_distance_snr_range
            if marginalize_distance_density:
                setup_args['density'] = marginalize_distance_density
            i = setup_distance_marg_interpolant(self.distance_marginalization,
                                                phase=self.marginalize_phase,
                                                **setup_args)
            self.distance_interpolator = i

        kwargs['static_params']['distance'] = dist_ref

        # Save marginalized parameters' name into one place,
        # coa_phase will be a static param if been marginalized
        if marginalize_distance:
            self.marginalized_params_name =\
                list(self.marginalize_vector_params.keys()) +\
                [marginalize_distance_param]

        return variable_params, kwargs

    def reset_vector_params(self):
        """ Redraw vector params from their priors
        """
        for param in self.marginalized_vector_priors:
            vprior = self.marginalized_vector_priors[param]
            values = vprior.rvs(self.vsamples)[param]
            self.marginalize_vector_params[param] = values

    def marginalize_loglr(self, sh_total, hh_total,
                          skip_vector=False, return_peak=False):
        """ Return the marginal likelihood

        Parameters
        -----------
        sh_total: float or ndarray
            The total <s|h> inner product summed over detectors
        hh_total: float or ndarray
            The total <h|h> inner product summed over detectors
        skip_vector: bool, False
            If true, and input is a vector, do not marginalize over that
            vector, instead return the likelihood values as a vector.
        """
        interpolator = self.distance_interpolator
        return_complex = False
        distance = self.distance_marginalization

        if self.reconstruct_vector:
            skip_vector = True

        if self.reconstruct_distance:
            interpolator = None
            skip_vector = True

        if self.reconstruct_phase:
            interpolator = None
            distance = False
            skip_vector = True
            return_complex = True

        return marginalize_likelihood(sh_total, hh_total,
                                      logw=self.marginalize_vector_weights,
                                      phase=self.marginalize_phase,
                                      interpolator=interpolator,
                                      distance=distance,
                                      skip_vector=skip_vector,
                                      return_complex=return_complex,
                                      return_peak=return_peak)

    def premarg_draw(self):
        """ Choose random samples from prechosen set"""

        # Update the current proposed times and the marginalization values
        logw = self.premarg['logw_partial']
        if self.vsamples == len(logw):
            choice = slice(None, None)
        else:
            choice = numpy.random.choice(len(logw), size=self.vsamples,
                                         replace=False)

        for k in self.snr_params:
            self.marginalize_vector_params[k] = self.premarg[k][choice]

        self._current_params.update(self.marginalize_vector_params)
        self.sample_idx = self.premarg['sample_idx'][choice]

        # Update the importance weights for each vector sample
        logw = self.marginalize_vector_weights + logw[choice]
        self.marginalize_vector_weights = logw - logsumexp(logw)
        return self.marginalize_vector_params

    def snr_draw(self, wfs=None, snrs=None, size=None):
        """ Improve the monte-carlo vector marginalization using the SNR time
        series of each detector
        """
        try:
            p = self.current_params
            set_scalar = numpy.isscalar(p['tc'])
        except:
            set_scalar = False

        if not set_scalar:
            if hasattr(self, 'premarg'):
                return self.premarg_draw()

            if snrs is None:
                snrs = self.get_snr(wfs)
            if ('tc' in self.marginalized_vector_priors and
                not ('ra' in self.marginalized_vector_priors
                     or 'dec' in self.marginalized_vector_priors)):
                return self.draw_times(snrs, size=size)
            elif ('tc' in self.marginalized_vector_priors and
                  'ra' in self.marginalized_vector_priors and
                  'dec' in self.marginalized_vector_priors):
                return self.draw_sky_times(snrs, size=size)
        else:
            # OK, we couldn't do anything with the requested monte-carlo
            # marginalizations.
            self.precalc_antenna_factors = None
            return None

    def draw_times(self, snrs, size=None):
        """ Draw times consistent with the incoherent network SNR

        Parameters
        ----------
        snrs: dist of TimeSeries
        """
        if not hasattr(self, 'tinfo'):
            # determine the rough time offsets for this sky location
            tcprior = self.marginalized_vector_priors['tc']
            tcmin, tcmax = tcprior.bounds['tc']
            tcave = (tcmax + tcmin) / 2.0
            ifos = list(snrs.keys())
            if hasattr(self, 'keep_ifos'):
                ifos = self.keep_ifos
            d = {ifo: Detector(ifo, reference_time=tcave) for ifo in ifos}
            self.tinfo = tcmin, tcmax, tcave, ifos, d
            self.snr_params = ['tc']

        tcmin, tcmax, tcave, ifos, d = self.tinfo
        vsamples = size if size is not None else self.vsamples

        # Determine the weights for the valid time range
        ra = self._current_params['ra']
        dec = self._current_params['dec']

        # Determine the common valid time range
        iref = ifos[0]
        dref = d[iref]
        dt = dref.time_delay_from_earth_center(ra, dec, tcave)

        starts = []
        ends = []

        delt = snrs[iref].delta_t
        tmin = tcmin + dt - delt
        tmax = tcmax + dt + delt
        if hasattr(self, 'tstart'):
            tmin = self.tstart[iref]
            tmax = self.tend[iref]

        # Make sure we draw from times within prior and that have enough
        # SNR calculated to do later interpolation
        starts.append(max(tmin, snrs[iref].start_time + delt))
        ends.append(min(tmax, snrs[iref].end_time - delt * 2))

        idels = {}
        for ifo in ifos[1:]:
            dti = d[ifo].time_delay_from_detector(dref, ra, dec, tcave)
            idel = round(dti / snrs[iref].delta_t) * snrs[iref].delta_t
            idels[ifo] = idel

            starts.append(snrs[ifo].start_time - idel)
            ends.append(snrs[ifo].end_time - idel)
        start = max(starts)
        end = min(ends)
        if end <= start:
            return

        # get the weights
        snr = snrs[iref].time_slice(start, end, mode='nearest')
        logweight = snr.squared_norm().numpy()
        for ifo in ifos[1:]:
            idel = idels[ifo]
            snrv = snrs[ifo].time_slice(snr.start_time + idel,
                                        snr.end_time + idel,
                                        mode='nearest')
            logweight += snrv.squared_norm().numpy()
        logweight /= 2.0
        logweight -= logsumexp(logweight)  # Normalize to PDF

        # Draw proportional to the incoherent likelihood
        # Draw first which time sample
        tci = draw_sample(logweight, size=vsamples)
        # Second draw a subsample size offset so that all times are covered
        tct = numpy.random.uniform(-snr.delta_t / 2.0,
                                   snr.delta_t / 2.0,
                                   size=vsamples)
        tc = tct + tci * snr.delta_t + float(snr.start_time) - dt

        # Update the current proposed times and the marginalization values
        # assumes uniform prior!
        logw = - logweight[tci] + numpy.log(1.0 / len(logweight))
        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['logw_partial'] = logw

        if self._current_params is not None:
            # Update the importance weights for each vector sample
            self._current_params.update(self.marginalize_vector_params)
            self.marginalize_vector_weights += logw

        return self.marginalize_vector_params

    @staticmethod
    def _draw_in_cells(logw, lo, hi, rng):
        """Draw one cell per sample from ``[lo, hi)`` of a shared law.

        Every sample has its own bounds, because the window each one may
        draw in follows from where the previous detector was drawn. The
        cumulative sum is over the whole series and is formed once; a
        sample's own normalisation is then the difference of two of its
        entries.

        Returns the cells and the log probability of drawing each within
        its own window, or ``-inf`` where the window holds no cell. The
        cell returned for such a sample is arbitrary; every caller drops
        it on the ``-inf``.
        """
        lmax = logw.max()
        cum = numpy.concatenate(([0.0], numpy.cumsum(numpy.exp(logw - lmax))))
        # A window holding no cell is rejected on the last line; what
        # these two do is keep the arithmetic between here and there
        # defined. Left alone, such a sample divides by a zero or negative
        # norm and carries an inf or a nan through to a value that is
        # discarded anyway, so the only thing lost would be two warnings
        # -- but the same inf from a window that is NOT empty would be a
        # real fault, and silencing the warning would hide it. So they are
        # pointed at the whole series, which is always drawable, and the
        # result thrown away.
        empty = hi <= lo
        lo = numpy.where(empty, 0, lo)
        hi = numpy.where(empty, len(logw), hi)
        norm = cum[hi] - cum[lo]
        target = cum[lo] + rng.random(len(lo)) * norm
        idx = numpy.clip(numpy.searchsorted(cum, target) - 1,
                         0, len(logw) - 1)
        mass = (logw[idx] - lmax) - numpy.log(norm)
        return idx, numpy.where(empty, -numpy.inf, mass)

    @staticmethod
    def _weighted_draw(weights, shape, rng):
        """Draw indices ``propto weights``, with replacement.

        How many draws each cell wins is exactly a multinomial, so the
        counts can be drawn in one call and expanded, which costs the
        length of the series plus the number of draws and searches for
        nothing. The expansion comes out ordered by cell, and the caller
        pairs each draw with a different sample, so it is shuffled back
        into an independent order.
        """
        size = int(numpy.prod(shape))
        counts = rng.multinomial(size, weights / weights.sum())
        idx = numpy.repeat(numpy.arange(len(weights)), counts)
        rng.shuffle(idx)
        return idx.reshape(shape)

    @staticmethod
    def _ellipse_slice(smat, dt12):
        """The dt13 the physical sky allows, once dt12 is fixed.

        The physical region is the ellipse ``dt^T S dt <= 1``. Fixing dt12
        cuts a segment from it whose ends are the roots of a quadratic in
        dt13. Returns the segment's centre and half width, and which
        samples have no segment because their dt12 is already outside.
        """
        quad = smat[1, 1]
        lin = 2.0 * smat[0, 1] * dt12
        disc = lin * lin - 4.0 * quad * (smat[0, 0] * dt12 * dt12 - 1.0)
        root = numpy.sqrt(numpy.maximum(disc, 0.0))
        return -lin / (2.0 * quad), root / (2.0 * quad), disc <= 0.0

    def _draw_second_delay(self, third, smat, epoch,
                           t_off, dt12):
        """Draw dt13 on the ellipse slice the first delay allows.

        On that slice the exact isotropic sky prior is the ARCSINE law,
        and substituting ``dt13 = mid + half*sin(theta)`` makes it uniform
        in theta -- which is what cancels the ``1/s`` Jacobian singularity
        at the ends of the slice, where the source lies in the detector
        plane and the two timing cones go tangent. A proposal uniform in
        dt13 instead has log-divergent weight variance and is worse than
        the map this replaces.

        The third detector's own SNR then tilts the draw, exactly as the
        second detector's does: the slice is a window of arrival times, a
        cell is drawn from the series inside it, and dt13 is placed within
        that cell by the arcsine. The window is the cells that OVERLAP the
        slice rather than those centred in it, because the arcsine density
        diverges at the ends and the two edge cells carry a fifth of the
        mass.

        The weight is then ``W P_j / sh_j`` for the drawn cell j, where W
        is the window's own share of the series. That is exact for the
        cell drawn, not an estimate, so nothing here has a sample count to
        choose.
        """
        mid, half, invalid = self._ellipse_slice(smat, dt12)
        # a sample whose slice has no width is already invalid; this only
        # keeps the divisions below finite until it is dropped
        half = numpy.maximum(half, 1e-300)

        snr, logl = third
        base, delta = float(snr.start_time - epoch), float(snr.delta_t)
        ncell = len(logl)
        # the slice, as a window of this detector's arrival times
        lo = numpy.clip(numpy.floor((t_off - mid - half - base) / delta),
                        0, ncell).astype(int)
        hi = numpy.clip(numpy.ceil((t_off - mid + half - base) / delta) + 1,
                        0, ncell).astype(int)
        cell, logq = self._draw_in_cells(logl, lo, hi, self._sky_rng)

        # where the drawn cell sits on the slice, and half a cell, both in
        # units of the half width, which is what the arcsine takes
        centre = (t_off - (base + cell * delta) - mid) / half
        step = delta / (2.0 * half)
        edge_lo = numpy.arcsin(numpy.clip(centre - step, -1.0, 1.0))
        edge_hi = numpy.arcsin(numpy.clip(centre + step, -1.0, 1.0))
        # arcsine is monotone and the step is positive, so hi >= lo always
        mass = (edge_hi - edge_lo) / numpy.pi
        invalid |= ~(mass > 0) | ~numpy.isfinite(logq)
        safe = numpy.maximum(mass, 1e-300)

        theta = (edge_lo
                 + self._sky_rng.random(len(dt12)) * (edge_hi - edge_lo))
        dt13 = mid + half * numpy.sin(theta)
        # q = (sh_j / W) p(dt13) / P_j  =>  p/q = W P_j / sh_j
        dens = logq - numpy.log(safe)
        return dt13, dens, invalid

    def analytic_sky_draw(self, snrs, ifos, vsamples):
        """Draw (tc, ra, dec) by inverting the inter-detector delays exactly.

        An alternative to the pregenerated delay-to-sky map. Each delay is one
        linear constraint on the unit source direction,

            nhat . d_ij = -c * dt_ij       d_ij = r_i - r_j

        so ``min(ndet - 1, 2)`` delays -- plus one azimuth when there is only a
        single delay -- determine the sky. The physical region is exactly the
        ellipse ``dt^T S dt <= 1``, so a drawn delay tuple cannot fail to
        correspond to a sky position, and the Jacobian is analytic, so the
        prior/proposal ratio is closed form rather than approximate.

        The isotropic sky prior in delay space is exactly uniform in the first
        delay and arcsine in the second, and drawing them that way makes every
        geometric factor collapse to the constant ``-log(2 t12max)``. The
        residual weight is then only the SNR tilt on the drawn times.

        Returns True on success. On failure the caller falls back to the map,
        which is reported so that a run cannot silently change estimator.
        """
        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        tcmin, tcmax = float(tcmin), float(tcmax)
        if self._sky_detectors is None:
            # the only part of this worth keeping between calls: Detector
            # scans the LAL table and builds an astropy EarthLocation, about
            # a millisecond for three of them, against microseconds for
            # everything else here
            self._sky_detectors = {
                i: Detector(i, reference_time=0.5 * (tcmin + tcmax))
                for i in self.data}
        dets = self._sky_detectors
        loc = {i: numpy.asarray(d.location) for i, d in dets.items()}
        gmst = next(iter(dets.values())).gmst_estimate(
            0.5 * (tcmin + tcmax))

        # each detector's SNR series, cut to the times the draw may use
        series = {}
        for ifo in ifos:
            snr = snrs[ifo]
            lo_t, hi_t = tcmin - EARTH_RADIUS, tcmax + EARTH_RADIUS
            if hasattr(self, 'tstart'):
                lo_t, hi_t = self.tstart[ifo], self.tend[ifo]
            snr = snr.time_slice(max(lo_t, snr.start_time + snr.delta_t),
                                 min(hi_t, snr.end_time - snr.delta_t * 2),
                                 mode='nearest')
            series[ifo] = (snr, snr.squared_norm().numpy() / 2.0)

        order = sorted(ifos, key=lambda i: -series[i][1].max())
        ndet = len(order)
        log_tcspan = numpy.log(tcmax - tcmin)

        snr0, logl0 = series[order[0]]
        delta0 = float(snr0.delta_t)
        l0max = logl0.max()
        shifted0 = numpy.exp(logl0 - l0max)
        i0 = self._weighted_draw(shifted0, vsamples, self._sky_rng)
        # Every time below is an OFFSET from this epoch, never an absolute
        # GPS. The delays are differences of times near 1e9 s, which a
        # double resolves to 2.4e-7 s -- a thousandth of a sample -- and
        # that error lands on the cell boundaries the weights are read
        # from, putting about 0.1% of samples in the neighbouring cell.
        # As offsets the differences are exact to ~1e-18.
        epoch = snr0.start_time
        epoch_t = float(epoch)
        t_off = (i0 * delta0
                 + self._sky_rng.uniform(-delta0 / 2.0, delta0 / 2.0,
                                         vsamples))
        dens = ((logl0[i0] - l0max) - numpy.log(shifted0.sum())
                - numpy.log(delta0))
        invalid = numpy.zeros(vsamples, dtype=bool)

        if ndet == 1:
            # no delay to invert, so the sky comes from the prior and only
            # the arrival time carries anything
            ra = self.marginalized_vector_priors['ra'].rvs(size=vsamples)['ra']
            dec = self.marginalized_vector_priors['dec'].rvs(
                size=vsamples)['dec']
            lon = ra - gmst
            nhat = numpy.array([numpy.cos(dec) * numpy.cos(lon),
                                numpy.cos(dec) * numpy.sin(lon),
                                numpy.sin(dec)])
            logw = -log_tcspan - dens
        else:
            # each delay is one linear constraint on the direction; the
            # pseudo-inverse of the baselines turns a delay pair back into
            # the component of the direction they fix
            baselines = [loc[order[0]] - loc[i] for i in order[1:3]]
            mpinv = numpy.linalg.pinv(numpy.vstack(baselines))
            t12max = dets[order[0]].light_travel_time_to_detector(
                dets[order[1]])

            snr1, logl1 = series[order[1]]
            base1 = float(snr1.start_time - epoch)
            delta1 = float(snr1.delta_t)
            # the second detector's arrival lies within one light crossing
            # of the first, intersected with the series we have: under peak
            # locking the region is narrower than the crossing, so this
            # clip is the usual case rather than an edge one. A cell is in
            # the window if it overlaps it at all, as in the third
            # detector below: a drawn cell is jittered across its width, so
            # excluding one that straddles an end would leave the times in
            # the overlap undrawable rather than merely unlikely.
            lo1 = numpy.clip(numpy.floor((t_off - t12max - base1) / delta1),
                             0, len(logl1)).astype(int)
            hi1 = numpy.clip(numpy.ceil((t_off + t12max - base1) / delta1)
                             + 1, 0, len(logl1)).astype(int)
            i1, mass1 = self._draw_in_cells(logl1, lo1, hi1, self._sky_rng)
            t_two = (base1 + i1 * delta1
                     + self._sky_rng.uniform(-delta1 / 2.0, delta1 / 2.0,
                                             vsamples))
            dens1 = mass1 - numpy.log(delta1)
            dens += dens1
            dt12 = t_off - t_two
            invalid |= ~numpy.isfinite(dens1) | (numpy.abs(dt12) > t12max)

            branch_corr = 0.0
            if ndet == 2:
                # One delay leaves the azimuth about the baseline free.
                # The Jacobian is constant there, so sweeping it uniformly
                # is exactly isotropic and needs no correction.
                dhat = baselines[0] / numpy.linalg.norm(baselines[0])
                uvec = numpy.cross(dhat, [0.0, 0.0, 1.0])
                uvec = uvec / numpy.linalg.norm(uvec)
                vvec = numpy.cross(dhat, uvec)
                cos_t = numpy.clip(-dt12 / t12max, -1.0, 1.0)
                sin_t = numpy.sqrt(numpy.maximum(1.0 - cos_t * cos_t, 0.0))
                az = self._sky_rng.uniform(0, 2 * numpy.pi, vsamples)
                nhat = (cos_t * dhat[:, None]
                        + sin_t * (numpy.cos(az) * uvec[:, None]
                                   + numpy.sin(az) * vvec[:, None]))
            else:
                normal = numpy.cross(*baselines)
                e3 = normal / numpy.linalg.norm(normal)
                dt13, dens2, bad2 = self._draw_second_delay(
                    series[order[2]], C_SI ** 2.0 * (mpinv.T @ mpinv),
                    epoch, t_off, dt12)
                dens += dens2
                invalid |= bad2
                npar = mpinv @ (-C_SI * numpy.vstack([dt12, dt13]))
                perp = numpy.sqrt(numpy.clip(1.0 - (npar * npar).sum(0),
                                             0.0, None))
                # The delays fix the direction only up to a reflection in
                # the detector plane. Both roots are physical but put the
                # coalescence at times differing by twice the reference
                # detector's offset from that plane, so for some geometries
                # one of them is outside the tc prior for every sample.
                # Testing a root is one dot product, so rather than draw
                # blind and reject, draw only among the roots that survive
                # and scale the weight by how many did: with both kept the
                # expectation is (w+ + w-)/2 either way.
                loc0 = loc[order[0]] / C_SI
                tc_mid = epoch_t + t_off + loc0 @ npar
                swing = perp * float(loc0 @ e3)
                ok_up = ((tc_mid + swing >= tcmin)
                         & (tc_mid + swing <= tcmax))
                ok_dn = ((tc_mid - swing >= tcmin)
                         & (tc_mid - swing <= tcmax))
                invalid |= ~(ok_up | ok_dn)
                up = ok_up & (~ok_dn | (self._sky_rng.random(vsamples) < 0.5))
                nhat = npar + numpy.where(up, perp, -perp) * e3[:, None]
                # both kept, the half cancels the two roots; one kept, it
                # does not and the survivor carries the pair
                branch_corr = numpy.where(ok_up & ok_dn, 0.0, -numpy.log(2.0))
            logw = (-log_tcspan - numpy.log(2.0 * t12max)
                    - dens + branch_corr)

        fplus, fcross, delay = {}, {}, {}
        for ifo, det in dets.items():
            fplus[ifo], fcross[ifo] = det.antenna_pattern_from_direction(nhat)
            delay[ifo] = det.time_delay_from_direction(nhat)
        tc = epoch_t + t_off - delay[order[0]]
        outside = (tc < tcmin) | (tc > tcmax)
        logw = numpy.where(invalid | outside, -numpy.inf, logw)
        if not numpy.isfinite(logw).any():
            return False

        self.snr_params = ['tc', 'ra', 'dec']
        self.sample_idx = numpy.arange(vsamples)
        self.precalc_antenna_factors = fplus, fcross, delay
        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['dec'] = numpy.arcsin(
            numpy.clip(nhat[2], -1.0, 1.0))
        self.marginalize_vector_params['ra'] = numpy.mod(
            numpy.arctan2(nhat[1], nhat[0]) + gmst, 2 * numpy.pi)
        self.marginalize_vector_params['logw_partial'] = logw
        if self._current_params is not None:
            self._current_params.update(self.marginalize_vector_params)
            self.marginalize_vector_weights += logw
        return True

    def draw_sky_times(self, snrs, size=None):
        """ Draw ra, dec, and tc together using SNR timeseries to determine
        monte-carlo weights.
        """
        # First setup
        # precalculate dense sky grid and make dict and or array of the results
        ifos = list(snrs.keys())
        if hasattr(self, 'keep_ifos'):
            ifos = self.keep_ifos
        ikey = ''.join(ifos)

        vsamples = size if size is not None else self.vsamples

        if self.marginalize_sky_analytic and len(ifos) > 0:
            # alternative path: invert the delays instead of looking them up.
            # Falls through to the map below if it cannot produce a usable
            # draw, and says so, so a run cannot silently change estimator.
            if self.analytic_sky_draw(snrs, ifos, vsamples):
                return self.marginalize_vector_params
            if not getattr(self, '_analytic_fellback', False):
                self._analytic_fellback = True
                logging.warning(
                    "analytic sky draw produced no usable samples; "
                    "falling back to the pregenerated map")

        # No good SNR peaks, go with prior draw
        if len(ifos) == 0:
            self.marginalize_vector_params['logw_partial'] = numpy.zeros(vsamples)
            return

        def make_init():
            self.snr_params = ['tc', 'ra', 'dec']
            size = self.marginalize_sky_initial_samples
            logging.info('drawing samples: %s', size)
            ra = self.marginalized_vector_priors['ra'].rvs(size=size)['ra']
            dec = self.marginalized_vector_priors['dec'].rvs(size=size)['dec']
            tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
            tcave = (tcmax + tcmin) / 2.0
            d = {ifo: Detector(ifo, reference_time=tcave) for ifo in self.data}

            # What data structure to hold times? Dict of offset -> list?
            logging.info('sorting into time delay dict')
            dts = []
            for i in range(len(ifos) - 1):
                dt = d[ifos[0]].time_delay_from_detector(d[ifos[i+1]],
                                                         ra, dec, tcave)
                dt = numpy.rint(dt / snrs[ifos[0]].delta_t)
                dts.append(dt)

            fp, fc, dtc = {}, {}, {}
            for ifo in self.data:
                fp[ifo], fc[ifo] = d[ifo].antenna_pattern(ra, dec, 0.0, tcave)
                dtc[ifo] = d[ifo].time_delay_from_earth_center(ra, dec, tcave)

            dmap = {}
            for i, t in enumerate(tqdm.tqdm(zip(*dts))):
                if t not in dmap:
                    dmap[t] = []
                dmap[t].append(i)

            if len(ifos) == 1:
                dmap[()] = numpy.arange(0, size, 1).astype(int)

            # Sky prior by bin
            bin_prior = {t: len(dmap[t]) / size for t in dmap}

            return dmap, tcmin, tcmax, fp, fc, ra, dec, dtc, bin_prior

        if not hasattr(self, 'tinfo'):
            self.tinfo = {}

        if ikey not in self.tinfo:
            logging.info('pregenerating sky pointings')
            self.tinfo[ikey] = make_init()

        dmap, tcmin, tcmax, fp, fc, ra, dec, dtc, bin_prior = self.tinfo[ikey]

        # draw times from each snr time series
        # Is it worth doing this if some detector has low SNR?
        sref = None
        iref = None
        idx = []
        dx = []
        mcweight = None
        for ifo in ifos:
            snr = snrs[ifo]
            tmin, tmax = tcmin - EARTH_RADIUS, tcmax + EARTH_RADIUS
            if hasattr(self, 'tstart'):
                tmin = self.tstart[ifo]
                tmax = self.tend[ifo]

            start = max(tmin, snr.start_time + snr.delta_t)
            end = min(tmax, snr.end_time - snr.delta_t * 2)
            snr = snr.time_slice(start, end, mode='nearest')

            w = snr.squared_norm().numpy() / 2.0
            i = draw_sample(w, size=vsamples)

            if sref is not None:
                mcweight += w[i]
                delt = float(snr.start_time - sref.start_time)
                i += round(delt / sref.delta_t)
                dx.append(iref - i)
            else:
                sref = snr
                iref = i
                mcweight = w[i]

            idx.append(i)
        mcweight -= logsumexp(mcweight)

        # check if delay is in dict, if not, throw out
        ti = []
        ix = []
        wi = []
        rand = numpy.random.uniform(0, 1, size=vsamples)
        for i in range(vsamples):
            t = tuple(x[i] for x in dx)
            if t in dmap:
                randi = int(rand[i] * (len(dmap[t])))
                ix.append(dmap[t][randi])
                wi.append(bin_prior[t])
                ti.append(i)

        # If we had really poor efficiency at finding a point, we should
        # give up and just use the original random draws
        if len(ix) < 0.05 * vsamples:
            self.marginalize_vector_params['logw_partial'] = numpy.zeros(vsamples)
            return

        # fill back to fixed size with repeat samples
        # sample order is random, so this should be OK statistically
        ix = numpy.resize(numpy.array(ix, dtype=int), vsamples)
        self.sample_idx = ix
        self.precalc_antenna_factors = fp, fc, dtc
        resize_factor = len(ti) / vsamples

        ra = ra[ix]
        dec = dec[ix]
        dtc = {ifo: dtc[ifo][ix] for ifo in dtc}

        ti = numpy.resize(numpy.array(ti, dtype=int), vsamples)
        wi = numpy.resize(numpy.array(wi), vsamples)

        # Second draw a subsample size offset so that all times are covered
        tct = numpy.random.uniform(-snr.delta_t / 2.0,
                                   snr.delta_t / 2.0,
                                   size=len(ti))

        tc = tct + iref[ti] * snr.delta_t + float(sref.start_time) - dtc[ifos[0]]

        # Update the current proposed times and the marginalization values
        # There's an overall normalization here which may introduce a constant
        # factor at the moment.
        logw_sky = -mcweight[ti] + numpy.log(wi) - numpy.log(resize_factor)

        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['ra'] = ra
        self.marginalize_vector_params['dec'] = dec
        self.marginalize_vector_params['logw_partial'] = logw_sky

        if self._current_params is not None:
            # Update the importance weights for each vector sample
            self._current_params.update(self.marginalize_vector_params)
            self.marginalize_vector_weights += logw_sky

        return self.marginalize_vector_params

    def get_precalc_antenna_factors(self, ifo):
        """ Get the antenna factors for marginalized samples if they exist """
        ix = self.sample_idx
        fp, fc, dtc = self.precalc_antenna_factors
        return fp[ifo][ix], fc[ifo][ix], dtc[ifo][ix]

    def setup_peak_lock(self,
                        sample_rate=4096,
                        snrs=None,
                        peak_lock_snr=None,
                        peak_lock_ratio=1e4,
                        peak_lock_region=4,
                        **kwargs):
        """ Determine where to constrain marginalization based on
        the observed reference SNR peaks.

        Parameters
        ----------
        sample_rate : float
            The SNR sample rate
        snrs : Dict of SNR time series
            Either provide this or the model needs a function
            to get the reference SNRs.
        marginalize_sky_analytic: bool
            Draw (tc, ra, dec) by inverting the inter-detector delays in closed
            form instead of looking them up in the pregenerated map. The
            physical region is the exact ellipse the delays allow, so a drawn
            delay tuple cannot fail to correspond to a sky position, and the
            prior/proposal ratio is analytic rather than approximate.
            Falls back to the map, with a warning, if it cannot produce a
            usable draw.
        peak_lock_snr: float
            The minimum SNR to bother restricting from the prior range
        peak_lock_ratio: float
            The likelihood ratio (not log) relative to the peak to
            act as a threshold bounding region.
        peak_lock_region: int
            Number of samples to inclue beyond the strict region
            determined by the relative likelihood
        """

        if 'tc' not in self.marginalized_vector_priors:
            return

        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        tstart = tcmin - EARTH_RADIUS
        tmax = tcmax - tcmin + EARTH_RADIUS * 2.0
        num_samples = int(tmax * sample_rate)
        self.tstart = {ifo: tstart for ifo in self.data}
        self.num_samples = {ifo: num_samples for ifo in self.data}

        if snrs is None:
            if not hasattr(self, 'ref_snr'):
                raise ValueError("Model didn't have a reference SNR!")
            snrs = self.ref_snr

        # Restrict the time range for constructing SNR time series
        # to identifiable peaks
        if peak_lock_snr is not None:
            peak_lock_snr = float(peak_lock_snr)
            peak_lock_ratio = float(peak_lock_ratio)
            peak_lock_region = int(peak_lock_region)

            for ifo in snrs:
                s = max(tstart, snrs[ifo].start_time)
                e = min(tstart + tmax, snrs[ifo].end_time)
                z = snrs[ifo].time_slice(s, e, mode='nearest')
                peak_snr, imax = z.abs_max_loc()
                times = z.sample_times
                peak_time = times[imax]

                logging.info('%s: Max Ref SNR Peak of %s at %s',
                             ifo, peak_snr, peak_time)

                if peak_snr > peak_lock_snr:
                    target = peak_snr ** 2.0 / 2.0 - numpy.log(peak_lock_ratio)
                    target = (target * 2.0) ** 0.5

                    region = numpy.where(abs(z) > target)[0]
                    ts = times[region[0]] - peak_lock_region / sample_rate
                    te = times[region[-1]] + peak_lock_region / sample_rate
                    self.tstart[ifo] = ts
                    self.num_samples[ifo] = int((te - ts) * sample_rate)

            # Check times are commensurate with each other
            for ifo in snrs:
                ts = self.tstart[ifo]
                te = ts + self.num_samples[ifo] / sample_rate

                for ifo2 in snrs:
                    if ifo == ifo2:
                        continue
                    ts2 = self.tstart[ifo2]
                    te2 = ts2 + self.num_samples[ifo2] / sample_rate
                    det = Detector(ifo)
                    dt = Detector(ifo2).light_travel_time_to_detector(det)

                    ts = max(ts, ts2 - dt)
                    te = min(te, te2 + dt)

                self.tstart[ifo] = ts
                self.num_samples[ifo] = int((te - ts) * sample_rate) + 1
                logging.info('%s: use region %s-%s, %s points',
                             ifo, ts, te, self.num_samples[ifo])

        self.tend = self.tstart.copy()
        for ifo in snrs:
            self.tend[ifo] += self.num_samples[ifo] / sample_rate

    def draw_ifos(self, snrs, peak_snr_threshold=4.0, log=True,
                  precalculate_marginalization_points=False,
                  **kwargs):
        """ Helper utility to determine which ifos we should use based on the
        reference SNR time series.
        """
        if 'tc' not in self.marginalized_vector_priors:
            return

        peak_snr_threshold = float(peak_snr_threshold)

        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        ifos = list(snrs.keys())
        keep_ifos = []
        psnrs = []
        for ifo in ifos:
            snr = snrs[ifo]
            start = max(tcmin - EARTH_RADIUS, snr.start_time)
            end = min(tcmax + EARTH_RADIUS, snr.end_time)
            snr = snr.time_slice(start, end, mode='nearest')
            psnr = abs(snr).max()
            if psnr > peak_snr_threshold:
                keep_ifos.append(ifo)
            psnrs.append(psnr)

        if log:
            logging.info("Ifos used for SNR based draws:"
                         " %s, snrs: %s, peak_snr_threshold=%s",
                         keep_ifos, psnrs, peak_snr_threshold)

        self.keep_ifos = keep_ifos

        if precalculate_marginalization_points:
            num_points = int(float(precalculate_marginalization_points))
            self.premarg = self.snr_draw(size=num_points, snrs=snrs).copy()
            self.premarg['sample_idx'] = self.sample_idx

        return keep_ifos

    @property
    def current_params(self):
        """ The current parameters

        If a parameter has been vector marginalized, the likelihood should
        expect an array for the given parameter. This allows transparent
        vectorization for many models.
        """
        params = self._current_params
        for k in self.marginalize_vector_params:
            if k not in params:
                params[k] = self.marginalize_vector_params[k]
        self.marginalize_vector_weights = - numpy.log(self.vsamples)
        return params

    def reconstruct(self, rec=None, seed=None, set_loglr=None):
        """ Reconstruct the distance or vectored marginalized parameter
        of this class.
        """
        if seed:
            numpy.random.seed(seed)

        if rec is None:
            rec = {}

        if set_loglr is None:
            def get_loglr():
                p = self.current_params.copy()
                p.update(rec)
                self.update(**p)
                return self.loglr
        else:
            get_loglr = set_loglr

        if self.marginalize_vector_params:
            logging.debug('Reconstruct vector')
            self.reconstruct_vector = True
            self.reset_vector_params()
            loglr = get_loglr()
            xl = draw_sample(loglr + self.marginalize_vector_weights)
            for k in self.marginalize_vector_params:
                rec[k] = self.marginalize_vector_params[k][xl]
            self.reconstruct_vector = False

        if self.distance_marginalization:
            logging.debug('Reconstruct distance')
            # call likelihood to get vector output
            self.reconstruct_distance = True
            _, weights = self.distance_marginalization
            loglr = get_loglr()
            xl = draw_sample(loglr + numpy.log(weights))
            rec['distance'] = self.dist_locs[xl]
            self.reconstruct_distance = False

        if self.marginalize_phase:
            logging.debug('Reconstruct phase')
            self.reconstruct_phase = True
            s, h = get_loglr()
            phasev = numpy.linspace(0, numpy.pi*2.0, int(1e4))
            # This assumes that the template was conjugated in inner products
            loglr = (numpy.exp(-2.0j * phasev) * s).real + h
            xl = draw_sample(loglr)
            rec['coa_phase'] = phasev[xl]
            self.reconstruct_phase = False

        rec['loglr'] = loglr[xl]
        rec['loglikelihood'] = self.lognl + rec['loglr']
        return rec


def setup_distance_marg_interpolant(dist_marg,
                                    phase=False,
                                    snr_range=(1, 50),
                                    density=(1000, 1000)):
    """ Create the interpolant for distance marginalization

    Parameters
    ----------
    dist_marg: tuple of two arrays
        The (dist_loc, dist_weight) tuple which defines the grid
        for integrating over distance
    snr_range: tuple of (float, float)
        Tuple of min, max SNR that the interpolant is expected to work
        for.
    density: tuple of (float, float)
        The number of samples in either dimension of the 2d interpolant

    Returns
    -------
    interp: function
        Function which returns the precalculated likelihood for a given
        inner product sh/hh.
    """
    dist_rescale, _ = dist_marg
    logging.info("Interpolator valid for SNRs in %s", snr_range)
    logging.info("Interpolator using grid %s", density)
    # approximate maximum shr and hhr values, assuming the true SNR is
    # within the indicated range (and neglecting noise fluctuations)
    snr_min, snr_max = snr_range
    smax = dist_rescale.max()
    smin = dist_rescale.min()
    shr_max = snr_max ** 2.0 / smin
    hhr_max = snr_max ** 2.0 / smin / smin

    shr_min = snr_min ** 2.0 / smax
    hhr_min = snr_min ** 2.0 / smax / smax

    shr = numpy.geomspace(shr_min, shr_max, density[0])
    hhr = numpy.geomspace(hhr_min, hhr_max, density[1])
    lvals = numpy.zeros((len(shr), len(hhr)))
    logging.info('Setup up likelihood interpolator')
    for i, sh in enumerate(tqdm.tqdm(shr)):
        for j, hh in enumerate(hhr):
            lvals[i, j] = marginalize_likelihood(sh, hh,
                                                 distance=dist_marg,
                                                 phase=phase)
    interp = RectBivariateSpline(shr, hhr, lvals)

    # said once, the first time it happens
    warned = [False]

    def warn_out_of_range():
        warned[0] = True
        logging.warning(
            "A likelihood evaluation asked for a signal to noise ratio "
            "outside marginalize_distance_snr_range %s; beyond it the "
            "distance marginalized likelihood is set to zero, which biases "
            "the result. Widen the range.", snr_range)

    def interp_wrapper(x, y, bounds_check=True):
        k = None
        if bounds_check:
            if isinstance(x, float):
                if x > shr_max or x < shr_min or y > hhr_max or y < hhr_min:
                    if not warned[0]:
                        warn_out_of_range()
                    return -numpy.inf
            else:
                k = (x > shr_max) | (x < shr_min)
                k = k | (y > hhr_max) | (y < hhr_min)
                # short circuits, so this costs nothing once said
                if not warned[0] and k.any():
                    warn_out_of_range()

        v = interp(x, y, grid=False)
        if k is not None:
            v[k] = -numpy.inf
        return v
    return interp_wrapper


def marginalize_likelihood(sh, hh,
                           logw=None,
                           phase=False,
                           distance=False,
                           skip_vector=False,
                           interpolator=None,
                           return_peak=False,
                           return_complex=False,
                           ):
    """ Return the marginalized likelihood.

    Apply various marginalizations to the data, including phase, distance,
    and brute-force vector marginalizations. Several options relate
    to how the distance marginalization is approximated and others allow for
    special return products to aid in parameter reconstruction.

    Parameters
    ----------
    sh: complex float or numpy.ndarray
        The data-template inner product
    hh: complex float or numpy.ndarray
        The template-template inner product
    logw:
        log weighting factors if vector marginalization is used, if not
        given, each sample is assumed to be equally weighted
    phase: bool, False
        Enable phase marginalization. Only use if orbital phase can be related
        to just a single overall phase (e.g. not true for waveform with
        sub-dominant modes)
    skip_vector: bool, False
        Don't apply marginalization of vector component of input (i.e. leave
        as vector).
    interpolator: function, None
        If provided, internal calculation is skipped in favor of a
        precalculated interpolating function which takes in sh/hh
        and returns the likelihood.
    return_peak: bool, False
        Return the peak likelihood and index if using passing an array as
        input in addition to the marginalized over the array likelihood.
    return_complex: bool, False
        Return the sh / hh data products before applying phase marginalization.
        This option is intended to aid in reconstucting phase marginalization
        and is unlikely to be useful for other purposes.

    Returns
    -------
    loglr: float
        The marginalized loglikehood ratio
    """
    if distance and not interpolator and not numpy.isscalar(sh):
        raise ValueError("Cannot do vector marginalization "
                         "and distance at the same time")

    if logw is None:
        if isinstance(hh, float):
            logw = 0
        else:
            logw = -numpy.log(len(sh))

    if return_complex:
        pass
    elif phase:
        sh = abs(sh)
    else:
        sh = sh.real

    if interpolator:
        # pre-calculated result for this function
        vloglr = interpolator(sh, hh)

        if skip_vector:
            return vloglr
    else:
        # explicit calculation
        if distance:
            # brute force distance path
            dist_rescale, dist_weights = distance
            sh = sh * dist_rescale
            hh = hh * dist_rescale ** 2.0
            logw = numpy.log(dist_weights)

        if return_complex:
            return sh, -0.5 * hh

        # Apply the phase marginalization
        if phase:
            sh = numpy.log(i0e(sh)) + sh

        # Calculate loglikelihood ratio
        vloglr = sh - 0.5 * hh

    if return_peak:
        maxv = vloglr.argmax()
        maxl = vloglr[maxv]

    # Do brute-force marginalization if loglr is a vector
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    elif not skip_vector:
        vloglr = float(logsumexp(vloglr, b=numpy.exp(logw)))

    if return_peak:
        return vloglr, maxv, maxl
    return vloglr
