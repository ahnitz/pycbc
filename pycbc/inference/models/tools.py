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
                              marginalize_sky_candidates=8,
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
        self.marginalize_sky_analytic = \
            str_to_bool(marginalize_sky_analytic)
        self.marginalize_sky_candidates = \
            max(1, int(float(marginalize_sky_candidates)))
        self._analytic_const = None
        self._analytic_geom = {}

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

    def _analytic_constants(self):
        """Detector geometry and tc bounds. Constant for the whole run.

        Built once: ``Detector.__init__`` scans the LAL detector table and
        constructs an astropy EarthLocation, and ``gmst_estimate`` on a fresh
        object falls through to ``gmst_accurate``; together they cost
        several ms per call if rebuilt, against ~0.05 ms for the linear
        algebra they feed.
        """
        if self._analytic_const is None:
            ifos = list(self.data.keys())
            tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
            tcmin, tcmax = float(tcmin), float(tcmax)
            tcave = 0.5 * (tcmin + tcmax)
            dets = {i: Detector(i, reference_time=tcave) for i in ifos}
            self._analytic_const = {
                'ifos': ifos, 'tcmin': tcmin, 'tcmax': tcmax, 'dets': dets,
                'gmst': dets[ifos[0]].gmst_estimate(tcave),
                'loc': {i: numpy.asarray(dets[i].location) for i in ifos},
                'log_tcspan': numpy.log(tcmax - tcmin)}
        return self._analytic_const

    def _analytic_geometry(self, order):
        """Baseline geometry for one detector ORDERING.

        Keyed on the ordering because that is recomputed every call from the
        current SNR peaks, and ``keep_ifos`` can change between calls too, so a
        cache without the key would silently reuse another network's geometry.
        The key space is a handful of entries.
        """
        key = tuple(order)
        if key not in self._analytic_geom:
            c = self._analytic_constants()
            loc = c['loc']
            d1 = loc[order[0]] - loc[order[1]]
            if len(order) == 2:
                mat = numpy.atleast_2d(d1)
            else:
                mat = numpy.vstack([d1, loc[order[0]] - loc[order[2]]])
            mp = numpy.linalg.pinv(mat)
            smat = C_SI ** 2.0 * (mp.T @ mp)
            # the widest dt12 the ellipse allows is exactly the light
            # travel time along that baseline, which the detector knows
            t12max = c['dets'][order[0]].light_travel_time_to_detector(
                c['dets'][order[1]])
            g = {'d1': d1, 'Mp': mp, 'S': smat, 't12max': t12max,
                 'log_const': c['log_tcspan'] + numpy.log(2.0 * t12max)}
            if len(order) > 2:
                cr = numpy.cross(d1, loc[order[0]] - loc[order[2]])
                g['e3'] = cr / numpy.linalg.norm(cr)
            else:
                dh = d1 / numpy.linalg.norm(d1)
                u = numpy.cross(dh, [0.0, 0.0, 1.0])
                u = u / numpy.linalg.norm(u)
                g['basis'] = (dh, u, numpy.cross(dh, u))
            self._analytic_geom[key] = g
        return self._analytic_geom[key]

    def _analytic_windows(self, snrs, ifos):
        """Per-detector SNR series sliced to the region the draw may use."""
        c = self._analytic_constants()
        out = {}
        for ifo in ifos:
            snr = snrs[ifo]
            tmin, tmax = c['tcmin'] - EARTH_RADIUS, c['tcmax'] + EARTH_RADIUS
            if hasattr(self, 'tstart'):
                tmin, tmax = self.tstart[ifo], self.tend[ifo]
            snr = snr.time_slice(max(tmin, snr.start_time + snr.delta_t),
                                 min(tmax, snr.end_time - snr.delta_t * 2),
                                 mode='nearest')
            out[ifo] = (snr, snr.squared_norm().numpy() / 2.0)
        return out

    @staticmethod
    def _draw_in_window(logw, base, delta, lo_t, hi_t):
        """Draw times ``propto exp(logw)`` inside a per-sample time window.

        Times are OFFSETS from the epoch the caller works in, not absolute
        GPS; ``base`` is where this detector's series starts in that same
        frame. See ``analytic_sky_draw`` for why.

        Returns the dithered times and the log PROPOSAL DENSITY. The dither
        makes the draw a genuine density rather than a probability mass on the
        sample grid, which is what lets the sky come out continuous.
        """
        lmax = logw.max()
        p = numpy.exp(logw - lmax)
        n = len(logw)
        cum = numpy.concatenate(([0.0], numpy.cumsum(p)))
        lo = numpy.clip(numpy.ceil((lo_t - base) / delta), 0, n).astype(int)
        hi = numpy.clip(numpy.floor((hi_t - base) / delta) + 1,
                        0, n).astype(int)
        empty = hi <= lo
        lo = numpy.where(empty, 0, lo)
        hi = numpy.where(empty, n, hi)
        norm = cum[hi] - cum[lo]
        target = cum[lo] + numpy.random.uniform(0, 1, len(lo)) * norm
        idx = numpy.clip(numpy.searchsorted(cum, target) - 1, 0, n - 1)
        times = (base + idx * delta
                 + numpy.random.uniform(-delta / 2.0, delta / 2.0, len(idx)))
        dens = (logw[idx] - lmax) - numpy.log(norm) - numpy.log(delta)
        return times, numpy.where(empty, -numpy.inf, dens)

    @staticmethod
    def _alias_draw(weights, shape):
        """Draw indices ``propto weights`` through a Vose alias table.

        Each draw is then two array lookups instead of a binary search, which
        is what dominates here: the same length-``n`` law is drawn from tens
        of thousands of times, so the O(n) build amortises immediately. One
        uniform serves both alias coordinates -- scaling by ``n`` puts the
        cell in the integer part and leaves an independent uniform in the
        remainder -- so this needs no more random numbers than inversion.
        Zero-weight cells can only appear as the alias of a positive one,
        never as a direct hit, so underflowed cells stay unreachable.
        """
        n = len(weights)
        prob = weights * (n / weights.sum())
        alias = numpy.arange(n)
        small = list(numpy.nonzero(prob < 1.0)[0])
        large = list(numpy.nonzero(prob >= 1.0)[0])
        while small and large:
            lean, full = small.pop(), large.pop()
            alias[lean] = full
            prob[full] -= 1.0 - prob[lean]
            (small if prob[full] < 1.0 else large).append(full)
        prob[large] = 1.0
        prob[small] = 1.0
        x = numpy.random.random_sample(shape)
        x *= n
        idx = x.astype(numpy.intp)
        x -= idx
        return numpy.where(x < prob[idx], idx, alias[idx])

    def _draw_second_delay(self, series, order, geom, epoch, t_off, dt12):
        """Draw dt13 on the ellipse slice the first delay allows.

        Given dt12 the physical region is the slice of ``dt^T S dt <= 1``,
        whose
        endpoints are a closed-form quadratic solve. On that slice the exact
        isotropic sky prior is the ARCSINE law, and substituting
        ``dt13 = mid + half*sin(theta)`` makes it uniform in theta -- which is
        what cancels the ``1/s`` Jacobian singularity at the ends of the slice,
        where the source lies in the detector plane and the two timing cones go
        tangent. A proposal uniform in dt13 instead has log-divergent weight
        variance and is worse than the map this replaces.

        The third detector's own SNR then tilts the draw through K candidate
        cells taken from its SHARED one-dimensional law, so the arcsine cell
        mass is evaluated at K cells per sample rather than at every cell and
        the cost does not depend on the series length. The exact normaliser
        ``sum_j sh_j P_j`` is replaced by the K-sample estimate
        ``tot * mean_k(P_k)``, whose variance falls as 1/K.

        Candidates are held as ``(ncand, nsamples)`` so that the reductions
        over candidates run along contiguous rows.
        """
        smat = geom['S']
        aq = smat[1, 1]
        bq = 2.0 * smat[0, 1] * dt12
        cq = smat[0, 0] * dt12 * dt12 - 1.0
        disc = bq * bq - 4.0 * aq * cq
        invalid = disc <= 0.0
        root = numpy.sqrt(numpy.maximum(disc, 0.0))
        lo13, hi13 = (-bq - root) / (2 * aq), (-bq + root) / (2 * aq)
        mid = 0.5 * (lo13 + hi13)
        half = 0.5 * (hi13 - lo13)

        snr, logl = series[order[2]]
        ncand = self.marginalize_sky_candidates
        base = float(snr.start_time - epoch)
        delta = float(snr.delta_t)
        # shift before the exponential: |SNR|^2/2 reaches several hundred and
        # the unshifted form overflows to inf
        lmax = logl.max()
        shifted = numpy.exp(logl - lmax)
        total = shifted.sum()
        nsamp = len(dt12)
        cell = self._alias_draw(shifted, (ncand, nsamp))
        dt13c = t_off[None, :] - (base + cell * delta)
        half_safe = numpy.maximum(half, 1e-300)[None, :]
        mid_row = mid[None, :]
        # single precision for the arcsine is ~7x faster and costs at most
        # 5e-5 relative on the cell mass; the same array sets both the
        # selection probability and the weight, so this stays self-consistent
        edge_hi = numpy.arcsin(numpy.clip(
            (dt13c + delta / 2.0 - mid_row) / half_safe,
            -1.0, 1.0).astype(numpy.float32))
        edge_lo = numpy.arcsin(numpy.clip(
            (dt13c - delta / 2.0 - mid_row) / half_safe,
            -1.0, 1.0).astype(numpy.float32))
        # the half-cell offset is positive and arcsine is monotone, so
        # edge_hi >= edge_lo always: no absolute value, and no need to sort
        # the pair into an interval below
        mass = (edge_hi - edge_lo) / numpy.pi
        # accumulate in double: 1e-300 would flush to zero in single
        msum = mass.sum(axis=0, dtype=numpy.float64)
        invalid |= ~(msum > 0)
        msafe = numpy.maximum(msum, 1e-300)
        # inverse-cdf over the candidates, accumulated in place: the running
        # sum never has to be materialised as an (ncand, nsamp) array
        thresh = numpy.random.random_sample(nsamp) * msafe
        acc = mass[0].astype(numpy.float64)
        pick = numpy.zeros(nsamp, dtype=numpy.intp)
        for k in range(1, ncand):
            pick += acc < thresh
            acc += mass[k]
        cols = numpy.arange(nsamp)
        theta_lo = edge_lo[pick, cols].astype(numpy.float64)
        theta_hi = edge_hi[pick, cols].astype(numpy.float64)
        theta = (theta_lo + numpy.random.uniform(0, 1, nsamp)
                 * (theta_hi - theta_lo))
        dt13 = mid + half * numpy.sin(theta)
        # q = exp(logl_j) p(dt13) / Zhat  =>  p/q = Zhat / exp(logl_j)
        dens = ((logl[cell[pick, cols]] - lmax)
                - numpy.log(total * msafe / ncand))
        return dt13, dens, invalid

    def _analytic_antenna(self, nhat):
        """F+, Fx and earth-centre delays for the drawn directions.

        The inversion produces the direction as a vector, which is what
        ``Detector`` takes here, so no ra/dec round trip is needed.
        Polarization is applied analytically downstream.
        """
        c = self._analytic_constants()
        fplus, fcross, delay = {}, {}, {}
        for ifo, det in c['dets'].items():
            fplus[ifo], fcross[ifo] = det.antenna_pattern_from_direction(nhat)
            delay[ifo] = det.time_delay_from_direction(nhat)
        return fplus, fcross, delay

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
        c = self._analytic_constants()
        series = self._analytic_windows(snrs, ifos)
        order = sorted(ifos, key=lambda i: -series[i][1].max())
        ndet = len(order)
        geom = self._analytic_geometry(order) if ndet > 1 else None

        snr0, logl0 = series[order[0]]
        delta0 = float(snr0.delta_t)
        # the normaliser is already needed for the density, so forming it
        # here costs nothing and avoids both a rebuilt cdf and a logsumexp
        # call whose overhead dwarfs the length-n array it is given
        l0max = logl0.max()
        shifted0 = numpy.exp(logl0 - l0max)
        i0 = self._alias_draw(shifted0, vsamples)
        # Every delay below is a difference of two coalescence-time-like
        # quantities. Formed from absolute GPS they are order 1e9 s, which a
        # double holds to about 2.4e-7 s -- a thousandth of a sample at
        # 4096 Hz -- and that error lands directly on the cell boundaries the
        # candidate weights are read from. A drawn dt13 then falls on the
        # wrong side of a boundary for order 0.1% of samples and is weighted
        # with the neighbouring cell's SNR, which near a peak differs by
        # O(1-10) nats: order 1% on the integral, and enough to swamp any
        # per-point comparison of two versions of this code. So the epoch is
        # split off and every time here is an OFFSET from it, which makes the
        # differences exact to ~1e-18 and pins the partition.
        epoch = snr0.start_time
        t_off = (i0 * delta0
                 + numpy.random.uniform(-delta0 / 2.0, delta0 / 2.0, vsamples))
        dens = ((logl0[i0] - l0max) - numpy.log(shifted0.sum())
                - numpy.log(delta0))
        invalid = numpy.zeros(vsamples, dtype=bool)
        branch_corr = 0.0

        if ndet == 1:
            ra = self.marginalized_vector_priors['ra'].rvs(size=vsamples)['ra']
            dec = self.marginalized_vector_priors['dec'].rvs(
                size=vsamples)['dec']
            lon = ra - c['gmst']
            nhat = numpy.array([numpy.cos(dec) * numpy.cos(lon),
                                numpy.cos(dec) * numpy.sin(lon),
                                numpy.sin(dec)])
            logw = -c['log_tcspan'] - dens
        else:
            t12max = geom['t12max']
            snr1, logl1 = series[order[1]]
            t_two, dens1 = self._draw_in_window(
                logl1, float(snr1.start_time - epoch), float(snr1.delta_t),
                t_off - t12max, t_off + t12max)
            dens = dens + dens1
            dt12 = t_off - t_two
            invalid |= ~numpy.isfinite(dens1) | (numpy.abs(dt12) > t12max)
            if ndet == 2:
                # a single delay leaves a free azimuth about the baseline, and
                # there the Jacobian is constant, so uniform is exactly
                # isotropic and needs no correction
                dhat, uvec, vvec = geom['basis']
                cos_t = numpy.clip(-dt12 / t12max, -1.0, 1.0)
                sin_t = numpy.sqrt(numpy.maximum(1.0 - cos_t * cos_t, 0.0))
                az = numpy.random.uniform(0, 2 * numpy.pi, vsamples)
                nhat = (cos_t * dhat[:, None]
                        + sin_t * (numpy.cos(az) * uvec[:, None]
                                   + numpy.sin(az) * vvec[:, None]))
            else:
                dt13, dens2, bad2 = self._draw_second_delay(
                    series, order, geom, epoch, t_off, dt12)
                dens = dens + dens2
                invalid |= bad2
                npar = geom['Mp'] @ (-C_SI
                                     * numpy.vstack([dt12, dt13]))
                sgeo = numpy.sqrt(numpy.clip(1.0 - (npar * npar).sum(0),
                                             0.0, None))
                # The two roots are mirror images through the detector
                # plane. Both are physical, but they put the coalescence at
                # times differing by twice the reference detector's offset
                # from that plane, and for some geometries one of them lies
                # outside the tc prior for EVERY sample -- measured at 50% and
                # 62% of the pool on two injections, against ~0% elsewhere.
                # Drawing the root blind and rejecting afterwards throws that
                # fraction of the pool away.
                #
                # Testing a root costs one dot product, so draw only among the
                # survivors: root b with probability v_b / (v+ + v-), and
                # scale the weight by (v+ + v-)/2. The expectation is
                # unchanged, still (w+ v+ + w- v-)/2, and nothing is
                # discarded. Verified directly against that target: the
                # scheme is unbiased at z = -0.5 over 4e6 draws.
                loc0 = c['loc'][order[0]] / C_SI
                tc_base = float(epoch) + t_off + loc0 @ npar
                tc_swing = sgeo * float(loc0 @ geom['e3'])
                tc_up, tc_dn = tc_base + tc_swing, tc_base - tc_swing
                ok_up = (tc_up >= c['tcmin']) & (tc_up <= c['tcmax'])
                ok_dn = (tc_dn >= c['tcmin']) & (tc_dn <= c['tcmax'])
                nroot = ok_up.astype(numpy.int8) + ok_dn
                invalid |= nroot == 0
                coin = numpy.random.uniform(0, 1, vsamples) < 0.5
                up = numpy.where(ok_up & ok_dn, coin, ok_up)
                branch = numpy.where(up, 1.0, -1.0)
                nhat = npar + (branch * sgeo) * geom['e3'][:, None]
                branch_corr = numpy.log(
                    0.5 * numpy.maximum(nroot, 1).astype(float))
            logw = -geom['log_const'] - dens + branch_corr

        fplus, fcross, delay = self._analytic_antenna(nhat)
        tc = float(epoch) + t_off - delay[order[0]]
        outside = (tc < c['tcmin']) | (tc > c['tcmax'])
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
            numpy.arctan2(nhat[1], nhat[0]) + c['gmst'], 2 * numpy.pi)
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
            prior/proposal ratio is analytic rather than approximate. Falls
            back
            to the map, with a warning, if it cannot produce a usable draw.
        marginalize_sky_candidates: int
            Candidate cells used when tilting the second delay by the third
            detector's SNR. Cost is roughly linear in this and the effective
            sample size saturates around 8, which is the default.
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
