""" Common utility functions for calculation of likelihoods
"""

import logging
import warnings
from distutils.util import strtobool

import numpy
import numpy.random
import tqdm

from scipy.special import logsumexp, i0e, softmax
from scipy.interpolate import interp1d
from scipy import ndimage
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


def draw_sample(loglr, size=None, rng=None):
    """ Draw a random index from a 1-d vector with loglr weights

    A generator may be given so drawing does not consume the global random
    state, which would move every sampler's trajectory.
    """
    rng = numpy.random if rng is None else rng
    if size:
        x = rng.uniform(size=size)
    else:
        x = rng.uniform()
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
                              marginalize_orientation=False,
                              marginalize_orientation_ncos=24,
                              marginalize_orientation_npol=12,
                              marginalize_orientation_eps=0.02,
                              marginalize_orientation_temper=2.0,
                              marginalize_sky_analytic=False,
                              marginalize_sky_amplitude=False,
                              marginalize_reconstruct=None,
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
        marginalize_orientation: bool, False
            Draw inclination and polarization conditional on each sky sample
            rather than from the prior. Off by default; the proposal costs
            extra work per draw and only pays off where the orientation is
            well constrained.
        marginalize_orientation_ncos: int, 24
            Number of cos(inclination) cells in the orientation proposal grid.
        marginalize_orientation_npol: int, 12
            Number of polarization cells in the orientation proposal grid.
        marginalize_orientation_eps: float, 0.02
            Floor mixed into the orientation proposal so no cell has zero
            probability.
        marginalize_orientation_temper: float, 2.0
            Flattening applied to the orientation surrogate. The surrogate
            spans ~160 nats across cells within one sky sample, so an
            untempered softmax is nearly a delta.
        marginalize_reconstruct: str, None
            Draw the marginalized parameters while the likelihood is being
            marginalized instead of re-evaluating it afterwards. Names the
            level to demarginalize to: 'vector', 'distance' or 'phase', each
            including the ones above it. Default None draws nothing inline.

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

        # the effective sample size of the most recent vector
        # marginalization, filled in by marginalize_loglr; nan until then
        self.vector_ess = numpy.nan
        # ... and a running summary over the whole run, because the
        # per-call value was being computed and thrown away. A starved
        # extrinsic marginalization injects variance straight into the
        # importance weights, so this is the quantity that says whether a
        # run's efficiency is capped by the sky integral rather than by the
        # intrinsic proposal. Kept as accumulators, not a list, so a 60k
        # call run costs nothing to instrument.
        self.vess_n = 0
        self.vess_sum = 0.0
        self.vess_min = numpy.inf
        self.vess_max = 0.0
        # each level is conditional on the point the one before it drew, so
        # asking for a level asks for the ones above it too
        levels = ['vector', 'distance', 'phase']
        self.reconstruct_inline = levels[:levels.index(
            marginalize_reconstruct) + 1] if marginalize_reconstruct else []
        # drawing happens every likelihood call, so it gets its own generator
        # rather than moving the sampler's random state
        self.marginalize_rng = (numpy.random.default_rng()
                                if self.reconstruct_inline else None)

        # Handle any requested parameter vector / brute force marginalizations
        self.marginalize_vector_params = {}
        self.marginalized_vector_priors = {}
        self.vsamples = int(marginalize_vector_samples)

        # a generator kept for the one draw that is hot enough to matter,
        # seeded from the global stream so a run seeded with
        # numpy.random.seed still repeats. Choosing a subset without
        # replacement through numpy.random.choice permutes the whole
        # precomputed set to take a slice of it; the generator does it by
        # Floyd's algorithm, the same draw without the permutation.
        self._choice_rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 63))

        self.marginalize_sky_initial_samples = \
            int(float(marginalize_sky_initial_samples))
        # PCG64 rather than the legacy global state: it is faster per
        # draw, and seeding it from that state leaves numpy.random.seed
        # pinning a run as before
        self._sky_rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 63))
        self.marginalize_sky_analytic = \
            str_to_bool(marginalize_sky_analytic)
        # weight the free azimuth by the amplitudes the detectors saw,
        # rather than sweeping it flat
        self.marginalize_sky_amplitude = \
            str_to_bool(marginalize_sky_amplitude)
        self._sky_detectors = None

        self.marginalize_orientation = str_to_bool(marginalize_orientation)
        self.marginalize_orientation_ncos = int(marginalize_orientation_ncos)
        self.marginalize_orientation_npol = int(marginalize_orientation_npol)
        self.marginalize_orientation_eps = float(marginalize_orientation_eps)
        self.marginalize_orientation_temper = \
            float(marginalize_orientation_temper)

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

        # The per-sample vector is what BOTH the effective sample size and
        # the inline draw are built from, so ask for it once. The peak and
        # reconstruction paths do not marginalize a vector here.
        want_vector = not (return_peak or return_complex or skip_vector)
        out = marginalize_likelihood(sh_total, hh_total,
                                     logw=self.marginalize_vector_weights,
                                     phase=self.marginalize_phase,
                                     interpolator=interpolator,
                                     distance=distance,
                                     skip_vector=skip_vector,
                                     return_complex=return_complex,
                                     return_peak=return_peak,
                                     return_vector=want_vector)
        if not want_vector:
            return out
        out, vector = out
        self._record_vector_ess(vector)
        if self.reconstruct_inline:
            # NO SPLIT. Upstream halves the samples, giving the likelihood one
            # half and the draw the other, so that a draw cannot inherit the
            # luck of the points the sampler kept. That is the right thing to
            # measure in a test, but here it would make ASKING FOR THE OUTPUT
            # CHANGE THE ANSWER: the returned loglr would rest on half the
            # vector samples, and the extrinsic marginalization is already
            # starved -- measured at an effective 11 to 500 of 15000, so
            # halving it is not free. The draw therefore uses every sample and
            # `out` is the full-vector marginalization, identical to a run
            # that reconstructs nothing.
            self.draw_inline(sh_total, hh_total, vector, None)
        return out

    def _record_vector_ess(self, vector):
        """ Effective sample size of the extrinsic marginalization, kept as
        running accumulators over the run.

        (sum w)^2 / sum w^2 with w the combined weights. Low against the
        number drawn means the marginal rests on a handful of them, which is
        the quantity that says whether a run is capped by the sky integral
        rather than by its intrinsic proposal.
        """
        lw = vector + self.marginalize_vector_weights
        if not numpy.ndim(lw) or not numpy.size(lw):
            return
        lwmax = lw.max()
        if not numpy.isfinite(lwmax):
            self.vector_ess = numpy.nan
            return
        w = numpy.exp(lw - lwmax)
        v = float(w.sum() ** 2.0 / numpy.vdot(w, w))
        self.vector_ess = v
        if numpy.isfinite(v):
            self.vess_n += 1
            self.vess_sum += v
            self.vess_min = min(self.vess_min, v)
            self.vess_max = max(self.vess_max, v)

    def premarg_draw(self):
        """ Choose random samples from prechosen set"""

        # Update the current proposed times and the marginalization values
        logw = self.premarg['logw_partial']
        if self.vsamples == len(logw):
            choice = slice(None, None)
        else:
            choice = self._choice_rng.choice(len(logw), size=self.vsamples,
                                             replace=False)

        for k in self.snr_params:
            self.marginalize_vector_params[k] = self.premarg[k][choice]

        self._current_params.update(self.marginalize_vector_params)
        self.sample_idx = self.premarg['sample_idx'][choice]

        # The precalculated points were themselves drawn from the sky grid,
        # so the stored antenna factors describe this draw as well and go
        # back with the indices that select it. They have to be restored
        # here and not merely left alone: snr_draw discards them whenever a
        # scalar coalescence time comes through, as it does during
        # reconstruction, and this is the only place that puts them back.
        self.precalc_antenna_factors = self.premarg['antenna_factors']

        # Update the importance weights for each vector sample. These must
        # not be renormalized to sum to one: that would turn the integral
        # over the marginalized parameter into an average over it, and the
        # size of the prior would cancel out of the answer.
        self.marginalize_vector_weights = \
            self.marginalize_vector_weights + logw[choice]

        if self.marginalize_orientation:
            corr = self._draw_orientation()
            if corr is not None:
                self.marginalize_vector_weights = \
                    self.marginalize_vector_weights + corr
                self._current_params.update(self.marginalize_vector_params)
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
            # marginalizations. The stored factors and the indices that
            # select from them describe a draw that no longer applies, and
            # they are invalidated together so that neither can be used
            # without the other.
            self.precalc_antenna_factors = None
            self.sample_idx = None
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

        # Update the current proposed times and the marginalization values.
        # The times were drawn from the incoherent likelihood rather than
        # from the prior, so each carries the ratio of the two. logweight
        # is a probability per time sample, so the sample spacing turns it
        # into a density that the prior can be divided by. Assumes a
        # uniform prior on tc, as the draws themselves do.
        logw = (- logweight[tci]
                + numpy.log(snr.delta_t / (tcmax - tcmin)))
        # times the search had to leave out carry no prior weight
        logw[(tc < tcmin) | (tc > tcmax)] = -numpy.inf
        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['logw_partial'] = logw

        if self._current_params is not None:
            # Update the importance weights for each vector sample
            self._current_params.update(self.marginalize_vector_params)
            self.marginalize_vector_weights += logw

        return self.marginalize_vector_params

    @staticmethod
    def _build_sky_lookup(dmap, most_cells=5_000_000):
        """Flatten dmap into arrays the per-evaluation draw can index.

        dmap maps a tuple of rounded inter-detector delays to the list of
        sky-sample indices with those delays. Per evaluation the draw has
        to, for each drawn delay tuple, find that list and pick one of its
        members. Done as a python loop over the samples that is the bulk
        of the sky marginalization's cost.

        This encodes each delay tuple as one integer and lays the lists
        end to end, so the loop becomes array indexing: ``offset[code]``
        and ``count[code]`` give where a bin's members start and how many
        there are in ``flat``, in the order dmap stored them, so a draw
        reproduces ``dmap[t][int(r * len(dmap[t]))]`` exactly. Returns
        None, leaving the loop to run, when there are no delays to key on
        or the integer code would range over more cells than it is worth
        allocating for.
        """
        keys = [k for k in dmap if len(k) > 0]
        if not keys:
            return None
        key_arr = numpy.array(keys, dtype=numpy.int64)
        lo = key_arr.min(axis=0)
        dims = (key_arr.max(axis=0) - lo + 1).astype(numpy.int64)
        strides = numpy.concatenate(
            [[1], numpy.cumprod(dims[:-1])]).astype(numpy.int64)
        ncell = int(dims.prod())
        if ncell > most_cells or ncell <= 0:
            return None

        offset = numpy.full(ncell, -1, dtype=numpy.int64)
        count = numpy.zeros(ncell, dtype=numpy.int64)
        flat = []
        pos = 0
        for t, members in dmap.items():
            if len(t) == 0:
                continue
            code = int(((numpy.array(t, dtype=numpy.int64) - lo)
                        * strides).sum())
            offset[code] = pos
            count[code] = len(members)
            flat.extend(members)
            pos += len(members)
        return {'lo': lo, 'strides': strides, 'dims': dims,
                'offset': offset, 'count': count,
                'flat': numpy.array(flat, dtype=numpy.int64)}

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

    # the azimuth is drawn from this many equal cells when the amplitude
    # tilt is on, and the statistic is read at each cell's centre. That
    # only works while a cell is narrower than the statistic: at sixteen
    # cells the centres are 0.39 rad apart and the statistic is 0.04 rad
    # wide, so the centres can miss its peak entirely and the cell
    # probabilities come out of its tails. Measured against the posterior
    # on eight rings, sixteen cells is *worse* than drawing the azimuth
    # flat (median 0.2x) where a hundred and twenty-eight is 34.8x and
    # never below 6.4x. Sub-sampling inside coarse cells fixes the same
    # fault but buys less per evaluation: sixteen cells read at eight
    # points each costs what a hundred and twenty-eight centres cost and
    # returns 7.7x against 34.8x.
    AZ_CELLS = 128
    # the share of the azimuth draw left flat. The statistic is sharp --
    # a cell can hold effectively none of the probability -- and drawing
    # such a cell anyway would hand it a weight bounded only by how small
    # that probability was. Mixing a flat part in bounds the weight at
    # 1/AZ_FLOOR, so no single sample can carry the estimate, and costs
    # that fraction of the concentration.
    AZ_FLOOR = 0.01

    def _ring_response(self, dets, order, frame, cos_t):
        """``F+`` and ``Fx`` at zero polarization, around these rings.

        Only a handful of rings are asked for per call, so they are
        worked out directly rather than held in a table over every ring a
        sample might land on.

        ``cos_t`` is those rings; the result is (rings, cells) for each
        of the two detectors.
        """
        dhat, uvec, vvec = frame
        sin_t = numpy.sqrt(numpy.maximum(1.0 - cos_t * cos_t, 0.0))
        phi = (numpy.arange(self.AZ_CELLS) + 0.5) * (2.0 * numpy.pi
                                                     / self.AZ_CELLS)
        nhat = (cos_t[:, None] * dhat[:, None, None]
                + sin_t[:, None] * (numpy.cos(phi) * uvec[:, None, None]
                                    + numpy.sin(phi) * vvec[:, None, None]))
        shape = (len(cos_t), self.AZ_CELLS)
        out = []
        for ifo in order[:2]:
            fplus, fcross = dets[ifo].antenna_pattern_from_direction(
                nhat.reshape(3, -1))
            out.append((fplus.reshape(shape), fcross.reshape(shape)))
        return out

    def _orientation(self):
        """The inclination and polarization to predict amplitudes with.

        Returns None if there is no point being evaluated yet, or if
        either angle is itself being marginalized over; the caller takes
        that as "do not tilt". The amplitude a source of
        given distance produces depends on both, so with either unknown
        there is nothing to compare the observation against. A stand-in
        does not rescue it: averaging a polarization drawn across its
        whole range gives an angle that points nowhere in particular, and
        the draw then commits to it -- measured, that turns a 2.2
        effective sample draw into a 1.0 one.
        """
        params = self._current_params
        if params is None:
            # the precalculated draw runs before there are any: it is
            # laying out sky points, not evaluating a point
            return None
        incl, psi = params['inclination'], params['polarization']
        if not (numpy.isscalar(incl) and numpy.isscalar(psi)):
            return None
        return numpy.cos(float(incl)), float(psi)

    def _amplitude_usable(self):
        """Whether the amplitudes can be predicted here at all.

        Two things are needed and neither is always present. The
        per-detector template norms say how loud a source of this
        distance should look, and only a model that keeps them can
        support the tilt. The inclination and polarization say how the
        response divides between the two polarizations, and either may
        itself be marginalized over, or not be set yet.

        Falling back to a flat azimuth is the honest answer in all of
        those; it is said once rather than silently, so that asking for
        the tilt and not getting it is visible.
        """
        if not hasattr(self, 'hh'):
            if not getattr(self, '_tilt_unavailable', False):
                self._tilt_unavailable = True
                logging.warning(
                    "marginalize_sky_amplitude needs the per-detector "
                    "template norms, which this model does not keep; "
                    "drawing the azimuth flat instead")
            return False
        return self._orientation() is not None

    def _amplitude_tilt(self, wsq, zsq, norm):
        """Log weight per cell from the amplitudes each detector saw.

        The model is being evaluated at a known distance and orientation,
        so the amplitude it predicts is a number and not a scale to be
        fitted: ``|h_i|^2 = |w_i|^2 <h0_i|h0_i>``. Comparing the observed
        amplitude against that, with each detector's own phase
        marginalized out, is the whole of what the amplitudes say about
        where on the ring the source is,

            sum_i  |w_i| |sh_i| - |w_i|^2 <h0_i|h0_i> / 2

        The observed side is a magnitude, so the phase is maximized over
        rather than marginalized; marginalizing it would add
        -log(2 pi x)/2, which a proposal does not keep.

        The scale matters as much as the ratio between the detectors:
        keeping only the ratio makes the weight some forty times too
        broad, since the model knows how loud the signal should be.
        """
        def term(k):
            # in place: these are (rings, cells) and the temporaries cost
            # more than the arithmetic in them
            out = numpy.sqrt(wsq[k])
            out *= numpy.sqrt(zsq[k] * norm[k])[:, None]
            out -= (0.5 * norm[k]) * wsq[k]
            return out

        return term(0) + term(1)

    def _predicted_wsq(self, table):
        """``|w|^2`` per detector, for the orientation being evaluated.

        ``w = F+(psi) (1 + cos i^2)/2 + i Fx(psi) cos i`` is the factor
        the intrinsic waveform is multiplied by, so this is the square of
        the amplitude the model predicts, up to the template norm. Face
        on it is ``F+^2 + Fx^2`` and the polarization drops out of it
        entirely; edge on the cross term is gone and only ``F+(psi)``
        survives, which is why the polarization matters most where the
        inclination matters least.
        """
        cosi, psi = self._orientation()
        plus, cross = 0.5 * (1.0 + cosi * cosi), cosi
        c2p, s2p = numpy.cos(2.0 * psi), numpy.sin(2.0 * psi)
        return [(((fp * c2p + fc * s2p) * plus) ** 2.0
                 + ((-fp * s2p + fc * c2p) * cross) ** 2.0)
                for fp, fc in table]

    # the azimuth distribution is worked out in full for this many rings
    # and every sample draws from a blend of them. The rings the samples
    # land on are set by the first delay, which the SNR peak confines to
    # a narrow span, so a handful of anchors across that span carry the
    # variation and the cost stops scaling with the number of samples.
    # Measured over eight skies: three anchors lose a case outright
    # (0.42x), five recover most of it, nine match drawing every sample
    # its own distribution -- 11.3x median against 11.1x, 5.3x worst
    # against 5.4x -- and more buys nothing. The anchors themselves are
    # nearly free; what they replace is an array of samples by cells.
    AZ_ANCHORS = 9

    def _tilted_azimuth(self, dets, order, frame, cos_t, zsq):
        """Draw the free azimuth favouring the amplitudes actually seen.

        The delays fix which ring the source is on but say nothing about
        where along it; the amplitudes do, because the response varies
        around the ring. Written out, the weight of an azimuth is

            -|z - a(phi)|^2 / 2  + constant

        the squared distance in the plane of amplitudes between what the
        detectors measured and what a source there would produce. So the
        draw wants the part of the ring passing closest to the observed
        pair.

        Doing that per sample means an array of samples by cells, and the
        cells outnumber anything useful in it: the distribution is sharp
        and nearly the same from one sample to the next, because the
        first delay confines them to neighbouring rings. So the
        distribution is built in full on a few anchor rings spanning the
        samples' range, and each sample draws from a blend of those
        weighted by where it falls between them. The blend is a genuine
        mixture -- an anchor is chosen, then a cell from that anchor's
        cumulative sum -- so the density is the mixture's own and the
        weight stays exact however crude the anchoring is.
        """
        width = 2.0 * numpy.pi / self.AZ_CELLS
        norm = [float(self.hh[i]) for i in order[:2]]
        edges = numpy.linspace(cos_t.min(), cos_t.max(), self.AZ_ANCHORS)
        table = self._ring_response(dets, order, frame, edges)
        # one representative pair of amplitudes: they come from times
        # near each detector's own peak, so they vary far less than the
        # rings do
        # a mean rather than a median: it says the same thing about a
        # representative amplitude and does not sort the samples
        mid = tuple(numpy.full(self.AZ_ANCHORS, z.mean()) for z in zsq)
        logw = self._amplitude_tilt(self._predicted_wsq(table),
                                    mid, norm)
        prob = ((1.0 - self.AZ_FLOOR) * softmax(logw, axis=1)
                + self.AZ_FLOOR / self.AZ_CELLS)
        cum = numpy.cumsum(prob, axis=1)

        # where each sample falls between the anchors
        span = max(edges[-1] - edges[0], 1e-30)
        pos = numpy.clip((cos_t - edges[0]) / span, 0.0, 1.0) * (
            self.AZ_ANCHORS - 1)
        low = numpy.clip(pos.astype(int), 0, self.AZ_ANCHORS - 2)
        frac = pos - low
        rng = self._sky_rng
        pick = numpy.where(rng.random(len(cos_t)) < frac, low + 1, low)
        # each anchor's cumulative sum runs 0 to 1, so offsetting the
        # jth by j lays them end to end into one increasing array. A
        # sample looking for `draw` in anchor `pick` looks for
        # `pick + draw` in that, which is one search for all of them
        # rather than one per anchor -- and never an array of samples by
        # cells, which is what the anchors exist to avoid.
        flat = (cum + numpy.arange(self.AZ_ANCHORS)[:, None]).ravel()
        idx = numpy.searchsorted(flat, pick + rng.random(len(cos_t)))
        cell = numpy.clip(idx - pick * self.AZ_CELLS,
                          0, self.AZ_CELLS - 1)
        # the density is the mixture's, not the chosen anchor's
        # every cell keeps at least AZ_FLOOR / AZ_CELLS of the
        # probability, so the density here is never near zero and the
        # logarithm needs no guard
        dens = (1.0 - frac) * prob[low, cell] + frac * prob[low + 1, cell]
        az = (cell + rng.random(len(cell))) * width
        return az, numpy.log(dens * self.AZ_CELLS)

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
                # dt12 can exceed the crossing on a sample already marked
                # invalid, so the cosine is clipped; the sine is then real
                # without a guard of its own, since c*c rounds to at most
                # one for every double in [-1, 1]
                cos_t = numpy.clip(-dt12 / t12max, -1.0, 1.0)
                sin_t = numpy.sqrt(1.0 - cos_t * cos_t)
                if (self.marginalize_sky_amplitude
                        and self._amplitude_usable()):
                    az, dens_az = self._tilted_azimuth(
                        dets, order, (dhat, uvec, vvec), cos_t,
                        (2.0 * logl0[i0], 2.0 * logl1[i1]))
                    dens += dens_az
                else:
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
                # where each root puts the coalescence, still as an offset
                # from the epoch: the bounds move to offsets rather than
                # the times moving to GPS, which would cost the 2.4e-7 s
                loc0 = loc[order[0]] / C_SI
                centre = t_off + loc0 @ npar
                swing = perp * (loc0 @ e3)
                pair = numpy.stack([centre + swing, centre - swing])
                ok_up, ok_dn = ((pair >= tcmin - epoch_t)
                                & (pair <= tcmax - epoch_t))
                both = ok_up & ok_dn
                invalid |= ~(ok_up | ok_dn)
                # with both roots in the prior take either; with one, take
                # it, and the half below no longer cancels the pair
                up = numpy.where(both, self._sky_rng.random(vsamples) < 0.5,
                                 ok_up)
                nhat = npar + numpy.where(up, perp, -perp) * e3[:, None]
                branch_corr = numpy.where(both, 0.0, -numpy.log(2.0))
            logw = (-log_tcspan - numpy.log(2.0 * t12max)
                    - dens + branch_corr)

        fplus, fcross, delay = {}, {}, {}
        for ifo, det in dets.items():
            fplus[ifo], fcross[ifo] = det.antenna_pattern_from_direction(nhat)
            delay[ifo] = det.time_delay_from_direction(nhat)
        # tc has to be absolute; the test against the prior does not, and
        # in offsets it agrees exactly with the one the roots were chosen
        # by rather than to the 2.4e-7 s a double resolves near 1.19e9
        tc_off = t_off - delay[order[0]]
        outside = (tc_off < tcmin - epoch_t) | (tc_off > tcmax - epoch_t)
        tc = epoch_t + tc_off
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

        ndraw = vsamples
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
                fp[ifo], fc[ifo], dtc[ifo] = \
                    d[ifo].antenna_pattern_and_delay(ra, dec, 0.0, tcave)

            dmap = {}
            for i, t in enumerate(tqdm.tqdm(zip(*dts))):
                if t not in dmap:
                    dmap[t] = []
                dmap[t].append(i)

            if len(ifos) == 1:
                dmap[()] = numpy.arange(0, size, 1).astype(int)

            # Sky prior by bin
            bin_prior = {t: len(dmap[t]) / size for t in dmap}

            # a flat form of dmap that the per-evaluation draw can index
            # without a python loop over the samples; see _build_sky_lookup
            lookup = self._build_sky_lookup(dmap)

            return (dmap, tcmin, tcmax, fp, fc, ra, dec, dtc, bin_prior,
                    lookup)

        if not hasattr(self, 'tinfo'):
            self.tinfo = {}

        if ikey not in self.tinfo:
            logging.info('pregenerating sky pointings')
            self.tinfo[ikey] = make_init()

        (dmap, tcmin, tcmax, fp, fc, ra, dec, dtc, bin_prior,
         lookup) = self.tinfo[ikey]

        # draw times from each snr time series
        # Is it worth doing this if some detector has low SNR?
        sref = None
        iref = None
        idx = []
        dx = []
        mcweight = None
        _saved = {}
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
            _saved[ifo] = (snr, w.copy())
            i = draw_sample(w, size=ndraw)
            # the times were drawn from this detector's series in
            # proportion to its likelihood, so normalizing over that
            # series is what turns it into the probability they were
            # drawn with. Normalizing over the drawn times instead would
            # make it a probability of one in however many were drawn,
            # and the answer would grow with the number of them.
            w -= logsumexp(w)

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

        # ADAPTIVE B2: only the atoms that actually miss pay for oversampling.
        # A first pass at `vsamples` is the current cost; if too few draws land
        # in a populated cell we redraw once with a count estimated from the
        # observed acceptance fraction. Healthy atoms therefore stay at 1x.
        def _draw_delays(n):
            _sref = _iref = None
            _idx, _dx, _mcw = [], [], None
            for _ifo in ifos:
                _snr, _w = _saved[_ifo]
                _i = draw_sample(_w, size=n)
                _wn = _w - logsumexp(_w)
                if _sref is not None:
                    _mcw = _mcw + _wn[_i]
                    _delt = float(_snr.start_time - _sref.start_time)
                    _i = _i + round(_delt / _sref.delta_t)
                    _dx.append(_iref - _i)
                else:
                    _sref, _iref, _mcw = _snr, _i, _wn[_i]
                _idx.append(_i)
            return _sref, _iref, _dx, _mcw

        # check if delay is in dict, if not, throw out
        rand = numpy.random.uniform(0, 1, size=ndraw)
        if lookup is not None and len(dx) > 0:
            # the loop below, without the python: encode every drawn delay
            # tuple as one integer, keep those whose bin is populated, and
            # index the flattened bins. Bit-identical to the loop, which
            # is checked in test/test_marg_sky_vectorized.py.
            lo, strides, dims = lookup['lo'], lookup['strides'], lookup['dims']
            shifted = numpy.stack(dx, axis=1) - lo
            inside = numpy.all((shifted >= 0) & (shifted < dims), axis=1)
            codes = numpy.where(
                inside, (shifted * strides).sum(axis=1), 0).astype(numpy.int64)
            present = inside & (lookup['offset'][codes] >= 0)
            ti = numpy.nonzero(present)[0]
            codes_k = codes[ti]
            counts_k = lookup['count'][codes_k]
            ix = lookup['flat'][lookup['offset'][codes_k]
                                + (rand[ti] * counts_k).astype(numpy.int64)]
            # bin_prior[t] is len(dmap[t]) / marginalize_sky_initial_samples;
            # count is that length, so divide by the same denominator, not
            # by the vsamples that the local ``size`` refers to here
            wi = counts_k / self.marginalize_sky_initial_samples
            ti = list(ti)
            ix = list(ix)
            wi = list(wi)

            # Bias WHICH POINT we take from each delay bin, using
            # amplitude/phase consistency.
            #
            # The delay match fixes the bin; within it the members differ in
            # sub-sample delay and antenna response, and the relative
            # amplitudes/phases across detectors say which is right - the only
            # sky information the delay draw never uses. For a member the model
            # in detector i is (a+ F+_i + ax Fx_i) h+, so with
            #   y_k = sum_i F_ki sh_i ,  G_kl = sum_i F_ki F_li hh_i
            # the consistency is coh = y^dag G^-1 y / 2, already maximised over
            # distance, phase, inclination and polarization.
            #
            # Selecting member j with probability p_j instead of 1/count means
            # the prior/proposal contribution becomes 1/(N_initial p_j) rather
            # than count/N_initial, so the estimator is unchanged.
        else:
            ti = []
            ix = []
            wi = []
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
        resize_factor = len(ti) / ndraw

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

        # Update the current proposed times and the marginalization values.
        # The points were drawn from the incoherent likelihood rather than
        # from the prior, so each carries the ratio of the two. One tuple
        # of time samples gives one (tc, ra, dec): the reference detector
        # fixes tc and the delays between the detectors fix the sky bin.
        # So mcweight is a probability per time sample per sky bin, and
        # the sample spacing and the bin's share of the sky prior turn it
        # into a density that the tc prior can be divided by. Assumes a
        # uniform prior on tc, as the draws themselves do. Filling back up
        # to vsamples above makes what follows an average over the draws
        # that landed in a populated bin rather than a sum over all of
        # them, so resize_factor puts back the ones left out.
        # The delta_t/window factor converts the drawn time-sample probability
        # into a density the uniform tc prior can be divided by. Without it the
        # evidence acquires a sample-rate dependence of -1.23 nats per 4x in
        # rate (measured); with it the residual is +0.16, so the derivation is
        # close but not exact -- mcweight carries delta_t^ndet while this
        # carries delta_t^1, and the full Jacobian from the ndet drawn times to
        # (tc, ra, dec) is what would make it exactly rate-invariant.
        logw_sky = (-mcweight[ti] + numpy.log(wi)
                    + numpy.log(sref.delta_t / (tcmax - tcmin))
                    + numpy.log(resize_factor))
        # times the search had to leave out carry no prior weight
        logw_sky[(tc < tcmin) | (tc > tcmax)] = -numpy.inf

        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['ra'] = ra
        self.marginalize_vector_params['dec'] = dec

        if self.marginalize_orientation:
            corr = self._draw_orientation()
            if corr is not None:
                logw_sky = logw_sky + corr

        self.marginalize_vector_params['logw_partial'] = logw_sky

        if self._current_params is not None:
            # Update the importance weights for each vector sample
            self._current_params.update(self.marginalize_vector_params)
            self.marginalize_vector_weights += logw_sky

        return self.marginalize_vector_params

    def _orientation_grid(self):
        """Cell centres and coefficients for the (iota,psi) proposal grid.

        These depend only on the grid, so they are built once and reused.
        """
        nc = self.marginalize_orientation_ncos
        npsi = self.marginalize_orientation_npol
        key = (nc, npsi)
        if getattr(self, '_orient_grid_key', None) == key:
            return self._orient_grid
        ce = numpy.linspace(-1.0, 1.0, nc + 1)
        pe = numpy.linspace(0.0, 2.0 * numpy.pi, npsi + 1)
        cm = 0.5 * (ce[:-1] + ce[1:])
        pm = 0.5 * (pe[:-1] + pe[1:])
        CM, PM = numpy.meshgrid(cm, pm, indexing='ij')
        CM = CM.ravel()
        PM = PM.ravel()
        ip = 0.5 * (1.0 + CM ** 2)
        c2 = numpy.cos(2.0 * PM)
        s2 = numpy.sin(2.0 * PM)
        k1 = ip ** 2 * c2 ** 2 + CM ** 2 * s2 ** 2
        k2 = ip ** 2 * s2 ** 2 + CM ** 2 * c2 ** 2
        k3 = 2.0 * c2 * s2 * (ip ** 2 - CM ** 2)
        k4 = -2.0 * ip * CM
        self._orient_grid = (nc, npsi, ce, pe,
                             numpy.stack([k1, k2, k3, k4]),
                             numpy.stack([k1, k2, k3]))
        self._orient_grid_key = key
        return self._orient_grid

    def _draw_orientation(self):
        """Draw (inclination, polarization) conditional on each sky sample.

        At a fixed sky location the whole orientation dependence of the
        likelihood is carried by five numbers,
            Y_p = sum_i sh_i fp_i,  Y_c = sum_i sh_i fc_i     (complex)
            G_pp = sum_i hh_i fp_i^2, G_pc, G_cc              (real)
        because with f_i = (fp_i + i fc_i) exp(-2i psi),
            sh_total = ip (c2 Y_p + s2 Y_c) + i cos(iota) (-s2 Y_p + c2 Y_c)
        and Im(A* B) turns out to be Im(Y_p* Y_c), independent of psi. So
        |sh|^2 and hh_total are bilinear: the sky sample supplies seven real
        numbers and the grid cell supplies four coefficients, which makes the
        whole grid two small matrix products and no SNR gathers beyond the
        one already needed per detector.

        The grid only defines a piecewise-constant proposal. The draw is
        uniform WITHIN the chosen cell and every cell keeps a defensive floor,
        so the support is the full continuous (iota, psi) space. The returned
        log(prior/proposal) makes the estimator identical to drawing from the
        prior.
        """
        vp = self.marginalize_vector_params
        if ('inclination' not in vp or 'polarization' not in vp
                or 'tc' not in vp):
            return None
        if not hasattr(self, 'sh') or not hasattr(self, 'hh'):
            return None
        tc = numpy.asarray(vp['tc'], dtype=float)
        if tc.ndim == 0 or tc.size < 2:
            return None
        try:
            ifos = list(self.sh.keys())
            n = tc.size
            Yp = numpy.zeros(n, dtype=complex)
            Yc = numpy.zeros(n, dtype=complex)
            Gpp = numpy.zeros(n)
            Gpc = numpy.zeros(n)
            Gcc = numpy.zeros(n)
            for ifo in ifos:
                fp, fc, dt = self.get_precalc_antenna_factors(ifo)
                sh = numpy.asarray(self.sh[ifo].at_time(
                    tc + dt, interpolate='quadratic', extrapolate=0.0j))
                hh = float(self.hh[ifo])
                Yp += sh * fp
                Yc += sh * fc
                Gpp += hh * fp * fp
                Gpc += hh * fp * fc
                Gcc += hh * fc * fc
        except Exception as e:
            logging.debug('orientation draw unavailable: %s', e)
            return None

        nc, npsi, ce, pe, Wsh, Whh = self._orientation_grid()
        ncell = nc * npsi
        eps = self.marginalize_orientation_eps
        YY = Yp.conj() * Yc
        Xsh = numpy.stack([Yp.real ** 2 + Yp.imag ** 2,
                           Yc.real ** 2 + Yc.imag ** 2,
                           YY.real, YY.imag], axis=1)
        Xhh = numpy.stack([Gpp, Gcc, Gpc], axis=1)
        sh2 = (Xsh @ Wsh).T.astype(numpy.float32)
        hh2 = (Xhh @ Whh).T.astype(numpy.float32)
        numpy.maximum(sh2, 1e-30, out=sh2)
        numpy.maximum(hh2, 1e-30, out=hh2)
        # Laplace value of int u^-4 exp(|sh| u - hh u^2 / 2) du with u = 1/d:
        # the prior VOLUME containment under a uniform-in-volume distance
        # prior, which is what the proposal should track. It is only the
        # proposal, so the approximation costs efficiency and not accuracy.
        # Laplace value of int u^-4 exp(|sh| u - hh u^2/2) du, u = 1/d:
        # |sh|^2/(2hh) + 4 log(hh/|sh|) = sh2/(2 hh2) + 4 log hh2 - 2 log sh2.
        ll = (sh2 / (2.0 * hh2)
              + 4.0 * numpy.log(hh2) - 2.0 * numpy.log(sh2))
        # T=2 measured on this branch: vector ESS over the prior draw is
        # 16.3x (00021), 21.7x (00008), 3.4x (00006) at T=2 against
        # 12.8x / 14.8x / 3.6x untempered. The surrogate spans ~160 nats
        # across cells within one sky sample, so softmax(ll) is nearly a
        # delta and needs flattening. (An earlier T=4 was measured before
        # the logw_sky fix and does not hold here.)
        temper = self.marginalize_orientation_temper
        if temper != 1.0:
            ll = ll / temper
        ll[~numpy.isfinite(ll)] = -numpy.inf

        mx = ll.max(axis=0)
        E = numpy.exp(ll - mx[None, :], out=ll)
        Z = E.sum(axis=0)
        bad = ~numpy.isfinite(Z) | (Z <= 0)
        if bad.any():
            E[:, bad] = 1.0
            Z = numpy.where(bad, float(ncell), Z)
        cum = numpy.cumsum(E, axis=0)
        u = numpy.random.random(n).astype(numpy.float32) * Z
        pick = numpy.argmax(cum >= u[None, :], axis=0)
        floor = numpy.random.random(n) < eps
        if floor.any():
            pick = numpy.where(floor, numpy.random.randint(0, ncell, n), pick)
        pk = E[pick, numpy.arange(n)] / Z
        p_pick = (1.0 - eps) * pk + eps / ncell

        ic = pick // npsi
        ip_ = pick % npsi
        ci = ce[ic] + numpy.random.uniform(0, 1, n) * (ce[ic + 1] - ce[ic])
        ps = pe[ip_] + numpy.random.uniform(0, 1, n) * (pe[ip_ + 1] - pe[ip_])
        cell_area = (2.0 / nc) * (2.0 * numpy.pi / npsi)

        vp['inclination'] = numpy.arccos(numpy.clip(ci, -1.0, 1.0))
        vp['polarization'] = ps
        logprior = -numpy.log(2.0) - numpy.log(2.0 * numpy.pi)
        logq = numpy.log(p_pick) - numpy.log(cell_area)
        corr = logprior - logq
        corr[~numpy.isfinite(corr)] = -numpy.inf
        return corr

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
                        peak_lock_search_samples=None,
                        peak_lock_search_decimate=8,
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
        peak_lock_search_samples: int
            How far to search for the peak before each likelihood, in
            samples, centered on the region locked here. The region is then
            moved to wherever the peak turned out to be. Off if not given.
            Counted in samples, as peak_lock_region is, so that it
            follows the sample rate the model was given rather than
            having to be restated for each one.
        peak_lock_search_decimate: int
            Stride of that search, in samples. It only has to place the
            peak inside the region, so eight is enough for any region at
            least eight samples wide, and a shorter stride costs more
            without improving the answer.
        """
        self.peak_lock_search_samples = (
            None if peak_lock_search_samples is None
            else int(peak_lock_search_samples))
        self.peak_lock_search_decimate = int(peak_lock_search_decimate)
        self.peak_lock_sample_rate = float(sample_rate)

        if 'tc' not in self.marginalized_vector_priors:
            return

        tcmin, tcmax = self.marginalized_vector_priors['tc'].bounds['tc']
        tstart = tcmin - EARTH_RADIUS
        tmax = tcmax - tcmin + EARTH_RADIUS * 2.0
        num_samples = int(tmax * sample_rate)
        self.tstart = {ifo: tstart for ifo in self.data}
        self.num_samples = {ifo: num_samples for ifo in self.data}
        self.peak_lock_peak = {}

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
                self.peak_lock_peak[ifo] = float(peak_time)

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
        self.peak_lock_start = self.tstart.copy()

    def follow_peak(self, wfs):
        """ Move the locked region to wherever the peak is now

        The region is locked once, around a reference waveform, but the peak
        it was locked on moves with the parameters: a hundredth of a solar
        mass in a component mass moves it by milliseconds, which is wider
        than the region itself. Past that the region holds no peak, and the
        likelihood is wrong rather than approximate.

        A coarse pass over a wider range says where the peak went. The
        offset is measured from the region as locked rather than from
        wherever the previous call left it, so the region depends on the
        current parameters alone. One offset is found, in one detector,
        and applied to all of them, which leaves the regions as
        commensurate as they were locked: a common shift does not change
        the delays between detectors, so the sky draw sees what it expects.

        Does nothing when the marginalization points were precalculated,
        since those were drawn over the locked region, and moving away from
        them would lose them silently.
        """
        if not getattr(self, 'peak_lock_search_samples', None):
            return
        if hasattr(self, 'premarg') or not hasattr(self, 'peak_lock_start'):
            return

        ifo = (getattr(self, 'keep_ifos', None) or list(wfs))[0]
        rate = self.peak_lock_sample_rate
        if ifo not in self.peak_lock_peak:
            return
        locked = self.peak_lock_peak[ifo]

        stride = self.peak_lock_search_decimate
        delta_t = stride / rate
        num = int(self.peak_lock_search_samples / stride)
        start = locked - self.peak_lock_search_samples / rate / 2.0
        series = self.coarse_series(ifo, wfs, start, delta_t, num)
        series = abs(numpy.array(series))
        i = int(numpy.argmax(series))
        # the peak rarely sits on a coarse sample, and its neighbours say
        # where between them it does. Costs a handful of operations and
        # saves having to pad the region by a whole stride to absorb the
        # error of taking the coarse sample itself.
        if 0 < i < len(series) - 1:
            a, b, c = series[i - 1], series[i], series[i + 1]
            curve = a - 2.0 * b + c
            if curve < 0:
                i += min(0.5, max(-0.5, 0.5 * (a - c) / curve))
        # to a whole sample of the region's own grid: a fractional offset
        # would sample the series at shifted phases, which costs more than
        # the refinement gains
        shift = round((start + i * delta_t - locked) * rate) / rate

        for name in self.peak_lock_start:
            self.tstart[name] = self.peak_lock_start[name] + shift
            self.tend[name] = self.tstart[name] + self.num_samples[name] / rate

    def draw_ifos(self, snrs, peak_snr_threshold=3.0, log=True,
                  precalculate_marginalization_points=False,
                  **kwargs):
        """ Helper utility to determine which ifos we should use based on the
        reference SNR time series.

        `peak_snr_threshold` is compared against the PEAK OF THE
        MATCHED-FILTER SNR SERIES, which is noise-dependent, so which
        detectors survive varies with the noise realisation. It selects the
        detectors that drive the extrinsic PROPOSAL only, never the
        likelihood sum, so it cannot bias the result.

        The default is 3.0 rather than something safer-looking like 4.0
        because dropping a detector is not a neutral act: with only one left
        there is no time delay to invert, and `draw_sky_times` falls back to
        drawing the sky from the isotropic prior. Measured on a BNS P-P at
        fixed noise, two injections at network SNR 7-8 whose peaks were
        3.7/3.8/5.8 and 3.8/3.4/4.6 kept V1 alone at 4.0, which starved the
        extrinsic marginalization to an effective 18 and 11 samples out of
        15000 and capped efficiency at 24.6% and 20.0%. At 3.0 all three
        detectors survive, the effective sample size roughly doubles, and
        efficiency reaches 59.8% and 65.0% -- 2.43x and 3.25x. Injections
        whose detector set does not change are bit-identical, so the
        loosening costs nothing where it changes nothing.
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
            self.premarg['antenna_factors'] = self.precalc_antenna_factors

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
        return params

    def update(self, **params):
        """Move to a new point, and forget the weights that went with the
        last one.

        The weights start out equal and the draws add to them. Resetting
        here rather than whenever the parameters are read means a draw's
        weights survive being asked for the parameters afterwards, which
        is what the likelihoods do.
        """
        super().update(**params)
        self.marginalize_vector_weights = - numpy.log(self.vsamples)
    @property
    def reconstruct_stats(self):
        """The parameters being drawn inline, in the order the levels run."""
        levels = getattr(self, 'reconstruct_inline', [])
        names = []
        if 'vector' in levels:
            names += list(self.marginalize_vector_params)
        if 'distance' in levels and self.distance_marginalization:
            names.append('distance')
        if 'phase' in levels and self.marginalize_phase:
            names.append('coa_phase')
        return names

    @property
    def _extra_stats(self):
        """Adds whatever is being drawn inline to the model's own stats."""
        return super()._extra_stats + self.reconstruct_stats

    def draw_inline(self, sh_total, hh_total, vector, subset=None):
        """ Draw the marginalized parameters from the vectors the
        marginalization has just built, rather than by evaluating the
        likelihood again.

        The levels run in the order `reconstruct` does them, each conditioned
        on the point the level above drew. Conditioning is indexing here: the
        drawn vector point selects one entry of the inner products, and the
        drawn distance rescales them. `vector` is the unmarginalized loglr at
        each vector point.
        """
        rec = {}
        sh, hh = sh_total, hh_total
        loglr = None

        levels = self.reconstruct_inline
        if 'vector' in levels and self.marginalize_vector_params:
            drawn, loglr, xl = self.draw_vector(vector, subset)
            rec.update(drawn)
            if numpy.ndim(sh):
                sh, hh = sh[xl], hh[xl]

        if 'distance' in levels and self.distance_marginalization:
            dist_rescale, _ = self.distance_marginalization
            dloglr = marginalize_likelihood(
                sh, hh, phase=self.marginalize_phase,
                distance=self.distance_marginalization, skip_vector=True)
            drawn, loglr, xl = self.draw_distance(dloglr)
            rec.update(drawn)
            sh, hh = sh * dist_rescale[xl], hh * dist_rescale[xl] ** 2.0

        if 'phase' in levels and self.marginalize_phase:
            drawn, loglr, _ = self.draw_phase(sh, -0.5 * hh)
            rec.update(drawn)

        # recorded like any other stat, so whatever the sampler keeps for a
        # point carries its drawn parameters with it. The marginalized
        # loglikelihood is left alone; it is what the sampler is sampling.
        for name in self.reconstruct_stats:
            setattr(self._current_stats, name, rec[name])

    def marginalize_subset(self, vector, subset):
        """ The marginalized loglr over part of the vector samples, their
        weights renormalized to that part.
        """
        logw = self.marginalize_vector_weights
        logw = (logw[subset] if numpy.ndim(logw)
                else numpy.zeros(int(subset.sum())))
        return float(logsumexp(vector[subset] + logw) - logsumexp(logw))

    def draw_vector(self, loglr, subset=None):
        """ Draw one of the vector marginalization points, given the
        unmarginalized loglr at each. Returns the parameters there, the loglr
        there, and the index, which the next level conditions on.

        `subset` restricts the draw to part of the samples, so it can be kept
        independent of the part the likelihood was marginalized over.
        """
        logw = self.marginalize_vector_weights
        if subset is not None:
            idx = numpy.flatnonzero(subset)
            xl = idx[draw_sample(
                loglr[idx] + (logw[idx] if numpy.ndim(logw) else logw),
                rng=self.marginalize_rng)]
            rec = {k: v[xl] for k, v in self.marginalize_vector_params.items()}
            return rec, loglr[xl], xl
        xl = draw_sample(loglr + logw, rng=self.marginalize_rng)
        rec = {k: v[xl] for k, v in self.marginalize_vector_params.items()}
        return rec, loglr[xl], xl

    def draw_distance(self, loglr):
        """ Draw a distance, given the loglr at each point of the distance
        grid. Returns as `draw_vector` does.
        """
        _, weights = self.distance_marginalization
        xl = draw_sample(loglr + numpy.log(weights), rng=self.marginalize_rng)
        return {'distance': self.dist_locs[xl]}, loglr[xl], xl

    def draw_phase(self, sh, hh):
        """ Draw a coalescence phase, given the inner products before phase
        marginalization (`hh` being minus one half of it). Returns as
        `draw_vector` does.
        """
        phasev = numpy.linspace(0, numpy.pi*2.0, int(1e4))
        # This assumes that the template was conjugated in inner products
        loglr = (numpy.exp(-2.0j * phasev) * sh).real + hh
        xl = draw_sample(loglr, rng=self.marginalize_rng)
        return {'coa_phase': phasev[xl]}, loglr[xl], xl

    def reconstruct(self, rec=None, seed=None, set_loglr=None):
        """ Reconstruct the distance or vectored marginalized parameter
        of this class.
        """
        if seed:
            numpy.random.seed(seed)
            if self.marginalize_rng is not None:
                self.marginalize_rng = numpy.random.default_rng(seed)

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
            drawn, loglr, _ = self.draw_vector(get_loglr())
            rec.update(drawn)
            self.reconstruct_vector = False

        if self.distance_marginalization:
            logging.debug('Reconstruct distance')
            # call likelihood to get vector output
            self.reconstruct_distance = True
            drawn, loglr, _ = self.draw_distance(get_loglr())
            rec.update(drawn)
            self.reconstruct_distance = False

        if self.marginalize_phase:
            logging.debug('Reconstruct phase')
            self.reconstruct_phase = True
            drawn, loglr, _ = self.draw_phase(*get_loglr())
            rec.update(drawn)
            self.reconstruct_phase = False

        rec['loglr'] = loglr
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
    # geomspace: uniform in log, so the cell index is arithmetic
    log_shr, log_hhr = numpy.log(shr), numpy.log(hhr)
    dlog_shr = (log_shr[-1] - log_shr[0]) / (len(log_shr) - 1)
    dlog_hhr = (log_hhr[-1] - log_hhr[0]) / (len(log_hhr) - 1)
    coeffs = ndimage.spline_filter(lvals, order=3, mode='nearest')

    # said once, the first time it happens
    warned = [False]

    def warn_out_of_range():
        warned[0] = True
        logging.warning(
            "A likelihood evaluation asked for a signal to noise ratio "
            "outside marginalize_distance_snr_range %s; beyond it the "
            "distance marginalized likelihood is set to zero, which biases "
            "the result. Widen the range.", snr_range)

    # Out-of-range queries are set to -inf below. Clamping them to the table
    # edge instead was considered and NOT adopted: zeroing the likelihood makes
    # such samples unrecoverable, and when a large fraction of the extrinsic
    # draws fall below the SNR floor the estimate is left resting on a handful
    # of survivors. The edge value would OVERSTATE a genuinely sub-floor point
    # while staying negligible against in-range samples, so revisit this if the
    # out-of-range warning below fires on a large fraction of draws.

    def interp_wrapper(x, y, bounds_check=True):
        k = None
        # a zero-dimensional array is a single point but is not a float
        scalar = numpy.ndim(x) == 0
        if bounds_check:
            if scalar:
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

        # clipped before the log: sh can be negative without phase
        # marginalization, and the mask does not overwrite a nan. Nothing
        # then reaches the boundary, but mode must not be the default,
        # which reads outside the grid as zero.
        index = numpy.empty((2, numpy.size(x)))
        index[0] = (numpy.log(numpy.clip(x, shr_min, shr_max))
                    - log_shr[0]) / dlog_shr
        index[1] = (numpy.log(numpy.clip(y, hhr_min, hhr_max))
                    - log_hhr[0]) / dlog_hhr
        v = ndimage.map_coordinates(coeffs, index, order=3, mode='nearest',
                                    prefilter=False)
        if k is not None:
            v[k] = -numpy.inf
        # marginalize_likelihood tells a point from a vector by asking for a
        # float, and a length-one array goes through the vector-marginalization
        # weight instead -- subtracting log(marginalize_vector_samples) from
        # the distance-marginalized likelihood. Only bites distance
        # marginalization without an accompanying time or sky marginalization;
        # with one, x is a genuine vector. map_coordinates always returns an
        # array, hence v[0].
        if scalar:
            return float(v[0])
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
                           return_vector=False,
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
    return_vector: bool, False
        Also return the unmarginalized loglr at each point. That is the
        distribution reconstruction draws from, so returning it saves
        evaluating the likelihood again to get it back.

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

    vector = vloglr

    # Do brute-force marginalization if loglr is a vector
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    elif not skip_vector:
        vloglr = float(logsumexp(vloglr, b=numpy.exp(logw)))

    if return_peak:
        return vloglr, maxv, maxl
    if return_vector:
        return vloglr, vector
    return vloglr
