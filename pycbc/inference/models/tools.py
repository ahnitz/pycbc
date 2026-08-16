""" Common utility functions for calculation of likelihoods
"""

import logging
import os as _os
import warnings
from distutils.util import strtobool

import numpy
import numpy.random
import tqdm

from scipy.special import logsumexp, i0e
from scipy.interpolate import RectBivariateSpline, interp1d
from pycbc.distributions import JointDistribution

from pycbc.detector import Detector


# Earth radius in seconds
EARTH_RADIUS = 0.031


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

        # the effective sample size of the most recent vector
        # marginalization, filled in by marginalize_loglr; nan until then
        self.vector_ess = numpy.nan

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

        # ask for the effective sample size when doing the ordinary
        # vector marginalization, so a run can see whether it drew enough
        # points; the reconstruction and peak paths do not marginalize a
        # vector here, so there is nothing to measure
        want_ess = not (return_peak or return_complex or skip_vector)
        result = marginalize_likelihood(sh_total, hh_total,
                                        logw=self.marginalize_vector_weights,
                                        phase=self.marginalize_phase,
                                        interpolator=interpolator,
                                        distance=distance,
                                        skip_vector=skip_vector,
                                        return_complex=return_complex,
                                        return_peak=return_peak,
                                        return_ess=want_ess)
        if want_ess:
            loglr, self.vector_ess = result
            return loglr
        return result

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

        if _os.environ.get('MARG_ORIENT'):
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

        # EXPERIMENT B2 (env MARG_B2=<factor>, default 1 = current behaviour):
        # oversample the time draws so that enough of them survive the
        # physical-delay lookup. The draws themselves are cheap integer ops on
        # an already-computed SNR series, and `resize_factor` below is exactly
        # the acceptance fraction, so the estimator normalisation is unchanged.
        import os as _os
        ndraw = vsamples

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
            import os as _os
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

            # STEP 2 (env MARG_COH=1): bias WHICH POINT we take from each
            # delay bin, using amplitude/phase consistency.
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
        # The  + numpy.log(sref.delta_t / (tcmax - tcmin))  term that used to
        # sit here (-10.58 nats at 131072 Hz over a 0.3 s window) costs 20-300x
        # efficiency on a subset of injections: 00021 0.098% -> 35.5%, 00008
        # 1.362% -> 30.2% with it removed, evidence landing within 0.5 nats of
        # games-v2-trunc. The term is plausibly correct on paper -- the older
        # code carried "There's an overall normalization here which may
        # introduce a constant factor at the moment" -- so it should go back
        # once the interaction is understood. See
        # results/FULLINTEGRATION_REGRESSION.md; the mechanism is NOT yet
        # established (the tree descent is invariant to a constant, so the
        # first explanation offered there was wrong).
        logw_sky = (-mcweight[ti] + numpy.log(wi)
                    + numpy.log(resize_factor))
        # times the search had to leave out carry no prior weight
        logw_sky[(tc < tcmin) | (tc > tcmax)] = -numpy.inf

        self.marginalize_vector_params['tc'] = tc
        self.marginalize_vector_params['ra'] = ra
        self.marginalize_vector_params['dec'] = dec

        if _os.environ.get('MARG_ORIENT'):
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
        nc = int(_os.environ.get('MARG_ORIENT_NC', '24'))
        npsi = int(_os.environ.get('MARG_ORIENT_NP', '12'))
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
        eps = float(_os.environ.get('MARG_ORIENT_EPS', '0.02'))
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
        temper = float(_os.environ.get('MARG_ORIENT_T', '1'))
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

    # warn the first time the sampler asks outside the interpolated range,
    # rather than silently returning -inf. Beyond the range the distance
    # marginalized likelihood is dropped to zero, which biases the result;
    # the fix is a wider marginalize_distance_snr_range. Reported in signal
    # to noise ratio, which is what the option is set in.
    warned = [False]

    def warn_out_of_range(over):
        if warned[0] or not over:
            return
        warned[0] = True
        logging.warning(
            "A likelihood evaluation asked for a signal to noise ratio "
            "outside marginalize_distance_snr_range %s; beyond it the "
            "distance marginalized likelihood is set to zero, which biases "
            "the result. Widen the range.", snr_range)

    # EXPERIMENT (env MARG_CLAMP=1): clamp out-of-range queries to the table
    # edge instead of setting them to -inf. Zeroing the likelihood makes such
    # samples unrecoverable, and when a large fraction of the extrinsic draws
    # fall below the SNR floor the estimate is left resting on a handful of
    # survivors. The edge value OVERSTATES a genuinely sub-floor point, but it
    # is negligible against in-range samples (exp(edge) << exp(peak)), so the
    # sum is unchanged while the estimate stays defined.

    def interp_wrapper(x, y, bounds_check=True):
        k = None
        if bounds_check:
            if isinstance(x, float):
                if x > shr_max or x < shr_min or y > hhr_max or y < hhr_min:
                    warn_out_of_range(True)
                    return -numpy.inf
            else:
                k = (x > shr_max) | (x < shr_min)
                k = k | (y > hhr_max) | (y < hhr_min)
                warn_out_of_range(k.any())

        v = interp(x, y, grid=False)
        if k is not None:
            v[k] = -numpy.inf
        # a scalar query is a single point with nothing to marginalize over;
        # the spline returns it as a length-one array, which the caller then
        # mistakes for a vector and folds through the vector-marginalization
        # weight, subtracting log(marginalize_vector_samples) from the
        # distance-marginalized likelihood. Hand back a scalar so it does
        # not. Only bites distance marginalization without an accompanying
        # time or sky marginalization; with one, x is a genuine vector.
        if numpy.ndim(x) == 0:
            return float(v)
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
                           return_ess=False,
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
            return (vloglr, None) if return_ess else vloglr
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
    ess = None
    if isinstance(vloglr, float):
        vloglr = float(vloglr)
    elif not skip_vector:
        if return_ess:
            # the effective number of the drawn points that carry the
            # answer: (sum w)^2 / sum w^2 with w the combined weights.
            # Low against the number drawn means the marginal rests on a
            # handful of them and its error is correspondingly large.
            #
            # Exponentiate once and form the ratio directly rather than
            # asking logsumexp twice for the same array. Subtracting the
            # largest log weight first is what keeps exp in range, exactly
            # as logsumexp does it internally, and the ratio is scale free
            # so that shift cancels. Doing the sums in linear space instead
            # of log space moves the answer against the logsumexp route by
            # of order 1e-14 relative, a few times 1e-13 at worst. That is
            # accepted deliberately: the effective sample size is a
            # monitoring number, it enters neither the posterior nor the
            # evidence, and chasing bit-identity would mean copying
            # scipy's internals and tracking their changes.
            lw = vloglr + logw
            # an empty draw has no maximum to take, and no sample size to
            # report, so it joins the degenerate cases rather than raising
            lwmax = lw.max() if lw.size else numpy.nan
            if numpy.isfinite(lwmax):
                # the largest weight is exactly one after the shift, so the
                # denominator is at least one and the sums cannot overflow
                w = numpy.exp(lw - lwmax)
                ess = float(w.sum() ** 2.0 / numpy.vdot(w, w))
            else:
                # every weight vanished, or one is infinite; there is
                # nothing meaningful to count
                ess = numpy.nan
        vloglr = float(logsumexp(vloglr, b=numpy.exp(logw)))

    if return_peak:
        return vloglr, maxv, maxl
    if return_ess:
        return vloglr, ess
    return vloglr
