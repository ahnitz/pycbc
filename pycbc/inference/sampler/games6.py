""" Importance sampling from a precomputed map of the parameter space,
combined by strata with an adaptive generative proposal.

Draws come from two kinds of proposal:

  pool  a discrete pool of precomputed prior points, taken from the
        leaves of the map that survive a likelihood cut. Efficient, but
        finite -- a leaf can only supply the points it holds.
  gen   a Gaussian mixture fitted to the accumulated posterior samples,
        with a separate covariance per centre estimated from its k
        nearest neighbours. Unlimited, and takes over once the pool
        exhausts.

The two are not combined into a single density: the pool proposal puts
mass on atoms and the generative one has a Lebesgue density, so the
measures are mutually singular. Instead each proposal -- and each
generative refit -- is its own stratum with its own self-normalised
estimate, and the strata are combined as

    I = sum_s beta_s I_s,    sum_s beta_s = 1

with beta proportional to the stratum's ESS. That is the
inverse-variance choice, and it makes the combined ESS the sum of the
strata's, so a stratum that is doing badly contributes little.

Every sample keeps the weight computed against the proposal in force
when it was drawn; weights are never recomputed against a later
proposal. Out-of-prior draws are dropped before any likelihood call.

The map is built by ``pycbc_inference_build_map``.
"""
import logging
import os
import tqdm
import h5py
import numpy
import numpy.random
from scipy.special import logsumexp
from scipy.stats import chi2
from pycbc.io import FieldArray

from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler
from .games import call_likelihood, call_tile_likelihood, OutOfSamples


def _pareto_k(logw):
    """ Tail index of the largest importance weights, by Hill's estimator.

    The statistical error on an importance-sampling estimate is only
    meaningful if the weights have finite variance, and importance weights
    frequently do not: a single draw landing where q is tiny and pi is not
    dominates the sum, and then the sample variance says nothing about the
    true one. The tail index says which regime you are in --

        k < 0.5   finite variance, the CLT applies and sigma is trustworthy
        0.5 - 0.7 marginal
        k > 0.7   sigma understates, treat the value as indicative
        k > 1     infinite mean, the estimate itself is not to be trusted

    Hill's estimator rather than a fitted generalised Pareto: it needs no
    optimisation, works in log space so nothing overflows, and this is a
    diagnostic rather than a correction.
    """
    w = numpy.sort(numpy.asarray(logw, float))
    w = w[numpy.isfinite(w)]
    n = len(w)
    if n < 50:
        return float('nan')
    m = int(min(n // 5, 3.0 * numpy.sqrt(n)))
    if m < 10:
        return float('nan')
    tail = w[-m:]
    return float(numpy.mean(tail - w[-m - 1]))


class LocalCovarianceProposal:
    """ Gaussian mixture with a separate covariance per centre.

    Each centre's kernel covariance is estimated from its `k` nearest
    neighbours, so the kernel follows the local orientation and thickness
    of the posterior rather than a single global shape.

    Centres, kernels and draws all live in globally whitened coordinates,
    where every parameter is O(1); `logpdf` returns the density in the
    original coordinates.
    """

    def __init__(self, centres, rng, k=40, scale=1.0, weights=None,
                 ridge=1e-3, bounds=None, nz=150):
        from scipy.spatial import cKDTree
        self.rng = rng
        self.centres = numpy.atleast_2d(centres)
        self.n, self.dim = self.centres.shape
        if weights is None:
            weights = numpy.full(self.n, 1.0 / self.n)
        self.weights = numpy.asarray(weights, float)
        self.weights /= self.weights.sum()
        self.logw = numpy.log(numpy.maximum(self.weights, 1e-300))

        g = numpy.cov(self.centres.T)
        ev, V = numpy.linalg.eigh(g)
        ev = numpy.maximum(ev, 1e-12 * max(ev.max(), 1e-300))
        self._W = V / numpy.sqrt(ev)          # x -> z
        self._Winv = (V * numpy.sqrt(ev)).T   # z -> x
        self._logdetW = -0.5 * numpy.log(ev).sum()

        cz = self.centres @ self._W
        self._cz = cz
        tree = cKDTree(cz)
        kk = min(max(k, self.dim + 2), self.n)
        _, idx = tree.query(cz, k=kk)

        f = kk ** (-1.0 / (self.dim + 4)) * scale
        self.chol = numpy.empty((self.n, self.dim, self.dim))
        self.inv = numpy.empty((self.n, self.dim, self.dim))
        self.lognorm = numpy.empty(self.n)
        eye = numpy.eye(self.dim)
        for i in range(self.n):
            c = numpy.cov(cz[idx[i]].T) * f ** 2
            # dimensionless in whitened coordinates, so the same floor
            # regularises every direction equally
            c = c + ridge * f ** 2 * eye
            try:
                L = numpy.linalg.cholesky(c)
            except numpy.linalg.LinAlgError:
                L = numpy.linalg.cholesky(c + 10 * ridge * f ** 2 * eye)
            self.chol[i] = L
            self.inv[i] = numpy.linalg.inv(L @ L.T)
            self.lognorm[i] = -0.5 * (self.dim * numpy.log(2 * numpy.pi)
                                      + 2.0 * numpy.log(numpy.diag(L)).sum())

        # Truncate every kernel to the prior box and renormalise it.
        #
        # Without this the mixture leaks: each parameter is only ~15% too
        # wide, but containment is a product over parameters, so at ten of
        # them well over half the draws land outside, and the density that
        # survives is DEPRESSED just inside every wall -- which is where
        # pi/q then explodes. Truncating the MIXTURE would only rescale it
        # by a global constant, identical to rejecting those draws;
        # truncating each COMPONENT redistributes its mass back inside.
        #
        # It reduces the wall deficit rather than removing it. On an
        # exactly uniform target the mixture sits at int phi(t)/Phi(t) dt
        # = 0.69 of the target AT the wall, against 0.50 untruncated
        # (measured 0.65 in six dimensions). Per dimension that is pi/q
        # 1.45 instead of 2.0, so at a ten dimensional corner 42 instead
        # of 1024 -- most of the problem, but not all of it.
        #
        self._lo = self._hi = None
        self.logz = numpy.zeros(self.n)
        if bounds is not None:
            self._lo, self._hi = bounds
            # Every component's mass inside the box, in chunks: at ten
            # bounded parameters essentially every kernel touches a wall,
            # so there is nothing to skip and the loop has to be vectorised
            # rather than avoided.
            #
            # The draws are INDEPENDENT per component, deliberately. Sharing
            # one set across components is cheaper and is the usual variance
            # reduction trick, but it correlates every frac error, so they
            # no longer average out across the mixture and the density is
            # left mis-normalised by a percent or two however large nz is.
            # Independent draws let that error fall as 1/sqrt(nz n) instead.
            frac = numpy.empty(self.n)
            step = 256
            for a in range(0, self.n, step):
                b = min(a + step, self.n)
                e = self.rng.normal(size=(b - a, nz, self.dim))
                # (chunk, nz, dim)
                pts = (self._cz[a:b, None, :]
                       + numpy.einsum('cij,cnj->cni', self.chol[a:b], e)) \
                    @ self._Winv
                frac[a:b] = numpy.all((pts >= self._lo) & (pts <= self._hi),
                                      axis=2).mean(axis=1)
            # frac is a Monte Carlo estimate and the density divides by it,
            # so Jensen biases 1/frac high by (1-Z)/(nz Z) -- a few percent
            # at nz=150. That cancels in self-normalised importance
            # sampling but not in the evidence, where a few percent is
            # larger than the estimator's demonstrated accuracy. Correcting
            # frac itself removes the leading term.
            # frac itself removes the leading term -- but only where the
            # expansion holds. It is first order in (1-Z)/(nz Z), so it
            # needs nz Z large, and a component that keeps little of its
            # mass inside does not have that: at ten parameters many sit
            # near nz Z ~ 8, where correcting makes things worse. Measured
            # by a bridge estimator, q then integrated to 0.955 rather
            # than 1.
            #
            # Correcting is only safe where nz Z is large, so leave the
            # small ones alone rather than correcting them wrongly. The
            # AGGREGATE error is what matters and it falls as 1/sqrt(n),
            # which holds only because the draws above are independent per
            # component: measured, the integral of q over the box is 0.96
            # at 2000 components and 0.99 at 8000. Anything wanting better
            # than that -- evidence, where a percent is 0.01 nats -- should
            # spend draws on the poorly resolved components, since for
            # posteriors the error cancels in the self-normalised estimate.
            frac = numpy.maximum(frac, 1.0 / nz) + (1.0 - frac) / nz
            self.logz = numpy.log(frac)
            logging.info('kernel truncation: %i of %i components touch a '
                         'boundary', int((frac < 1.0).sum()), self.n)

    def sample(self, size):
        if self._lo is None:
            idx = self.rng.choice(self.n, size=size, p=self.weights)
            z0 = self.rng.normal(size=(size, self.dim))
            z = self._cz[idx] + numpy.einsum('nij,nj->ni', self.chol[idx], z0)
            return z @ self._Winv
        out = numpy.empty((size, self.dim))
        got = 0
        while got < size:
            idx = self.rng.choice(self.n, size=size, p=self.weights)
            z0 = self.rng.normal(size=(size, self.dim))
            x = (self._cz[idx]
                 + numpy.einsum('nij,nj->ni', self.chol[idx], z0)) @ self._Winv
            ok = numpy.all((x >= self._lo) & (x <= self._hi), axis=1)
            take = min(int(ok.sum()), size - got)
            out[got:got + take] = x[ok][:take]
            got += take
        return out

    def logpdf(self, x):
        z = x @ self._W
        out = numpy.full(len(z), -numpy.inf)
        for i in range(self.n):
            d = z - self._cz[i]
            m = numpy.einsum('ij,jk,ik->i', d, self.inv[i], d)
            out = numpy.logaddexp(out, -0.5 * m + self.lognorm[i]
                                  + self.logw[i] - self.logz[i])
        # change of variables from whitened z back to physical x
        out = out + self._logdetW
        if self._lo is not None:
            inside = numpy.all((x >= self._lo) & (x <= self._hi), axis=1)
            out = numpy.where(inside, out, -numpy.inf)
        return out


class FactorizedProposal:
    """ A KDE over some parameters, the prior over the rest.

    Directions the likelihood does not constrain are better left alone
    than modelled: the posterior there IS the prior, and a Gaussian
    mixture asked to represent a flat bounded density returns something
    slightly too wide with rounded shoulders, which is pure weight
    variance for no information. Measured on a BNS with four unconstrained
    in-plane spin components, modelling them cost a factor of ~3.7 in
    efficiency against deferring them.

    The deferred parameters supply exactly their prior, and the target
    carries the same factor, so the two cancel in a self-normalised
    POSTERIOR. They do not cancel in the evidence, which multiplies by the
    full prior explicitly, so `logpdf` must report the density of what it
    actually samples: the KDE over the modelled subset TIMES the prior over
    the deferred one. Returning the KDE alone made log Z high by the
    deferred prior density -- measured as 6.4 nats with four deferred
    in-plane spins, 4 ln(5), while the aligned cases with nothing deferred
    were exact.

    All of this holds only while the deferred parameters are
    prior-INDEPENDENT of the modelled ones.
    """

    def __init__(self, inner, model_idx, defer_idx, prior_draw, dim,
                 defer_logpdf=None):
        self.inner = inner
        self.model_idx = numpy.asarray(model_idx, dtype=int)
        self.defer_idx = numpy.asarray(defer_idx, dtype=int)
        self.prior_draw = prior_draw
        self.defer_logpdf = defer_logpdf
        self.dim = int(dim)
        self.n = inner.n
        self.centres = inner.centres
        self.weights = inner.weights

    def sample(self, size):
        out = numpy.empty((size, self.dim))
        out[:, self.model_idx] = self.inner.sample(size)
        out[:, self.defer_idx] = self.prior_draw(size)
        return out

    def logpdf(self, x):
        # The density of what `sample` produces: KDE over the modelled
        # subset times the prior over the deferred one. The deferred factor
        # cancels against the target in the posterior either way, but the
        # evidence needs it present.
        lp = self.inner.logpdf(x[:, self.model_idx])
        if self.defer_logpdf is not None:
            lp = lp + self.defer_logpdf(x[:, self.defer_idx])
        return lp


class GameSampler6(DummySampler):
    """ Stratified importance sampler over a precomputed parameter map.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Flat map file, as written by pycbc_inference_build_map.
    treefile : str, optional
        Hierarchical tree file; tiles without an entry fall back to the flat
        pool.
    loglr_region : int
        Only use regions within this loglr of the best node.
    target_likelihood_calls : int
        Likelihood calls to aim for per round.
    rounds : int
        Number of rounds.
    gen_start_round : int
        First round (1-based) at which the generative proposal may be used. It
        needs enough accumulated posterior samples to fit.
    gen_min_ess : float
        Minimum accumulated ESS before fitting the generative proposal.
    gen_max_centres : int
        Cap on mixture centres taken from the posterior samples.
    gen_refit : int
        Refit the generative proposal every this many rounds (1 = every round,
        0 = fit once and freeze). Refitting is safe because each sample retains
        the weight computed against the proposal in force when it was drawn.
    gen_bandwidth : float
        Extra scale on the kernel covariance. With bounded kernels the
        default 1.45 is chosen by worst case over seven SNR rungs, not by
        best case: the optimum per rung ranges from 1.15 to 1.8, and below
        1.45 the result stops being reproducible -- one rung gave 1.00 and
        64.5 on two seeds at 1.15. Truncation renormalises a cut kernel and
        so concentrates it, which is why bounded kernels want MORE
        inflation, not less.
    gen_local_k : int
        Neighbours used to estimate each centre's own kernel covariance.
    tree_accept_fraction : float
        Stop refining a tree node when at least this fraction of its
        children clear the likelihood cut, and take their subtrees whole
        instead. The descent only earns its calls by pruning; below SNR ~10
        it prunes nothing and enumerates the map. 0 disables.
    gen_defer : float
        Loading below which a parameter is ELIGIBLE to be handed to its exact
        uniform prior instead of being modelled by the kernel. The loading is
        a direction cosine -- |V[i,j]| for eigenvector j of the posterior
        covariance in prior units -- so 0.2 means the parameter axis sits
        within 12 degrees of perpendicular to everything the likelihood
        constrained. ELIGIBILITY ONLY: whether an eligible parameter is
        actually deferred is decided on measured efficiency, see
        `_defer_dimensions`. 0 disables deferral entirely.
    gen_defer_ratio : float
        Only attempt deferral when the generative round's efficiency has
        fallen below the pool stratum's by this factor. Deferral is a remedy
        for a kernel that is measurably doing badly; a kernel already beating
        the pool needs no remedy and must not be given one speculatively.
    gen_defer_patience : int
        Consecutive poor rounds required before attempting deferral, and the
        window over which "recovering on its own" is judged. The kernel is
        legitimately poor in its first rounds while it has few centres, and
        recovers unaided; acting on that transient is what has to be avoided.
    gen_defer_margin : float
        A trial deferral is kept only if it improves efficiency by at least
        this factor. Otherwise it is reverted and no further parameter is
        tried, so a wrong guess costs one round rather than the whole run.
    gen_switch_ess : float
        Accumulated ESS at which the generative proposal takes over from the
        pool, WITHOUT waiting for the pool to exhaust. The kernel beat the
        pool at every seed size measured -- by 2.4x when the seed was under
        2000 ESS -- and converged to 68-79% regardless of it, so the handover
        is worth making early. 0 keeps the old behaviour of waiting for
        exhaustion.
    gen_switch_patience : int
        Consecutive rounds the kernel must run BELOW the pool's most recent
        round rate before handing back. One bad round is noise.
    gen_switch_backoff : float
        After handing back, the next attempt waits until the accumulated ESS
        has grown by this factor. Multiplicative, so the number of retries is
        logarithmic in the budget rather than linear. Nothing here is sticky:
        an early switch that fails is retried later with more to fit, which
        is deliberately the opposite of the `gen_defer` latch.
    gen_truncate : int
        Truncate each kernel to the prior box and renormalise it. Requires
        every sampled parameter to have a bounded prior; inert otherwise.
    """
    name = 'games6'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None, treefile=None, dagfile=None,
                 dag_cull=0,
                 loglr_region=0,
                 loglr_coverage=1e-4,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 gen_start_round=1,
                 gen_min_ess=100,
                 gen_max_centres=20000,
                 gen_refit=1,
                 gen_bandwidth=1.45,
                 gen_local_k=160,
                 gen_truncate=1,
                 gen_defer=0.2,
                 gen_defer_ratio=2.0,
                 gen_defer_patience=2,
                 gen_defer_margin=1.1,
                 gen_switch_ess=0,
                 gen_switch_patience=2,
                 gen_switch_backoff=2.0,
                 tree_start_level=0,
                 tree_accept_fraction=0.9,
                 min_active_points=100,
                 tile_point_cap=200000,
                 gen_seed_leaves=0,
                 gen_seed_range=4.0,
                 leaf_weight_samples=0,
                 pool_adapt_mean=0,
                 **kwargs):
        super().__init__(model, *args)

        self.meta = {}
        self.mapfile = mapfile
        self.dagfile = dagfile
        # Opt-in, like gen_switch_ess: a cut that prunes is a behaviour
        # change, and one that is on by default cannot be validated against
        # the version without it.
        self.dag_cull = int(dag_cull)
        self.treefile = treefile
        self.rounds = int(rounds)
        self.dmap = {}
        self.draw = {}

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        # A d-dimensional posterior encloses 1-eps of its mass within a
        # loglr drop of chi2_d(1-eps)/2, so the region follows from the
        # coverage asked for and the number of parameters. Set
        # loglr_region explicitly to override.
        if float(loglr_region) > 0:
            self.loglr_region = float(loglr_region)
        else:
            d = len(model.variable_params)
            self.loglr_region = float(
                chi2.ppf(1.0 - float(loglr_coverage), d) / 2.0)
            logging.info('loglr_region=%.1f from %i parameters and '
                         'coverage 1-%g', self.loglr_region, d,
                         float(loglr_coverage))
        self.gen_start_round = int(gen_start_round)
        self.gen_min_ess = float(gen_min_ess)
        self.gen_max_centres = int(gen_max_centres)
        self.gen_refit = int(gen_refit)
        self.gen_bandwidth = float(gen_bandwidth)
        self.gen_local_k = int(gen_local_k)
        self.gen_truncate = int(gen_truncate)
        # Loading on informative directions below which a parameter is
        # handed to the prior exactly. 0 disables.
        self.gen_defer = float(gen_defer)
        self.gen_defer_ratio = float(gen_defer_ratio)
        self.gen_defer_patience = int(gen_defer_patience)
        self.gen_defer_margin = float(gen_defer_margin)
        self.gen_switch_ess = float(gen_switch_ess)
        self.gen_switch_patience = int(gen_switch_patience)
        self.gen_switch_backoff = float(gen_switch_backoff)
        self._deferred = set()
        self._defer_trial = None
        self._defer_eff_before = None
        self._defer_frozen = False
        self._pool_eff = None
        self._gen_eff = []
        self._switched = False
        self._switch_at = float(gen_switch_ess)
        self._switch_bad = 0
        self.tree_start_level = int(tree_start_level)
        self.tree_accept_fraction = float(tree_accept_fraction)
        self.tree_shortcut = 0
        self._dropped = []
        # Fall back to start-level tiles only if the cut leaves fewer
        # points than this in total. A pool that saturates is survivable,
        # since its budget passes to the generative stratum; a pool with
        # almost nothing in it is not.
        self.min_active_points = int(min_active_points)
        # Cap on points kept per start-level tile on the coarse path. A
        # tile's subtree can hold the whole map. Its true population is
        # still used as its prior mass, so the cap changes only what is
        # available to draw.
        self.tile_point_cap = int(tile_point_cap)
        self.gen_seed_leaves = int(gen_seed_leaves)
        self.leaf_weight_samples = int(leaf_weight_samples)
        self.pool_adapt_mean = int(pool_adapt_mean)
        self.gen_seed_range = float(gen_seed_range)
        self.leaf_loglr = {}
        self.extra_params = []
        self.extra_reference = {}
        # seeded from the global RNG that pycbc_inference --seed sets, so
        # generative draws are reproducible alongside the pool shuffles
        self._rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 31 - 1))
        self.ncalls = 0
        self.tree_ncalls = 0
        self.stratum_calls = {}

    # ---------------- discrete pool stratum ---------------------------------

    def draw_samples_from_bin(self, i, size):
        if i not in self.draw:
            self.draw[i] = numpy.arange(0, len(self.dmap[i]))
        if size > len(self.draw[i]):
            raise OutOfSamples
        numpy.random.shuffle(self.draw[i])
        selected = self.draw[i][:size]
        self.draw[i] = self.draw[i][size:]
        return self.dmap[i][selected]

    def pool_round(self, bin_weight, node_idx, lengths, budget,
                   prior_mass=None):
        """ One round of the discrete-pool draw.

        Returns (psamp, loglr, unnormalised logw, bin_id), or raises
        OutOfSamples when no bin can supply another point.

        `lengths` is how many points a bin holds, and so caps how many
        can be drawn from it. `prior_mass` is how much prior mass it
        represents, which is what enters the importance weight. The two
        are equal when a bin stores every point assigned to it, and
        differ when it stores a subsample of a larger population, as the
        coarse-tile path does. Defaults to `lengths`.
        """
        if prior_mass is None:
            prior_mass = lengths
        drawcount = (bin_weight * budget).astype(int)
        dorder = bin_weight.argsort()[::-1]
        remainder = 0
        for i in dorder:
            bincount, binlen = drawcount[i], lengths[i]
            if bincount > binlen:
                drawcount[i] = binlen
                remainder += bincount - binlen
            elif bincount < binlen:
                asize = min(binlen - bincount, remainder)
                drawcount[i] += asize
                remainder -= asize
        with numpy.errstate(divide='ignore', invalid='ignore'):
            drawweight = bin_weight / drawcount
        total = drawcount.sum()
        if total <= 0:
            raise OutOfSamples
        psamp = FieldArray(total, dtype=self.samp_dtype)
        pweight = numpy.zeros(total)
        bin_id = numpy.zeros(total, dtype=int)
        j = 0
        for i, (c, w) in enumerate(zip(drawcount, drawweight)):
            bdraw = self.draw_samples_from_bin(node_idx[i], c)
            psamp[j:j+len(bdraw)] = self._as_samples(bdraw)
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        # The map's points were drawn from the MAP's prior, which need not
        # be the model's: parameters the map does not carry are filled from
        # the prior here, and a constraint coupling them to the mapped ones
        # (|chi| <= a_max across all three spin components, say) can reject
        # the combination. Drop those before paying for the likelihood, as
        # the generative round does. A no-op when the two priors agree.
        if self.extra_params:
            keep = self._in_prior(
                numpy.column_stack([psamp[p_] for p_ in
                                    self.model.variable_params]),
                list(self.model.variable_params))
            if not keep.all():
                logging.info('%i of %i pool draws fell outside the model '
                             'prior', int((~keep).sum()), len(keep))
                psamp, pweight, bin_id = psamp[keep], pweight[keep], \
                    bin_id[keep]
                if not len(psamp):
                    raise OutOfSamples

        loglr = self._likelihoods(psamp, stratum='pool')
        logw = loglr + numpy.log(prior_mass[bin_id]) - pweight
        ratio = self._log_prior_ratio(psamp)
        if not numpy.isscalar(ratio):
            logw = logw + ratio
        # For the evidence, sum(w) over the whole map's prior draws is
        # M_total * Z, so the normalisation is the map's total size rather
        # than anything about this round.
        self._logz_terms.setdefault('pool', []).append(
            (logsumexp(logw), self.map_size))
        return psamp, loglr, logw, bin_id

    # ---------------- generative stratum -----------------------------------

    def gen_round(self, proposal, budget):
        """ One round from the generative proposal. Out-of-prior draws are
        dropped before any likelihood call.
        """
        names = list(self.model.variable_params)
        xall = proposal.sample(budget)
        keep = self._in_prior(xall, names)
        if keep.sum() == 0:
            return None
        x = xall[keep]
        logq = proposal.logpdf(x)

        psamp = FieldArray(len(x), dtype=self.samp_dtype)
        for k, p in enumerate(names):
            psamp[p] = x[:, k]
        loglr = self._likelihoods(psamp, stratum='gen')
        # The prior cancels in a self-normalised posterior estimate but NOT
        # in the evidence, so carry it explicitly.
        lpri = numpy.array([self.model.prior_distribution(
            **{p_: x[i, k_] for k_, p_ in enumerate(names)})
            for i in range(len(x))])
        logw = loglr + lpri - logq
        # Z = E_q[L pi / q], so the denominator is what was DRAWN, not what
        # survived the prior cut: rejected draws carry zero weight, they do
        # not disappear from the average.
        good = numpy.isfinite(logw)
        if good.any():
            self._logz_terms.setdefault('gen', []).append(
                (logsumexp(logw[good]), int(budget)))
        return psamp, loglr, logw

    def _in_prior(self, x, names):
        ok = numpy.zeros(len(x), dtype=bool)
        for k in range(len(x)):
            pset = {p: x[k, i] for i, p in enumerate(names)}
            ok[k] = numpy.isfinite(self.model.prior_distribution(**pset))
        return ok

    def _likelihoods(self, psamp, stratum=None):
        args = [{p: psamp[p][i] for p in self.model.variable_params}
                for i in range(len(psamp))]
        vals = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                              total=len(args)))
        self.ncalls += len(args)
        if stratum is not None:
            self.stratum_calls[stratum] = \
                self.stratum_calls.get(stratum, 0) + len(args)
        return numpy.array(vals)

    # ---------------- driver ----------------------------------------------

    @staticmethod
    def _ess(logw):
        if logw is None or len(logw) == 0:
            return 0.0
        w = numpy.exp(logw - logsumexp(logw))
        return float(1.0 / (w ** 2).sum())

    def _nodes_at_depth(self, group, depth):
        """ Groups exactly `depth` levels below `group`; a branch that
        bottoms out earlier contributes its leaf, so no points are lost.
        """
        if depth <= 0 or bool(group.attrs['is_leaf']):
            return [group]
        out = []
        for c in group['children']:
            out.extend(self._nodes_at_depth(group['children'][c], depth - 1))
        return out

    def _dag_cut(self):
        """ Reach the cut by descending a multi-parent structure.

        The flat cut evaluates a likelihood at every tile. Here the levels
        are template banks, coarsest first, and each node lists EVERY parent
        within reach of it rather than the one it was filed under. A level is
        evaluated, pruned, and only the children of survivors are evaluated
        next, so most tiles are never touched.

        Two things make the pruning sound, and both were got wrong before:

        The reference is the best representative seen AT THIS LEVEL. The best
        leaf is not known until the descent ends, so comparing against it is
        not something the algorithm can do; comparing against the level's own
        maximum also means the maximum always survives, so a level can never
        come back empty.

        The allowance is rho^2 sqrt(1 - M^2) / 2, the most the loglr can RISE
        from a representative to a point within match M of it: with
        loglr = rho^2 cos^2(theta)/2 that maximum is rho^2 sin(acos(M))/2.
        The obvious rho^2 (1 - M^2)/2 is the drop from a perfect match down
        to M, which is smaller by sqrt(1 - M^2) and prunes nodes holding the
        peak.

        A node survives if ANY of its parents did. That is what the extra
        links buy: a node is judged by its best parent, not by the accident
        of which one it hangs from, and it is what stopped the strict tree
        losing 178 of 1363 tiles across 18 of 39 injections.

        Returns (loglrs, nevaluated) with -inf at every tile never reached,
        so the caller's region cut sees only what was actually evaluated.
        """
        with h5py.File(self.dagfile, 'r') as dag:
            nlev = len([k for k in dag if k.startswith('level_')])
            thr, params, links = [], [], [None]
            for d in range(nlev):
                g = dag['level_%i' % d]
                # descent_radius, not threshold. The builder links a node to
                # every parent within a CUMULATIVE cone, deliberately wider
                # than the bank's own match level, and that cone is what
                # bounds how far a node's descendants can lie. Using
                # threshold here understates the reach and prunes nodes
                # holding the peak: at level 1 it gave 0.97 where the real
                # reach is 0.918, an allowance of 52 nats against a measured
                # rise of 200.
                thr.append(float(g.attrs.get('descent_radius',
                                             g.attrs['threshold'])))
                params.append({p: g['param_%s' % p][:]
                               for p in self.dtype.names})
                if d:
                    off = g['parent_offset'][:]
                    par = g['parents'][:]
                    links.append((par, off))
        logging.info('dag cut: %s nodes per level, descent radii %s',
                     [len(next(iter(p.values()))) for p in params], thr)

        alive, total, rho2 = None, 0, None
        last_loglr = None
        lognl = float(self.model.lognl)
        perlev = []
        for d in range(nlev):
            n = len(next(iter(params[d].values())))
            if alive is None:
                idx = numpy.arange(n)
            else:
                par, off = links[d]
                idx = numpy.array(
                    [j for j in range(n)
                     if alive.intersection(par[off[j]:off[j + 1]].tolist())],
                    dtype=int)
            if not len(idx):
                logging.warning('dag cut: level %i reached no nodes', d)
                break
            args = [dict(getattr(self, 'extra_reference', {}),
                         **{p: float(params[d][p][j])
                            for p in self.dtype.names}) for j in idx]
            vals = numpy.array(list(self.pool.imap(call_tile_likelihood,
                                                   args)))
            total += len(idx)
            here = numpy.nanmax(vals)
            # call_tile_likelihood returns a loglikelihood, so the noise term
            # has to come off before this is a loglr and 2 * loglr a rho^2.
            # Without that subtraction rho^2 is max(-5e5, 0) = 0 and the
            # match allowance is identically zero at every level.
            if rho2 is None or 2.0 * (here - lognl) > rho2:
                rho2 = 2.0 * max(here - lognl, 0.0)
            allow = rho2 * numpy.sqrt(max(0.0, 1.0 - thr[d] ** 2)) / 2.0
            keep = idx[vals >= here - allow - self.loglr_region]
            logging.info('dag cut: level %i evaluated %i kept %i '
                         '(best loglr %.2f, allowance %.1f + region %.1f)',
                         d, len(idx), len(keep), here - lognl, allow,
                         self.loglr_region)
            alive = set(keep.tolist())
            perlev.append((idx.copy(), vals.copy(), keep.copy(),
                           float(here - lognl), float(allow)))
            if d == nlev - 1:
                last_loglr = numpy.full(n, -numpy.inf)
                last_loglr[idx] = vals
        if last_loglr is None:
            raise RuntimeError('dag cut never reached the finest level')
        logging.info('dag cut: %i likelihoods against %i tiles (%.1fx fewer)',
                     total, len(last_loglr),
                     len(last_loglr) / max(total, 1))
        dump = os.environ.get('GAMES6_DUMP_LEVELS')
        if dump:
            numpy.savez(dump, lognl=lognl,
                        **{'l%i_%s' % (d, k): v
                           for d, row in enumerate(perlev)
                           for k, v in zip(('idx', 'vals', 'keep', 'best',
                                            'allow'), row)})
            logging.info('wrote the per-level cut to %s', dump)
        return last_loglr, total

    def _evaluate_start_nodes(self, treefile):
        """ Likelihood at the representative of every start node.

        Returns (start_nodes, loglrs), where `start_nodes` is None if the
        cut runs over the flat map's root tiles rather than a tree depth.
        """
        start_nodes = None
        if self.dag_cull and self.dagfile:
            loglrs, nev = self._dag_cut()
            self.meta['setup_ncalls'] = nev
            return None, loglrs
        if treefile is not None and self.tree_start_level > 0:
            start_nodes = []
            for k in sorted(treefile['tree'], key=int):
                start_nodes.extend(self._nodes_at_depth(
                    treefile['tree'][k], self.tree_start_level))
            args = [dict(getattr(self, 'extra_reference', {}),
                         **{p: float(n.attrs['param_%s' % p])
                            for p in self.dtype.names})
                    for n in start_nodes]
            logging.info('evaluating the cut at tree depth %i: %i nodes '
                         '(%i at depth 0)', self.tree_start_level,
                         len(start_nodes), len(treefile['tree']))
        else:
            with h5py.File(self.mapfile, 'r') as mapfile:
                bank = {p: mapfile['bank'][p][:] for p in self.dtype.names}
            args = [dict(getattr(self, 'extra_reference', {}),
                         **{p: bank[p][i] for p in bank})
                    for i in range(len(next(iter(bank.values()))))]
        logging.info('Calculating likelihood at nodes')
        loglrs = numpy.array(list(tqdm.tqdm(
            self.pool.imap(call_tile_likelihood, args), total=len(args))))
        self.meta['setup_ncalls'] = len(args)
        return start_nodes, loglrs

    def _leaves_as_bins(self, treefile, passed, start_nodes, node_loglrs,
                        bound):
        """ Descend into each passed tile and keep the leaves that clear
        `bound`. Fills self.dmap and self.leaf_loglr.
        """
        self.dmap, self.leaf_loglr = {}, {}
        lengths, mass, key = [], [], 0
        with h5py.File(self.mapfile, 'r') as mapfile:
            for tile in passed:
                if start_nodes is not None:
                    leaves = self._find_active_leaves(
                        start_nodes[tile], bound,
                        root_loglr=node_loglrs[tile])
                elif treefile is not None and 'tree' in treefile \
                        and str(tile) in treefile['tree']:
                    leaves = self._find_active_leaves(
                        treefile['tree'][str(tile)], bound,
                        root_loglr=node_loglrs[tile])
                else:
                    leaves = []
                if not leaves and start_nodes is None:
                    # no tree entry for this tile: use the flat pool
                    pts = mapfile['map'][str(tile)][:]
                    leaves = [(pts, node_loglrs[tile], float(len(pts)))]
                for pts, ll, m in leaves:
                    self.dmap[key] = pts
                    self.leaf_loglr[key] = ll
                    lengths.append(len(pts))
                    mass.append(m)
                    key += 1
        return (numpy.arange(key), numpy.array(lengths),
                numpy.array(mass, dtype=float))

    def _active_bins(self, treefile, start_nodes, node_loglrs):
        """ The bins to draw from, as (ids, lengths, prior_mass).

        Keeps every node within `loglr_region` of the best, then resolves
        those to leaves. Two cases need something else:

        coarse   the cut cannot discriminate, because the whole spread of
                 node likelihoods is smaller than the region. Descending
                 would enumerate the tree and still drop leaves on
                 representative-level noise, so use the tiles directly --
                 their subtree clouds union to the prior.
        starved  so few tiles pass, holding so few points, that the pool
                 has nothing to give. Widen the region until it does.
                 Free, since every start node was already evaluated, and
                 unbiased, since the extra tiles are low-likelihood and
                 weighted accordingly.
        """
        finite = node_loglrs[numpy.isfinite(node_loglrs)]
        if not len(finite):
            raise RuntimeError(
                'no start node had a finite likelihood; the map does not '
                'cover this signal, or the model is misconfigured')
        bound = finite.max() - self.loglr_region
        # NaN and -inf both compare False, so they are excluded here
        passed = numpy.where(node_loglrs > bound)[0]
        # the tile set is what the dag cut has to reproduce, so it is logged
        # in full rather than just counted, and GAMES6_DUMP_CUT names a file
        # to write the per-tile likelihoods to for a value-level comparison
        logging.info('cut kept tiles: %i of %i, best %.4f, bound %.4f',
                     len(passed), len(node_loglrs), finite.max(), bound)
        logging.info('cut kept tiles: %s', ','.join(map(str, passed)))
        dump = os.environ.get('GAMES6_DUMP_CUT')
        if dump:
            numpy.save(dump, node_loglrs)
            logging.info('wrote the cut likelihoods to %s', dump)

        if start_nodes is not None \
                and finite.max() - finite.min() < self.loglr_region:
            logging.info('%i of %i start nodes passed: the cut cannot '
                         'discriminate, using tiles directly',
                         len(passed), len(node_loglrs))
            return self._tiles_as_bins(passed, start_nodes, node_loglrs)

        logging.info('...resolving tree leaves for %i passed tiles',
                     len(passed))
        ids, lengths, mass = self._leaves_as_bins(
            treefile, passed, start_nodes, node_loglrs, bound)
        logging.info('...resolved to %i active bins (%i points)',
                     len(ids), int(lengths.sum()))
        if start_nodes is None or lengths.sum() >= self.min_active_points:
            return ids, lengths, mass

        logging.warning('active set starved (%i bins, %i points): falling '
                        'back to start-level tiles', len(ids),
                        int(lengths.sum()))
        region = self.loglr_region
        for _ in range(6):
            sel = numpy.where(node_loglrs > finite.max() - region)[0]
            ids, lengths, mass = self._tiles_as_bins(sel, start_nodes,
                                                     node_loglrs)
            if lengths.sum() >= self.min_active_points:
                break
            region *= 2.0
            logging.warning('still starved (%i points); widening the '
                            'region to %.0f', int(lengths.sum()), region)
        logging.info('...fallback gave %i tiles (%i points)', len(ids),
                     int(lengths.sum()))
        return ids, lengths, mass

    def run(self):
        with h5py.File(self.mapfile, 'r') as mapfile:
            self.dtype = mapfile['map']['0'].dtype
        self._setup_extra_params()
        self._setup_map_prior()

        treefile = h5py.File(self.treefile, 'r') if self.treefile else None
        try:
            start_nodes, node_loglrs = self._evaluate_start_nodes(treefile)
            active_ids, lengths, prior_mass = self._active_bins(
                treefile, start_nodes, node_loglrs)
        finally:
            if treefile is not None:
                treefile.close()

        if self.leaf_weight_samples > 1:
            self._refine_leaf_weights(active_ids, self.leaf_weight_samples)
        weight = self._round1_weights(len(active_ids), lengths)
        # per-stratum accumulators
        strata = {'pool': {'samp': None, 'loglr': None, 'logw': None}}
        self.stratum_calls = {}
        self._logz_terms = {}
        self._deferred = set()
        self._defer_trial = None
        self._defer_eff_before = None
        self._defer_frozen = False
        self._pool_eff = None
        self._gen_eff = []
        self._switched = False
        self._switch_at = self.gen_switch_ess
        self._switch_bad = 0
        proposal = self._seed_proposal_from_leaves() \
            if self.gen_seed_leaves else None
        pool_dead = False
        self._diag_binid = None
        self._leaf_acc = {}

        # The generative proposal is held in reserve: the pool keeps the
        # whole budget each round until it exhausts and hands it over.
        for rnd in range(1, self.rounds + 1):
            # The generative proposal no longer waits for the pool to
            # exhaust: `gen_switch_ess` hands over as soon as there is
            # enough to fit a kernel with, and `_switch_bad` hands back if
            # that turns out to be premature. See `gen_switch_ess`.
            use_gen = pool_dead or self._switched
            if not use_gen:
                try:
                    ps, pl, pw, bid = self.pool_round(
                        weight / weight.sum(), active_ids, lengths,
                        self.target_likelihood_calls, prior_mass=prior_mass)
                    strata['pool'] = self._accumulate(strata['pool'],
                                                      ps, pl, pw)
                    # This round's own rate, measured exactly as a
                    # generative round's is. The CUMULATIVE pool average
                    # would flatter it, because the pool decays as it
                    # drains, and every comparison against the kernel would
                    # inherit that bias.
                    self._pool_eff = self._ess(pw) / max(len(pw), 1)
                    # _accumulate PREPENDS the new round, so match that order
                    self._diag_binid = bid if self._diag_binid is None else \
                        numpy.concatenate([bid, self._diag_binid])
                    # concentrate later rounds on the leaves that are
                    # actually carrying posterior weight
                    if self.pool_adapt_mean:
                        # Re-weight leaves by their running MEAN likelihood
                        # rather than by the posterior weight they happened
                        # to realise. softmax(logw) is top-heavy: measured,
                        # it moved 100% of the weight mass in one round with
                        # 64% of it landing in the single leaf holding the
                        # single highest-weight atom, so the update chased
                        # per-atom noise, drained that leaf, and killed the
                        # pool by round 3. The mean over a leaf's drawn
                        # atoms estimates the same quantity the weight is
                        # meant to track, at a fraction of the variance.
                        self._update_leaf_means(bid, pl)
                        upd = self._temper_leaf_weights(
                            self._leaf_mean_ll(len(active_ids)), lengths)
                        if upd is not None:
                            weight = upd
                    else:
                        wn = numpy.exp(pw - logsumexp(pw))
                        for j, v in zip(bid, wn):
                            weight[j] += v
                except OutOfSamples:
                    logging.info('pool exhausted at round %i; handing its '
                                 'budget to the generative proposal', rnd)
                    pool_dead = True
                    if proposal is None:
                        proposal = self._fit_proposal(strata)
                    if proposal is None:
                        logging.info('no generative proposal available; stop')
                        break

            if use_gen and proposal is not None:
                got = self.gen_round(proposal,
                                     self.target_likelihood_calls)
                if got is not None:
                    gs, gl, gw = got
                    # one stratum PER ROUND: each round used a different
                    # refitted proposal, so its samples are only commensurate
                    # with themselves
                    rk = 'gen:r%i' % rnd
                    strata[rk] = {'samp': gs, 'loglr': gl, 'logw': gw}
                    self.stratum_calls[rk] = len(gw)
                    # max normalised weight alongside the ESS: a round
                    # whose ESS collapses because ONE draw dominates is a
                    # different failure from one that is uniformly poor,
                    # and only the first is a proposal-tail problem.
                    wn = numpy.exp(gw - logsumexp(gw))
                    # Per-round efficiency, and the pool's, are what the
                    # deferral gate is judged on -- see `_defer_gate`.
                    self._gen_eff.append(
                        self._ess(gw) / max(len(gw), 1))
                    logging.info('round %i generative: ESS %.1f from %i '
                                 'calls (%.2f%%), max weight %.3f',
                                 rnd, self._ess(gw), len(gw),
                                 100.0 * self._ess(gw) / max(len(gw), 1),
                                 wn.max())

            # Refit as posterior information accumulates. Safe because
            # each sample already carries the weight computed against the
            # proposal in force when it was drawn.
            if rnd >= self.gen_start_round and self.gen_refit > 0 and \
                    (rnd - self.gen_start_round) % self.gen_refit == 0:
                new = self._fit_proposal(strata)
                if new is not None:
                    proposal = new

            ess_p = self._ess(strata['pool']['logw'])
            gkeys = [q for q in strata if q.startswith('gen')]
            gsum = sum(self._ess(strata[q]['logw']) for q in gkeys)
            self._switch_check(ess_p + gsum, proposal, pool_dead)
            logging.info('round %i: ESS pool=%.1f gen[%i rounds]=%.1f '
                         'total=%.1f ncalls=%i', rnd, ess_p, len(gkeys),
                         gsum, ess_p + gsum, self.ncalls)

        self._finalise(strata)

    def _refine_leaf_weights(self, active_ids, k):
        """ Re-estimate each active leaf's representative loglr as the MEAN
        likelihood over `k` of its own atoms.

        A leaf's weight should be proportional to the likelihood INTEGRATED
        over it, i.e. E[LR] * volume, but `leaf_loglr` holds the likelihood
        AT one representative. Those differ by however much loglr varies
        inside the leaf, and the extrinsic marginalization makes that
        variation larger than the tiling's match threshold bounds: measured
        on a hard injection, one representative sits 2.8 nats (std) from its
        own leaf's mean, so exp() of it misallocates the round-1 draws by
        factors of e^3 between leaves that deserve equal weight.

        Averaging k atoms cuts that noise as sqrt(k) and costs k calls per
        ACTIVE leaf only (~100s of calls), not per cut node (~40k). This
        changes only the draw probability, never the estimator: the
        importance weight carries log(drawcount) and self-corrects, so this
        is a pure variance choice.
        """
        args, owners = [], []
        for key in active_ids:
            pts = self.dmap[int(key)]
            n = len(pts)
            if not n:
                continue
            take = min(int(k), n)
            sel = self._rng.choice(n, size=take, replace=False)
            for i in sel:
                args.append(dict(self.extra_reference,
                                 **{p: float(pts[p][i])
                                    for p in self.dtype.names}))
                owners.append(int(key))
        if not args:
            return
        lls = list(self.pool.imap(call_tile_likelihood, args))
        self.tree_ncalls += len(args)
        got = {}
        for key, ll in zip(owners, lls):
            if numpy.isfinite(ll):
                got.setdefault(key, []).append(ll)
        moved = []
        for key, v in got.items():
            new = logsumexp(v) - numpy.log(len(v))
            old = self.leaf_loglr.get(key, numpy.nan)
            if numpy.isfinite(old):
                moved.append(new - old)
            self.leaf_loglr[key] = new
        logging.info('refined %i leaf weights with %i calls (%i atoms each): '
                     'mean loglr moved %.2f nats (std %.2f) from the single '
                     'representative', len(got), len(args), int(k),
                     float(numpy.mean(moved)) if moved else 0.0,
                     float(numpy.std(moved)) if moved else 0.0)

    def _round1_weights(self, key, lengths):
        """ Draw probability for round 1, over `key` active bins.

        Leaf size alone is prior volume, so round 1 would spread draws by
        volume and start concentrating them only from round 2. Seeding
        with exp(loglr_rep) * volume reuses the representative
        likelihoods the cut already paid for. This changes the variance
        of the importance weights, not their expectation: logw carries
        log(drawcount), so they self-correct for any draw probability.

        The exponent is tempered because exp(loglr_rep) is the likelihood
        AT the representative, while a leaf's weight is the likelihood
        integrated OVER it. The two diverge as the within-leaf likelihood
        spread grows, and an untempered exponent piles every draw into a
        handful of leaves and drains them at once. Tempering holds the
        effective dynamic range near `gen_seed_range` nats.
        """
        volume = lengths / lengths.sum()
        if not self.gen_seed_leaves:
            return volume
        ll = numpy.array([self.leaf_loglr.get(k, numpy.nan)
                          for k in range(key)], dtype=float)
        w = self._temper_leaf_weights(ll, lengths)
        if w is None:
            return volume
        return w

    def _temper_leaf_weights(self, ll, lengths):
        """ Draw probability from per-leaf loglr estimates `ll`.

        exp(loglr) * volume, with the exponent tempered to hold the
        effective dynamic range near `gen_seed_range` nats. Returns None if
        no leaf has a usable estimate. Shared by the round-1 seed and the
        between-round update so both weight leaves on the same scale.
        """
        ok = numpy.isfinite(ll)
        if not ok.any():
            return None
        spread = float(ll[ok].max() - ll[ok].min())
        temper = max(1.0, spread / max(self.gen_seed_range, 1e-6))
        w = lengths.astype(float).copy()
        w[ok] *= numpy.exp((ll[ok] - ll[ok].max()) / temper)
        w[~ok] = 0.0 if ok.all() else w[~ok]
        if w.sum() <= 0:
            return None
        logging.info('leaf weights: temper %.2f over a loglr spread of %.1f '
                     'nats in %i leaves', temper, spread, int(ok.sum()))
        return w / w.sum()

    def _update_leaf_means(self, bid, loglr):
        """ Accumulate a running mean likelihood per leaf from drawn atoms.

        The leaf's weight should track the likelihood INTEGRATED over it,
        which the atoms drawn from it estimate directly:
        logsumexp(loglr) - log(n) is an unbiased estimate of log E[LR].
        Accumulated across rounds, so a leaf's estimate keeps sharpening
        with every atom it has ever given up.
        """
        for j in numpy.unique(bid):
            sel = loglr[bid == j]
            sel = sel[numpy.isfinite(sel)]
            if not len(sel):
                continue
            ls, n = self._leaf_acc.get(int(j), (-numpy.inf, 0))
            # pool_adapt_mean == 2 accumulates the mean of loglr (geometric
            # mean likelihood) instead of the log of the mean likelihood.
            # log E[LR] is the quantity the weight should track, but its
            # estimator logsumexp(loglr) - log(n) is dominated by the single
            # largest atom, so it keeps most of the per-atom variance we are
            # trying to average out: measured, it GREW the between-leaf
            # spread from 15 to 36 nats. The arithmetic mean of loglr is
            # biased low but every atom contributes equally, and for a DRAW
            # probability the bias is free -- the importance weight carries
            # log(drawcount) and self-corrects.
            if self.pool_adapt_mean == 2:
                add = float(numpy.sum(sel))
                self._leaf_acc[int(j)] = ((add if n == 0 else ls + add),
                                          n + len(sel))
            else:
                self._leaf_acc[int(j)] = (numpy.logaddexp(ls, logsumexp(sel)),
                                          n + len(sel))

    def _leaf_mean_ll(self, nleaves):
        """ Per-leaf log-mean-likelihood, falling back to the leaf's own
        representative where no atom has been drawn yet. """
        out = numpy.full(nleaves, numpy.nan)
        for k in range(nleaves):
            if k in self._leaf_acc and self._leaf_acc[k][1] > 0:
                ls, n = self._leaf_acc[k]
                # mode 2 stores sum(loglr); mode 1 stores logsumexp(LR)
                out[k] = ls / n if self.pool_adapt_mean == 2 \
                    else ls - numpy.log(n)
            else:
                out[k] = self.leaf_loglr.get(k, numpy.nan)
        return out

    def _setup_extra_params(self):
        """ Parameters the model samples but the map does not carry.

        A map only needs the parameters that change the waveform enough to
        matter for its match-based tiling. Others -- in-plane spin
        components for a BNS, say, where precession is weak -- can be left
        out of the map and drawn from the prior instead, which keeps the
        map from growing (its node count goes as mismatch^(-d_eff/2)).

        This is exact, not an approximation. The pool proposal becomes
        q = q_map(map params) * prior(extra), while the target is
        pi = L * prior(map) * prior(extra), so the prior(extra) factors
        cancel and the weight is unchanged. Support is complete because
        the draw covers the whole prior.

        Where those parameters are genuinely unconstrained the prior draw
        is close to optimal. Where they are not, the pool stratum loses
        efficiency and the generative proposal -- which is fitted over
        ALL sampled parameters -- learns them and takes over, which is
        the handover that already exists.
        """
        names = list(self.model.variable_params)
        self.extra_params = [p for p in names if p not in self.dtype.names]
        self.samp_dtype = numpy.dtype([(p, numpy.float64) for p in names])
        if not self.extra_params:
            return
        # A fixed reference for the cut. The cut compares nodes by the
        # likelihood at their representative, and the map says nothing
        # about these parameters, so giving each node a different random
        # value would add noise to a comparison that is already the
        # weakest link. The median of the prior is the neutral choice.
        draws = self.model.prior_distribution.rvs(size=20000)
        self.extra_reference = {p: float(numpy.median(draws[p]))
                                for p in self.extra_params}
        logging.info('map does not carry %s; drawing them from the prior '
                     '(cut evaluated at their prior median)',
                     ', '.join(self.extra_params))
        # The cancellation above needs q_extra = prior(extra | map), and a
        # draw from the joint prior gives the MARGINAL. Those agree only
        # when the missing parameters are prior-independent of the mapped
        # ones. A constraint spanning both -- |chi| <= a_max over all three
        # spin components, with z mapped and x,y not -- breaks that, and
        # would bias the result rather than merely cost efficiency.
        con = getattr(self.model.prior_distribution, '_constraints', None)
        if con:
            logging.warning(
                'the prior carries %i constraint(s). If any couples %s to '
                'the mapped parameters, the pool weights are WRONG: the '
                'missing parameters are drawn from their marginal, which '
                'equals the conditional only when they are independent of '
                'the map. Use independent priors for them, or put them in '
                'the map.', len(con), ', '.join(self.extra_params))

    def _setup_map_prior(self):
        """ The prior the MAP's points were drawn from.

        The pool's points are prior draws, so they represent the map's
        prior, not necessarily the model's. Where the two differ each
        point's weight picks up pi_model/pi_map -- which is what lets one
        map serve a different prior, and what makes the evidence
        normalisation meaningful.

        Silently inert for maps built before the builder recorded its
        configuration; those are assumed to share the model's prior, which
        is what the sampler assumed unconditionally until now.
        """
        self.map_prior = None
        self.map_size = 0
        with h5py.File(self.mapfile, 'r') as f:
            cfg = f.attrs.get('config')
            self.map_size = int(f.attrs.get('sample_size', 0))
            if not self.map_size:
                # Recover it from the map itself. The builder draws exactly
                # `sample_size` prior points and stores all of them, so the
                # total stored count IS the normalisation the pool evidence
                # needs. Only the dataset shapes are read, not the points.
                # Without this the pool stratum is dropped for every map
                # built before the builder recorded its metadata, and the
                # evidence then rests on the generative stratum alone with
                # no cross-check.
                self.map_size = int(sum(f['map'][k].shape[0]
                                        for k in f['map']))
                logging.info('map carries no sample_size; derived %i prior '
                             'draws from its contents', self.map_size)
        if self.map_size:
            logging.info('map holds %i prior draws', self.map_size)
        if cfg is None:
            logging.warning('map carries no build configuration: assuming '
                            'its prior matches the model\'s, and the '
                            'evidence cannot be normalised')
            return
        import tempfile
        from pycbc.workflow import WorkflowConfigParser
        from pycbc.distributions import read_distributions_from_config
        from pycbc.distributions.joint import JointDistribution
        with tempfile.NamedTemporaryFile('w', suffix='.ini',
                                         delete=False) as fh:
            fh.write(cfg if isinstance(cfg, str) else cfg.decode())
            path = fh.name
        try:
            cp = WorkflowConfigParser([path])
            dists = read_distributions_from_config(cp, 'prior')
            names = [p_ for d in dists for p_ in d.params]
            self.map_prior = JointDistribution(names, *dists)
            self.map_prior_params = [p_ for p_ in names
                                     if p_ in self.dtype.names]
            logging.info('map prior read for %s',
                         ', '.join(self.map_prior_params))
        except Exception as exc:
            logging.warning('could not rebuild the map prior (%s); assuming '
                            'it matches the model', exc)
            self.map_prior = None

    def _log_prior_ratio(self, psamp):
        """ log pi_model - log pi_map over the MAPPED parameters.

        Zero when the two priors agree, which is the case the sampler has
        always assumed. Non-zero lets one map serve a differently shaped
        prior over the same parameterisation.
        """
        if self.map_prior is None:
            return 0.0
        names = self.map_prior_params
        # The model prior is a joint over ALL sampled params, so calling
        # it with only the mapped subset raises "must provide all
        # parameters". We need
        # its MARGINAL over the mapped params, which -- because the priors
        # factorise -- is the sum of the independent sub-distributions confined
        # to those params. Identify them once.
        nameset = set(names)
        sub = [d for d in self.model.prior_distribution.distributions
               if set(d.params) <= nameset]
        covered = (set().union(*[set(d.params) for d in sub])
                   if sub else set())
        if covered != nameset:
            # a mapped param shares a sub-distribution with an unmapped one, so
            # its marginal is not separable here; fall back to the assumption
            # the sampler has always made (priors agree over mapped params).
            logging.warning('evidence: model prior over mapped params '
                            'not separable (%s uncovered); assuming it '
                            'matches the map prior',
                            ', '.join(sorted(nameset - covered)))
            return 0.0
        n = len(psamp)
        lm = numpy.empty(n)
        lp = numpy.empty(n)
        for i in range(n):
            kw = {p_: float(psamp[p_][i]) for p_ in names}
            lm[i] = self.map_prior(**kw)
            lp[i] = sum(float(d(**{q: kw[q] for q in d.params}))
                        for d in sub)
        r = lp - lm
        finite = r[numpy.isfinite(r)]
        if len(finite) and numpy.allclose(finite, finite[0]):
            c = float(finite[0])
            if abs(c) < 1e-12:
                return 0.0                  # priors identical: a true no-op
            # A non-zero CONSTANT ratio (e.g. a narrow map prior against a broad
            # model prior) cancels in the self-normalised posterior, so it does
            # not touch parameters or ESS -- but the pool evidence is NOT
            # self-normalised, so dropping it offsets log Z by exactly this
            # constant (the source of the "strata disagree" warning). Keep it,
            # preserving non-finite entries so out-of-model-prior points stay
            # zero-weight.
            return numpy.where(numpy.isfinite(r), c, r)
        return r

    def _draw_extra(self, psamp):
        """ Fill the not-in-map parameters with prior draws, in place. """
        if not self.extra_params:
            return
        d = self.model.prior_distribution.rvs(size=len(psamp))
        for p in self.extra_params:
            psamp[p] = d[p]

    def _as_samples(self, pts):
        """ Map-parameter records -> a FieldArray over ALL sampled
        parameters, with the missing ones drawn from the prior.
        """
        out = FieldArray(len(pts), dtype=self.samp_dtype)
        for p in self.dtype.names:
            out[p] = pts[p]
        self._draw_extra(out)
        return out

    def _accumulate(self, st, samp, loglr, logw):
        if st['samp'] is None:
            return {'samp': samp, 'loglr': loglr, 'logw': logw}
        return {'samp': numpy.concatenate([samp, st['samp']]),
                'loglr': numpy.concatenate([loglr, st['loglr']]),
                'logw': numpy.concatenate([logw, st['logw']])}

    def _seed_proposal_from_leaves(self, max_points=200000):
        """ Fit a generative proposal to the active leaves, before any
        posterior samples exist.

        A leaf's posterior mass is approximately exp(loglr_rep) times its
        prior volume, and its point count is proportional to that volume,
        so weighting every point in a leaf by exp(loglr_rep) turns the
        union of active leaves into an estimate of the posterior. The
        representative likelihoods were already paid for by the cut.

        Built as a synthetic stratum and handed to `_fit_proposal`, so it
        goes through the same centre subsampling and kernel construction
        as every other fit.
        """
        names = list(self.model.variable_params)
        keys = [k for k, v in self.leaf_loglr.items()
                if numpy.isfinite(v) and len(self.dmap[k])]
        if not keys:
            return None
        total = sum(len(self.dmap[k]) for k in keys)
        frac = min(1.0, max_points / max(total, 1))
        xs, lw = [], []
        for k in keys:
            pts = self.dmap[k]
            n = len(pts)
            if frac < 1.0:
                take = max(1, int(n * frac))
                sel = self._rng.choice(n, size=take, replace=False)
                pts = pts[sel]
            xs.append(pts)
            lw.append(numpy.full(len(pts), float(self.leaf_loglr[k])))
        samp = numpy.concatenate(xs)
        logw = numpy.concatenate(lw)
        # Temper the seed weights so the KDE receives a non-degenerate
        # spread of centres. Untempered, exp(loglr_rep) collapses onto the
        # single best leaf in the high-SNR limit -- the fit then subsamples
        # centres by weight without replacement and gets almost no support,
        # so a sharp posterior with a thin seed cannot bootstrap. Holding
        # the effective range near gen_seed_range nats (the same recipe as
        # _round1_weights) balloons the seed shape just enough to continue;
        # later rounds refit against real posterior samples, which carry
        # their own weights and let the proposal shrink back to shape.
        spread = float(logw.max() - logw.min())
        temper = max(1.0, spread / max(self.gen_seed_range, 1e-6))
        if temper > 1.0:
            logging.info('tempered seed proposal: temper %.2f over a loglr '
                         'spread of %.1f nats', temper, spread)
        logw = logw / temper
        fa = self._as_samples(samp)
        stratum = {'samp': fa, 'loglr': logw, 'logw': logw}
        prop = self._fit_proposal({'seed': stratum})
        if prop is not None:
            logging.info('seeded generative proposal from %i active leaves '
                         '(%i points, subsampled to %i)', len(keys), total,
                         len(samp))
        return prop

    def _fit_proposal(self, strata):
        """ Fit the generative proposal to the accumulated posterior samples
        from every stratum, resampled to equal weight.
        """
        names = list(self.model.variable_params)
        # Weight each stratum's contribution by its own ESS. Normalising
        # every stratum to sum to 1 and concatenating would give each
        # equal total mass regardless of quality, letting a weak stratum
        # drag the fit toward wherever it happened to draw. This is the
        # same beta-proportional-to-ESS choice `_finalise` uses.
        xs, ws = [], []
        for st in strata.values():
            if st['samp'] is None:
                continue
            wv = numpy.exp(st['logw'] - logsumexp(st['logw']))
            ess = self._ess(st['logw'])
            if ess <= 0:
                continue
            xs.append(numpy.column_stack([st['samp'][p] for p in names]))
            ws.append(wv * ess)
        if not xs:
            return None
        x = numpy.vstack(xs)
        w = numpy.concatenate(ws)
        if w.sum() <= 0:
            return None
        w = w / w.sum()
        ess = 1.0 / (w ** 2).sum()
        if ess < self.gen_min_ess:
            logging.debug('accumulated ESS %.1f below gen_min_ess %.1f; '
                          'not fitting yet', ess, self.gen_min_ess)
            return None
        # Centres are subsampled by weight WITHOUT replacement when there
        # are more than gen_max_centres. That keeps most of the pool's
        # weight while never picking the same point twice, which matters
        # because duplicated centres make the kernel covariance
        # degenerate at low ESS.
        # Sampling without replacement needs at least `n` entries with
        # non-zero probability. In the exhausted regime most weights underflow
        # to exactly zero -- at SNR 100 the loglr range spans thousands of
        # nats -- so the usable count can fall below gen_max_centres and
        # numpy raises "Fewer non-zero entries in p than size". Whether it
        # does is seed-dependent, which is how this survived until a
        # precessing SNR 100 run happened to land on the wrong side.
        nz = int(numpy.count_nonzero(w))
        n = min(self.gen_max_centres, len(x), nz)
        if n < 2:
            logging.warning('only %i accumulated draws carry non-zero weight, '
                            'too few to fit a proposal; the pool is exhausted '
                            'at this SNR', nz)
            return None
        if nz < min(self.gen_max_centres, len(x)):
            logging.info('only %i of %i accumulated draws carry non-zero '
                         'weight, so the fit uses %i centres rather than %i '
                         '-- a sign of exhaustion, not an error', nz, len(x),
                         n, min(self.gen_max_centres, len(x)))
        if n < len(x):
            # w is already normalised to sum to 1 above
            order = self._rng.choice(len(x), size=n, replace=False, p=w)
        else:
            order = numpy.arange(len(x))
        centres = x[order]
        cw = w[order]
        if cw.sum() <= 0:
            return None
        # The EFFECTIVE number of mixture components. Nominal centre count
        # says little: the component weights are the posterior weights, and
        # when those are peaked a nominally large mixture behaves like a
        # handful of Gaussians. This, not gen_max_centres, is what limits
        # how closely the proposal can approach the posterior.
        logging.info('centres=%i effective components=%.1f (of %i drawn '
                     'from %i)', len(cw), cw.sum() ** 2 / (cw ** 2).sum(),
                     n, len(x))
        defer = self._defer_dimensions(x, w, names) if self.gen_defer else []
        model_idx = [i for i in range(len(names)) if i not in defer]
        if defer:
            logging.info('deferring %s to the prior this round; modelling %s',
                         ', '.join(names[i] for i in defer),
                         ', '.join(names[i] for i in model_idx))
        try:
            box = self._box() if self.gen_truncate else None
            if box is not None and defer:
                box = (box[0][model_idx], box[1][model_idx])
            prop = LocalCovarianceProposal(
                centres[:, model_idx] if defer else centres,
                self._rng, k=self.gen_local_k,
                scale=self.gen_bandwidth, weights=cw,
                bounds=box)
            if defer:
                di = list(defer)

                def _draw(size, _di=di, _names=names):
                    d = self.model.prior_distribution.rvs(size=size)
                    return numpy.column_stack([d[_names[i]] for i in _di])

                dlp = self._defer_logpdf([names[i] for i in di])
                prop = FactorizedProposal(prop, model_idx, di, _draw,
                                          len(names), defer_logpdf=dlp)
        except numpy.linalg.LinAlgError:
            logging.info('generative proposal fit failed (singular kernel)')
            return None
        logging.info('fitted proposal: k=%i bw=%.2f centres=%i '
                     'accumulated ESS=%.1f', self.gen_local_k,
                     self.gen_bandwidth, prop.n, ess)
        return prop

    def _defer_logpdf(self, defer_names):
        """ log prior density of the deferred parameters alone, or None.

        The deferred factor cancels in the posterior but not in the
        evidence, so the proposal has to report it. The model's prior is a
        JointDistribution of independent pieces, so the deferred density is
        the sum of the pieces whose parameters are entirely deferred.

        Returns None if any deferred parameter is not covered that way --
        a piece spanning both deferred and modelled parameters cannot be
        factorised, which is the same prior-independence condition the
        deferral itself needs. The caller then gets the old behaviour, so
        the posterior is unaffected and only the evidence is missing a
        constant, which is warned about.
        """
        want = list(defer_names)
        pieces = getattr(self.model.prior_distribution, 'distributions', None)
        if not pieces:
            logging.warning('deferral: the prior does not expose independent '
                            'pieces, so the deferred prior density cannot be '
                            'factorised; the posterior is unaffected but '
                            'log Z will be high by that density')
            return None
        use, covered = [], set()
        for d in pieces:
            ps = list(d.params)
            if ps and all(q in want for q in ps):
                use.append(d)
                covered.update(ps)
        if covered != set(want):
            logging.warning('deferral: %s not covered by independent prior '
                            'pieces, so the deferred prior density cannot be '
                            'factorised; the posterior is unaffected but '
                            'log Z will be high by that density',
                            ', '.join(sorted(set(want) - covered)))
            return None

        def _lp(xd, _use=use, _want=want):
            out = numpy.zeros(len(xd))
            for i in range(len(xd)):
                kw = {q: float(xd[i, _want.index(q)]) for q in _want}
                out[i] = sum(float(d(**{q: kw[q] for q in d.params}))
                             for d in _use)
            return out
        return _lp

    def _defer_dimensions(self, x, w, names):
        """ Parameters the likelihood has said nothing about, by DIRECTION.

        A parameter can be handed to the prior exactly -- a uniform density
        over the prior range has that range as its support, so there is no
        boundary deficit and no contained-mass estimate for it at all. That
        removes the kernel's shape error in those coordinates rather than
        patching it, and measured on the ten-parameter case it is worth 2.5x
        on top of truncation (relative cost 2.53/2.82/2.59 -> 1.00/1.06/1.00
        over three seeds).

        The test has to be on DIRECTIONS, not marginals. In prior units the
        prior covariance is the identity, so the eigenvalues of the posterior
        covariance say how much each direction was constrained, and a
        parameter is deferrable only if it loads on none of the constrained
        ones. Testing each marginal against its prior instead is unsound and
        actively dangerous: `spin2z` has a prior-like marginal purely because
        the chi_eff ridge spans the box in that projection, while being one of
        the best measured directions in the problem. A KS marginal test
        defers it, destroys the ridge, and costs 22x (efficiency 22.4% ->
        0.99%).

        This method only MEASURES; `_defer_decide` chooses. The split may
        change freely between rounds because each round is its own stratum,
        but it is only ever changed on measured efficiency -- the old code
        changed it on this loading alone and latched, see `_defer_decide`.
        """
        box = self._box()
        if box is None or self.gen_defer <= 0:
            return []
        lo, hi = box
        psd = (hi - lo) / numpy.sqrt(12.0)
        u = (x - 0.5 * (lo + hi)) / psd
        try:
            lam, V = numpy.linalg.eigh(numpy.cov(u.T, aweights=w))
        except (numpy.linalg.LinAlgError, ValueError):
            return []
        informative = lam < 0.8
        if not informative.any():
            return []
        load = numpy.abs(V[:, informative]).max(axis=1)
        logging.info('information per direction: %s',
                     ' '.join('%.3f' % v for v in numpy.sort(lam)))
        logging.info('loading on informative directions: %s (eligible below '
                     '%.2f)', ' '.join('%s=%.2f' % (p, load[i])
                                       for i, p in enumerate(names)),
                     self.gen_defer)
        return self._defer_decide(load, names)

    def _defer_decide(self, load, names):
        """ Which eligible parameters to actually defer, decided on MEASURED
        efficiency rather than on the loading.

        The loading is a suggestion, not a trigger. Applying it speculatively
        is what broke: measured at idx 223, the fifth eigenvalue sits at
        0.757 while the informative cut is 0.8, the round-1 estimate has
        bootstrap scatter +-0.031 AND is biased high because the seed
        proposal has already deferred, and the resulting state LATCHES --
        deferring inflates the very eigenvalue that holds it deferred, so the
        old 1.75x release could never fire. Two seeds of one identical
        configuration gave 36.0% and 6.4%. Across the ten worst rows of a
        300-injection P-P the same latch cost 1.5-9.7x, every one of which
        recovered to 28-40% with deferral simply switched off.

        So: model everything unless the kernel has DEMONSTRATED it is sick.

        The reference is the pool stratum, not an absolute threshold, because
        it is measured on the same problem and so self-calibrates across
        SNR and dimension. Un-deferred, the kernel beat the pool by
        1.52-4.28x on all ten of those rows, while the ten-parameter
        precessing case deferral was built for runs ~2.5x WORSE than its
        deferred self. A 2x-worse-than-pool gate is a factor ~3 clear of
        both populations, which is why one rule serves both ends.

        Parameters are then added ONE AT A TIME, lowest loading first, and an
        addition is kept only if it pays. A wrong guess costs one round.
        """
        # Judge the outcome of the previous trial addition before anything
        # else: a trial that did not pay is reverted, and no more are tried.
        if self._defer_trial is not None:
            got = self._gen_eff[-1] if self._gen_eff else 0.0
            base = self._defer_eff_before or 0.0
            if got > base * self.gen_defer_margin:
                logging.info('deferring %s paid off (%.2f%% -> %.2f%%); '
                             'keeping it', names[self._defer_trial],
                             100 * base, 100 * got)
                self._defer_trial = None
            else:
                logging.info('deferring %s did not pay (%.2f%% -> %.2f%%); '
                             'reverting and stopping',
                             names[self._defer_trial], 100 * base, 100 * got)
                self._deferred.discard(self._defer_trial)
                self._defer_trial = None
                self._defer_frozen = True
                return sorted(self._deferred)

        if self._defer_frozen:
            return sorted(self._deferred)

        # The gate. Deferral is a remedy; do not treat a healthy kernel.
        if not self._defer_gate():
            return sorted(self._deferred)

        cand = [i for i in range(len(names))
                if load[i] < self.gen_defer and i not in self._deferred]
        # Never defer so much that the modelled block cannot support a
        # covariance: at low SNR one parameter can carry the single
        # informative direction on its own, and deferring the rest makes the
        # kernel fit raise, which costs the round.
        if not cand or len(names) - len(self._deferred) - 1 < 2:
            return sorted(self._deferred)
        pick = min(cand, key=lambda i: load[i])
        self._deferred.add(pick)
        self._defer_trial = pick
        self._defer_eff_before = self._gen_eff[-1] if self._gen_eff else 0.0
        logging.info('kernel is %.2f%% against pool %.2f%%; trying deferral '
                     'of %s (loading %.2f)', 100 * self._defer_eff_before,
                     100 * (self._pool_eff or 0.0), names[pick], load[pick])
        return sorted(self._deferred)

    def _switch_check(self, ess, proposal, pool_dead):
        """ Hand the budget to the kernel early, and hand it back if wrong.

        The pool was previously kept until it exhausted. Measured at
        handover, the kernel beat the pool at EVERY seed size -- 2.38x when
        the seed was under 2000 ESS, 1.39x above 20000 -- and converged to
        68-79% regardless of it. So waiting is not free, and 37% of a
        300-injection population never built a kernel at all.

        Nothing here is sticky. A switch that turns out premature hands back
        and is retried once the accumulated ESS has grown by
        `gen_switch_backoff`, so the number of attempts is logarithmic in the
        budget. That is deliberately the opposite of the `gen_defer` latch,
        where a single early decision was inherited by every later round.
        """
        if self.gen_switch_ess <= 0 or pool_dead:
            return
        if not self._switched:
            if ess >= self._switch_at and proposal is not None:
                self._switched = True
                self._switch_bad = 0
                logging.info('switching to the generative proposal at ESS '
                             '%.0f (pool was running at %.2f%%); the pool is '
                             'not exhausted', ess,
                             100.0 * (self._pool_eff or 0.0))
            return
        # switched: is it actually beating the pool it displaced?
        if self._pool_eff is None or not self._gen_eff:
            return
        if self._gen_eff[-1] < self._pool_eff:
            self._switch_bad += 1
            if self._switch_bad >= self.gen_switch_patience:
                self._switched = False
                self._switch_bad = 0
                self._switch_at = ess * self.gen_switch_backoff
                logging.info('the kernel ran below the pool (%.2f%% vs '
                             '%.2f%%) for %i rounds; handing back, will retry '
                             'at ESS %.0f', 100.0 * self._gen_eff[-1],
                             100.0 * self._pool_eff, self.gen_switch_patience,
                             self._switch_at)
        else:
            self._switch_bad = 0

    def _defer_gate(self):
        """ True when the kernel has measurably earned a remedy.

        Poor relative to the pool for `gen_defer_patience` consecutive rounds
        AND not improving over that window. The second condition matters: the
        kernel is legitimately poor while it has few centres and climbs on
        its own -- on a healthy run it went 50%% to 78%% over 13 rounds -- so
        a dip that is recovering must not trigger a remedy.
        """
        if self._pool_eff is None or self.gen_defer_ratio <= 0:
            return False
        recent = self._gen_eff[-self.gen_defer_patience:]
        if len(recent) < self.gen_defer_patience:
            return False
        if not all(e < self._pool_eff / self.gen_defer_ratio for e in recent):
            return False
        return recent[-1] <= recent[0]


    def _box(self):
        """ Per-parameter prior bounds, or None if any is unbounded. """
        names = list(self.model.variable_params)
        lo = numpy.full(len(names), -numpy.inf)
        hi = numpy.full(len(names), numpy.inf)
        try:
            dists = self.model.prior_distribution.distributions
        except AttributeError:
            return None
        for d in dists:
            for k, v in (getattr(d, 'bounds', None) or {}).items():
                if k in names:
                    i = names.index(k)
                    lo[i], hi[i] = float(v[0]), float(v[1])
        if not numpy.all(numpy.isfinite(lo) & numpy.isfinite(hi)):
            return None
        return lo, hi

    def _finalise(self, strata):
        """ Combine strata with beta proportional to ESS, which is the
        inverse-variance choice and makes the combined ESS additive.
        """
        parts = [(k, v) for k, v in strata.items() if v['samp'] is not None]
        esss = numpy.array([self._ess(v['logw']) for _, v in parts])
        total_ess = float(esss.sum())
        beta = esss / esss.sum() if esss.sum() > 0 else \
            numpy.ones(len(parts)) / len(parts)
        for (k, _), e, b in zip(parts, esss, beta):
            c = self.stratum_calls.get(k, 0)
            logging.info('stratum %s: ESS=%.1f calls=%i eff=%.2f%% beta=%.3f',
                         k, e, c, 100.0 * e / c if c else 0.0, b)
            self.meta[f'ess_{k}'] = float(e)
            self.meta[f'ncalls_{k}'] = int(c)
        logging.info('combined ESS=%.1f ncalls=%i (setup %i, tree cut %i)',
                     total_ess, self.ncalls,
                     self.meta.get('setup_ncalls', 0), self.tree_ncalls)
        # The extrinsic marginalization's own effective sample size, summed
        # over the run. If this is small the importance weights carry the
        # sky integral's Monte-Carlo noise and no intrinsic proposal can
        # recover the loss, so it distinguishes "the proposal is bad" from
        # "the likelihood we handed the proposal is noisy".
        n = getattr(self.model, 'vess_n', 0)
        if n:
            mean = self.model.vess_sum / n
            self.meta['vector_ess_mean'] = float(mean)
            self.meta['vector_ess_min'] = float(self.model.vess_min)
            self.meta['vector_ess_max'] = float(self.model.vess_max)
            self.meta['vector_ess_calls'] = int(n)
            logging.info('extrinsic marginalization ESS: mean %.1f, range '
                         '%.1f-%.1f of %i samples over %i calls (%.2f%%)',
                         mean, self.model.vess_min, self.model.vess_max,
                         getattr(self.model, 'vsamples', 0), n,
                         100.0 * mean / max(getattr(self.model, 'vsamples', 1),
                                            1))

        # Setup and tree-cut calls are deliberately NOT folded into ncalls:
        # both are one-time costs paid before/while resolving which leaves
        # are active, independent of how many posterior samples are then
        # drawn, so they amortize away at the ESS this sampler targets.
        # Reported separately because they are still a useful diagnostic --
        # a tree whose root level is too coarse prunes weakly and descends
        # a large fraction of its nodes, which shows up here.
        self._finalise_evidence(parts, beta)
        self.meta['ncalls'] = self.ncalls
        self.meta['ess'] = total_ess
        self.meta['setup_ncalls'] = self.meta.get('setup_ncalls', 0)
        self.meta['tree_ncalls'] = self.tree_ncalls

        # The file stores the WEIGHTED draws, following the convention
        # dynesty set: raw samples plus a `logwt` dataset and a
        # `log_evidence` attribute, with the resampling to equal weight done
        # on READ by the io class. Downstream tools all go through
        # `read_samples`, so nothing else has to know about weights.
        #
        # It matters because resampling discards information that cannot be
        # recovered: the effective sample size, the evidence and its error,
        # and the Pareto tail index of the weights all need the individual
        # weights. Those diagnostics used to be computable only inside a
        # running process.
        #
        # Strata are put on a common scale first. Each carries its own
        # self-normalised weights, and beta -- proportional to its ESS -- is
        # what says how much of the posterior it speaks for, so the stored
        # weight is beta_s * w_i / sum_s(w).
        outs, ols, ows = [], [], []
        for (k, v), b in zip(parts, beta):
            lw = numpy.asarray(v['logw'], float)
            if not len(lw) or b <= 0:
                continue
            outs.append(v['samp'])
            ols.append(numpy.asarray(v['loglr'], float))
            ows.append(lw - logsumexp(lw) + numpy.log(b))
        if not outs:
            # Every stratum contributed nothing: total ESS collapsed to ~0,
            # so no stratum's beta rounded up to even one output sample.
            # Reachable with a very permissive cut -- e.g. loglr_region=60
            # on a coarse start level admits nearly the whole prior, and
            # round-1 draws (which are uniform over the active region) then
            # almost never land anywhere with meaningful likelihood. Fail
            # loudly with the diagnosis rather than dying in concatenate.
            raise RuntimeError(
                'no posterior samples: combined ESS=%.3g across %i strata. '
                'The active region is almost certainly far too large -- '
                'reduce loglr_region or use a finer tree_start_level.'
                % (total_ess, len(parts)))
        fsamp = numpy.concatenate(outs)
        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = numpy.concatenate(ols)
        self._samples['logwt'] = numpy.concatenate(ows)
        logging.info('stored %i weighted draws (ESS %.0f); the io class '
                     'resamples them on read', len(fsamp), total_ess)

    @property
    def io(self):
        from pycbc.inference.io.games import GamesFile
        return GamesFile

    def finalize(self):
        """ Write the draws, then the evidence the io class needs to weight
        them. `log_evidence` is where every reader looks for it.
        """
        super().finalize()
        # Add any parameters the model drew inline. These are part of
        # default_stats, so they are read back from the model's cache;
        # uncached points are recomputed.
        if getattr(self.model, 'reconstruct_stats', None):
            self.write_cached_stats()
        if 'logz' in self.meta:
            with self.io(self.checkpoint_file, 'a') as fp:
                # the TOTAL, not the within-proposal floor: this is the
                # dlog_evidence every reader will treat as the uncertainty
                fp.write_logevidence(self.meta['logz'],
                                     self.meta.get('logz_sigma_total',
                                                   self.meta.get(
                                                       'logz_sigma',
                                                       float('nan'))))

    def _evidence_stats(self, parts):
        """ Per stratum: the within-proposal error, the Pareto tail index of
        the weights, the draw count and the ESS.

        Keyed by the AGGREGATE stratum name to match `_logz_terms`. The
        per-round parts are named `gen:r7` while the evidence accumulates under
        `gen`, so keying on part names silently missed every generative stratum
        and left its diagnostics NaN.
        """
        # Statistical error, per stratum, from the weights themselves.
        #
        # Z is the unnormalised MEAN of the weights, so Var(Z)/Z^2 = CV^2/N,
        # and since ESS = N/(1+CV^2) exactly, the relative error collapses to
        # sqrt(1/ESS - 1/N). Nothing extra has to be stored to get it.
        #
        # It is a WITHIN-PROPOSAL error and understates the total. Measured
        # at aligned 13.5: six repeats at one seed spread 0.010 nats (sd
        # 0.004, matching this formula), six at different seeds spread 0.115
        # (sd 0.044) -- 11x wider. The extra term is the proposal's own
        # randomness across seeds, which centres are subsampled, which draws
        # the contained-mass estimates used, where the pool landed. This
        # formula cannot see it, because each run's weights are computed
        # against the proposal that run actually built. Quote the seed-to-seed
        # spread when it matters.
        #
        # It is only meaningful if the weights have finite variance, so the
        # Pareto tail index is reported next to it. A large k means the
        # sample variance understates the true one and sigma is a lower
        # bound, which is the honest thing to say rather than quoting a
        # precise-looking number.
        # Keyed by the AGGREGATE stratum name, matching `_logz_terms`. The
        # per-round parts are named `gen:r7` while the evidence is
        # accumulated under `gen`, so keying on the part names silently
        # missed every generative stratum and left its diagnostics as NaN.
        pooled = {}
        for name, v in parts:
            pooled.setdefault(name.split(':')[0], []).append(
                numpy.asarray(v['logw'], float))
        stats = {}
        for name, lws in pooled.items():
            lw = numpy.concatenate(lws)
            n = len(lw)
            ess = self._ess(lw)
            var = max(1.0 / ess - 1.0 / n, 0.0) if ess > 0 and n else 0.0
            stats[name] = (numpy.sqrt(var), _pareto_k(lw), n, ess)
        return stats

    def _proposal_scale(self, parts):
        """ Rough proposal-scale error from one run, early rounds against late.
        """
        # Proposal-scale error, from ONE run: early rounds against late.
        #
        # sqrt(1/ESS - 1/N) is a WITHIN-proposal error, and the run-to-run
        # spread is larger: at aligned 13.5 six repeats at one seed spread
        # 0.010 nats while six at different seeds spread 0.115.
        #
        # Most of that excess turned out NOT to be the proposal. The model's
        # sky pool is drawn once per construction and offsets every likelihood
        # call in the run identically, contributing 0.046 on its own at this
        # SNR against a 0.043 observed logZ spread (MODEL_ISSUES.md 2a). This
        # term measures what the proposal adds beyond that, so expect it to be
        # small, and do not read a small value as the estimator failing.
        #
        # A run already contains many proposals, since every refit is a new one
        # and every round is its own stratum with its own unbiased estimate of
        # the SAME Z. But estimating the spread round BY round does not work,
        # and was measured returning exactly 0.000 against a true 0.043:
        # adjacent refits share almost all their accumulated samples, overlap
        # ~0.9 by round nine, so the estimates are correlated, the observed
        # scatter is deflated by (1 - rho), and recovery needs dividing by
        # ~0.1. Splitting the rounds in half instead leaves ~0.5 overlap, so
        # the same correction is ~2x, and the halves use genuinely different
        # proposals since the early ones are fitted on a fraction of the data.
        splits = []
        for name, terms in self._logz_terms.items():
            per = [t[0] - numpy.log(float(t[1])) for t in terms if t[1]]
            rounds = [(nm, v) for nm, v in parts if nm.split(':')[0] == name]
            if len(per) < 4 or len(rounds) != len(per):
                continue
            s2 = []
            for (_, v), _z in zip(rounds, per):
                lw = numpy.asarray(v['logw'], float)
                n, e = len(lw), self._ess(lw)
                s2.append(max(1.0 / e - 1.0 / n, 1e-12) if e > 0 and n else
                          1e-12)
            z = numpy.array(per)
            w = 1.0 / numpy.array(s2)

            # Each half is unbiased, so a difference beyond the two
            # within-half errors is proposal variability.
            half = len(z) // 2
            if half >= 2:
                def _blk(sl):
                    ww = w[sl]
                    return ((ww * z[sl]).sum() / ww.sum(), 1.0 / ww.sum())
                za, va = _blk(slice(0, half))
                zb, vb = _blk(slice(half, None))
                excess = (za - zb) ** 2 - (va + vb)
                self.meta['logz_split_%s' % name] = float(
                    numpy.sqrt(max(excess, 0.0)))
                splits.append(numpy.sqrt(max(excess, 0.0)))
                logging.info('evidence: %s early rounds %.3f vs late %.3f, '
                             'difference %.4f against %.4f expected from '
                             'their own errors -> proposal term %.4f',
                             name, za, zb, abs(za - zb),
                             numpy.sqrt(va + vb),
                             numpy.sqrt(max(excess, 0.0)))
        return splits

    def _descent_uncertainty(self, logz, sigma_l=0.15):
        """ What the likelihood cut cost the evidence, from the cut itself.

        The dominant evidence uncertainty is not the proposal, it is which
        region the descent declared active: each node is judged on one
        stochastic likelihood evaluation against a fixed threshold, so nodes
        near it are decided by noise, and a flip includes or excludes a whole
        subtree. Measured across six seeds, both strata shift together by
        0.043 nats while the within-proposal error is 0.004.

        No repeats are needed to quantify it, because the descent already knows
        everything required. A dropped node carries `n_points`, its subtree's
        share of the map, and its representative's loglr is already paid for,
        so its forfeited share of Z = E_prior[LR] is
        (n_i / M_total) exp(loglr_i):

            bias      sum over dropped nodes -- the cut only ever loses mass
            variance  the same terms weighted by p_i (1 - p_i), with
                      p_i = Phi((loglr_i - bound) / sigma_L) the chance the
                      decision would have gone the other way

        `sigma_l` is the likelihood's own stochastic scatter, ~0.15 nats at
        SNR 13.5 and ~0.68 at SNR 100 (see MODEL_ISSUES.md); it is the one
        number this needs from outside. A representative also understates its
        subtree, so this is a floor rather than a bound -- tightening it needs
        the bank's match and the SNR, which belong to the map.
        """
        if not self._dropped or not self.map_size:
            return
        from scipy.stats import norm
        ll = numpy.array([d[0] for d in self._dropped])
        n = numpy.array([d[1] for d in self._dropped])
        # each dropped subtree's forfeited contribution to Z, in logs
        # relative to Z from the start: loglr carries a large offset, so
        # exp(2*terms) overflows to inf and inf*0 is nan when a decision is
        # certain (p = 0 or 1).  Everything stays in logs.
        terms = ll + numpy.log(n) - numpy.log(float(self.map_size)) - logz
        # against the cut the descent used, not the region's width:
        # the docstring above specifies (loglr_i - bound).
        cut = getattr(self, '_last_bound', None)
        if cut is None:
            cut = float(numpy.max(ll)) - self.loglr_region
        p = norm.cdf((ll - cut) / max(sigma_l, 1e-6))
        # a Bernoulli per subtree: kept or not, weighted by what it carries
        with numpy.errstate(divide='ignore'):
            logvar = 2.0 * terms + numpy.log(p) + numpy.log1p(-p)
        self.meta['logz_cut_bias'] = float(numpy.exp(logsumexp(terms)))
        self.meta['logz_cut_sigma'] = float(
            numpy.exp(0.5 * logsumexp(logvar)) if numpy.isfinite(logvar).any()
            else 0.0)
        logging.info('evidence: the cut dropped %i subtrees holding a '
                     'fraction %.2e of Z (a floor on the bias, since a '
                     'representative understates its subtree), and marginal '
                     'decisions among them contribute %.2e relative -- '
                     'compare the quoted statistical error',
                     len(ll), self.meta['logz_cut_bias'],
                     self.meta['logz_cut_sigma'])

    def _finalise_evidence(self, parts, beta):
        """ log Z from the unnormalised weights each stratum accumulated.

        Importance sampling gives the evidence for free: Z = E_q[L pi / q],
        so the unnormalised MEAN of the weights is an estimate of it, and
        self-normalising -- which every stratum does for the posterior --
        is exactly what throws it away.

        The pool's denominator is the map's total prior draw count, since
        summing its weights over every point in the map would give
        M_total * Z. The generative rounds use the number of draws they
        requested.

        Strata are combined with the same beta the posterior uses. That is
        the inverse-variance choice for the POSTERIOR, not necessarily for
        the evidence, so treat the spread across strata as the honest
        uncertainty rather than the quoted one.
        """
        stats = self._evidence_stats(parts)

        zs = {}
        for name, terms in self._logz_terms.items():
            num = [t[0] for t in terms]
            den = [t[1] for t in terms]
            if not num or not all(den):
                # Say so rather than dropping it quietly. The pool's
                # denominator is the map's total prior draw count, which maps
                # built before the builder recorded its metadata do not carry
                # -- and then only the generative stratum contributes while
                # `logz_spread` reads 0.0, which looks like two strata
                # agreeing perfectly instead of one stratum being alone.
                logging.warning('evidence: dropping the %s stratum, its '
                                'normalisation is unavailable (for the pool '
                                'this means the map has no sample_size; '
                                'rebuild it with pycbc_inference_build_map '
                                'to get the cross-check)', name)
                continue
            # mean over rounds of (sum w / N)
            per = numpy.array(num) - numpy.log(numpy.array(den, dtype=float))
            zs[name] = logsumexp(per) - numpy.log(len(per))
        if not zs:
            logging.info('evidence not available (map size or prior missing)')
            return
        for k, v in zs.items():
            if k in stats:
                sig, pk, n, ess = stats[k]
                logging.info('log Z from the %s stratum: %.3f +- %.3f '
                             '(ESS %.0f of %i draws, Pareto k %.2f%s)',
                             k, v, sig, ess, int(n), pk,
                             '' if not pk > 0.7 else ' -- HEAVY TAIL, sigma '
                             'is a lower bound')
            else:
                sig = pk = float('nan')
                logging.info('log Z from the %s stratum: %.3f', k, v)
            self.meta['logz_%s' % k] = float(v)
            self.meta['logz_sigma_%s' % k] = float(sig)
            self.meta['logz_paretok_%s' % k] = float(pk)
        splits = self._proposal_scale(parts)

        # Combine the strata by INVERSE VARIANCE, and drop any whose weights
        # have no finite mean.
        #
        # Equal weights were the first cut and are indefensible once a stratum
        # degrades: at aligned SNR 100 the pool is exhausted (ESS 51, Pareto k
        # 1.88, so no finite mean at all) yet carried the same weight as a
        # healthy generative stratum, dragging log Z by 0.19 nats -- larger
        # than the model's own floor. The posterior has always used ESS
        # weighting; this brings the evidence into line.
        keys = list(zs)
        wts, use = [], []
        for k in keys:
            sig = stats[k][0] if k in stats else float('nan')
            pk = stats[k][1] if k in stats else float('nan')
            if numpy.isfinite(sig) and sig > 0 and not pk > 1.0:
                use.append(k)
                wts.append(1.0 / sig ** 2)
        if not use:
            # every stratum is unusable; fall back to equal weights over all
            # of them rather than reporting nothing, and say so
            logging.warning('evidence: no stratum has both a finite error and '
                            'a finite mean, so log Z is an equal-weight '
                            'average of unreliable estimates -- treat it as '
                            'indicative only')
            use, wts = keys, [1.0] * len(keys)
        elif len(use) < len(keys):
            logging.warning('evidence: dropping the %s stratum from log Z, '
                            'its weights have no finite mean (Pareto k > 1)',
                            ', '.join(k for k in keys if k not in use))
        vals = numpy.array([zs[k] for k in use])
        w = numpy.array(wts) / numpy.sum(wts)
        logz = float(numpy.sum(w * vals))
        self._descent_uncertainty(logz)
        # Combined statistical error: the strata are averaged with equal
        # weight above, so their variances add as 1/n^2 sum sigma_i^2.
        # Inverse-variance weights, so the combined error is 1/sqrt(sum 1/s^2).
        #
        # That formula assumes the strata are INDEPENDENT, and they are not
        # quite: the generative proposal is fitted to accumulated pool samples,
        # so the two share information and the combined error here is
        # optimistic. The same assumption was buried in the previous
        # equal-weight version; it is more visible now because inverse-variance
        # combining shrinks the number further. Treat the cross-stratum
        # disagreement (`logz_spread`) as the honest check, as the warning
        # below says.
        comb = (float(1.0 / numpy.sqrt(numpy.sum(wts)))
                if len(use) > 1 or numpy.isfinite(wts[0]) else float('nan'))
        # Report a TOTAL that carries the proposal scale, not just the
        # within-proposal floor.
        #
        # Everything here is the SAMPLER's own error. The run-to-run scatter is
        # larger and mostly not ours: the model's marginalization pools are
        # drawn once per construction, giving a frozen per-run offset measured
        # at 0.064 nats (0.049 - 0.091) at aligned 13.5 against a 0.008
        # statistical error (MODEL_ISSUES.md 2a). That floor cannot be
        # estimated from inside one run, which sees a single frozen draw, so it
        # is warned about rather than folded in -- inventing a number for it
        # would be worse than naming it.
        #
        # The early-versus-late term is rough and one degree of freedom and can
        # come back zero when unresolved; the total then says so rather than
        # passing the floor off as the total.
        prop = max(splits) if splits else 0.0
        # The cut's marginal decisions are a VARIANCE and belong in the total.
        # Its missing mass is a one-sided BIAS -- the cut can only lose mass --
        # so it is reported separately and NOT added in quadrature; Z is low by
        # that fraction rather than uncertain by it.
        cut = self.meta.get('logz_cut_sigma', 0.0)
        total = float(numpy.sqrt(comb ** 2 + prop ** 2 + cut ** 2))
        self.meta['logz'] = float(logz)
        self.meta['logz_sigma'] = comb
        self.meta['logz_sigma_total'] = total
        self.meta['logz_sigma_proposal'] = float(prop)
        self.meta['logz_nstrata'] = int(len(vals))
        # Spread over EVERY stratum, not just the ones combined. Dropping a
        # stratum for k > 1 must not destroy the cross-check -- that is
        # precisely the case where the disagreement is the informative number,
        # and reporting nan there loses the only honest diagnostic.
        allz = numpy.array(list(zs.values()))
        self.meta['logz_spread'] = float(allz.max() - allz.min()) \
            if len(allz) > 1 else float('nan')
        logging.info('evidence: the quoted error is the SAMPLER\'s. The '
                     'model adds a per-run offset of its own, drawn once with '
                     'its marginalization pools, measured at 0.064 nats at '
                     'aligned SNR 13.5 and growing with SNR; one run cannot '
                     'estimate it, so it is NOT included here (MODEL_ISSUES '
                     '2a)')
        if len(vals) > 1:
            logging.info('log Z = %.3f +- %.3f  [within-proposal %.4f, '
                         'proposal %.4f, cut %.4f%s]; spread across %i strata '
                         '%.3f', logz, total, comb, prop, cut,
                         ' UNRESOLVED, treat the total as a floor'
                         if prop <= 0 else '', len(vals),
                         self.meta['logz_spread'])
            if self.meta['logz_spread'] > 3 * comb:
                logging.warning('evidence: the strata disagree by %.3f, well '
                                'beyond the %.3f statistical error -- trust '
                                'the spread, not the error bar',
                                self.meta['logz_spread'], comb)
        else:
            logging.info('log Z = %.3f +- %.3f  [within-proposal %.4f, '
                         'proposal %.4f, cut %.4f%s] from the %s stratum '
                         'ALONE; no cross-check available', logz, total, comb,
                         prop, cut,
                         ' UNRESOLVED, treat the total as a floor'
                         if prop <= 0 else '', list(zs)[0])

    def _subtree_points(self, group, cap=None):
        """ Points under `group`, subsampled to `cap`, and the true total.

        The counterpart to `_find_active_leaves` for the case where the
        cut cannot prune: that one pays a likelihood call per node to
        decide what to keep, this one keeps everything below a node that
        already passed, so its cost is HDF5 reads rather than likelihood
        calls.

        The cap matters because a start-level tile's subtree can be the
        whole map, which a run has no use for. Returns (points, total) so
        the caller can keep prior mass and physical size apart: `total`
        is the tile's share of the prior, `len(points)` is only how many
        are available to draw.
        """
        chunks, total = [], 0

        def gather(g):
            nonlocal total
            if g.attrs['is_leaf']:
                n = int(g.attrs['n_points'])
                total += n
                if n:
                    pts = numpy.zeros(n, dtype=self.dtype)
                    # only the MAPPED params are stored in the leaf; extras
                    # (e.g. in-plane spins for a BNS) are drawn from the prior
                    # downstream by _as_samples, exactly as the descend() path
                    # does. Iterating self.variable_params here instead read
                    # points_<extra> and KeyError'd on broad (low-SNR)
                    # posteriors that fall into this gather path.
                    for p in self.dtype.names:
                        pts[p] = g[f'points_{p}'][:]
                    chunks.append(pts)
                return
            for c in g['children']:
                gather(g['children'][c])

        gather(group)
        if not chunks:
            return numpy.zeros(0, dtype=self.dtype), 0
        pts = numpy.concatenate(chunks)
        if cap is not None and len(pts) > cap:
            # uniform subsample: the stored cloud stays prior-distributed
            # inside the tile, which is what the pool draw assumes
            sel = numpy.random.choice(len(pts), size=cap, replace=False)
            pts = pts[sel]
        return pts, total

    def _tiles_as_bins(self, tiles, start_nodes, node_loglrs):
        """ Use whole start-level tiles as the pool, no likelihood calls.

        Shared by the two coarse paths. Each tile becomes one bin holding
        its subtree's points (capped) and its already-evaluated rep
        loglr; the bin's PRIOR MASS is the true subtree population, which
        is why the cap does not bias it.
        """
        self.dmap, self.leaf_loglr = {}, {}
        lengths, mass, key = [], [], 0
        for tile in tiles:
            pts, total = self._subtree_points(start_nodes[tile],
                                              self.tile_point_cap)
            if not len(pts):
                continue
            self.dmap[key] = pts
            self.leaf_loglr[key] = node_loglrs[tile]
            lengths.append(len(pts))
            mass.append(total)
            key += 1
        return (numpy.arange(key), numpy.array(lengths),
                numpy.array(mass, dtype=float))

    def _find_active_leaves(self, tree_root, loglr_bound, root_loglr=numpy.nan):
        """ Prune the subtree, evaluating the likelihood at each child's
        representative and descending only where it clears the bound.

        These are real likelihood calls and are counted into
        `self.tree_ncalls`, which is reported separately from `ncalls`.

        `root_loglr` is the loglr already measured at `tree_root` by the
        cut. It seeds the leaf representative when `tree_root` is itself a
        leaf (a tree exactly `tree_start_level` deep), which would
        otherwise store nan and silently disable the loglr-weighted
        round-1 seed and the leaf-seeded generative proposal.
        """
        results = []

        def collect(group):
            """ Every point and prior mass under a node, with NO likelihood
            calls. Used when the cut has stopped discriminating. """
            if group.attrs['is_leaf']:
                pts = numpy.zeros(int(group.attrs['n_points']),
                                  dtype=self.dtype)
                for p in self.dtype.names:
                    pts[p] = group[f'points_{p}'][:]
                return [pts], [float(group.attrs.get('prior_mass',
                                                     len(pts)))]
            pp, mm = [], []
            for c in group['children']:
                a, b = collect(group['children'][c])
                pp += a
                mm += b
            return pp, mm

        def descend(group, ll_here):
            if group.attrs['is_leaf']:
                pts = numpy.zeros(int(group.attrs['n_points']),
                                  dtype=self.dtype)
                for p in self.dtype.names:
                    pts[p] = group[f'points_{p}'][:]
                # Third element is the leaf's PRIOR MASS. A refined leaf
                # holds more points than its mass warrants, so the weight
                # must use the builder's unbiased count or oversampling
                # would inflate that leaf; equals the point count for
                # unrefined maps. The loglr is the leaf's own
                # representative, already paid for during the descent and
                # what makes the leaf-weighted cloud an estimate of the
                # posterior (see _seed_proposal_from_leaves).
                results.append((pts, ll_here,
                                float(group.attrs.get('prior_mass',
                                                      len(pts)))))
                return
            # A node's children are evaluated as ONE batch through the
            # pool. They used to go one at a time, and that is what made the
            # descent expensive: measured over 60 P-P injections, the three
            # slowest runs were exactly the three that descended, paying ~20
            # extra minutes for only 2059 calls against 70000 sampled ones --
            # about 0.6 s per descent call against 0.009 s per sampled one,
            # because serial calls get none of the pool's parallelism. The
            # call COUNT was never the problem.
            kids, args = [], []
            for c in group['children']:
                child = group['children'][c]
                kids.append(child)
                args.append(dict(self.extra_reference,
                                 **{q: float(child.attrs[f'param_{q}'])
                                    for q in self.dtype.names}))
            if not kids:
                return
            lls = (list(self.pool.imap(call_tile_likelihood, args))
                   if len(args) > 1 else [call_tile_likelihood(args[0])])
            self.tree_ncalls += len(args)
            ok = [numpy.isfinite(v) and v > loglr_bound for v in lls]

            # What the cut throws away, and how surely.
            #
            # This is where the evidence's real uncertainty comes from. Every
            # node is judged on ONE stochastic likelihood evaluation against a
            # fixed threshold, so nodes near it are decided by noise: measured,
            # tree_ncalls varies 3.3% between seeds and both strata then shift
            # in lockstep by 0.043 nats, an order of magnitude beyond the
            # within-proposal error.
            #
            # It does not need repeating runs to quantify. A dropped node
            # carries `n_points`, its own subtree's share of the map, and its
            # representative's loglr was just paid for, so its forfeited
            # contribution to Z = E_prior[LR] is (n_i / M_total) exp(loglr_i).
            # Summed over dropped nodes that is the BIAS -- the cut can only
            # lose mass, never invent it. Weighted by how nearly each decision
            # went the other way it is the VARIANCE, the same quantities
            # giving both.
            #
            # The representative understates its subtree, which is where the
            # bank's match enters: within a tile loglr spans about
            # (1 - m^2) rho^2 / 2, so exp(loglr_i + that) bounds what the
            # subtree could have held. Recorded raw here and combined at the
            # end, since rho and m belong to the map rather than the descent.
            for child, ll, good in zip(kids, lls, ok):
                if good or not numpy.isfinite(ll):
                    continue
                n = float(child.attrs.get('n_points', 0))
                if n > 0:
                    self._dropped.append((ll, n))

            # Stop descending where the cut has stopped discriminating.
            #
            # The descent only earns its likelihood calls by pruning. When
            # almost every child of a node passes, refining that node buys
            # nothing and simply pays for the whole subtree: below SNR ~10
            # nothing is pruned anywhere and the descent enumerates every
            # node, 150k-450k calls against ~4k for a sampling round. Taking
            # the passing children whole costs no calls at all.
            #
            # It is also the regime where the pool needs no help: a broad
            # posterior is close to the prior, and the pool is prior
            # distributed, so accepting unpruned points loses little.
            # Deliberately per node rather than per level, since the passing
            # fraction varies across a level.
            if (self.tree_accept_fraction > 0 and len(ok) > 1
                    and numpy.mean(ok) >= self.tree_accept_fraction):
                for child, ll, good in zip(kids, lls, ok):
                    if not good:
                        continue
                    pp, mm = collect(child)
                    if pp:
                        results.append((numpy.concatenate(pp), ll,
                                        float(numpy.sum(mm))))
                self.tree_shortcut += 1
                return

            for child, ll, good in zip(kids, lls, ok):
                if good:
                    descend(child, ll)

        descend(tree_root, root_loglr)
        if self.tree_shortcut:
            logging.info('descent: accepted %i subtrees whole where the cut '
                         'passed at least %.0f%% of children',
                         self.tree_shortcut, 100 * self.tree_accept_fraction)
        return results
