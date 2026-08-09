""" games3 plus inexhaustible generative proposals, combined by strata.

This keeps games3 exactly as-is (hierarchical tree of discrete pools, which
reaches ~32% efficiency at low SNR) and adds one or more *generative*
proposals that can supply unlimited draws once the discrete pools start to run
dry. The point is narrow: fix exhaustion without giving up games3's efficiency.

Two generative proposals, both fitted to the *accumulated posterior samples*
rather than to a tile's prior pool:

  'kde'      a Gaussian mixture centred on the current posterior samples, with
             a full sample-covariance bandwidth.
  'stretch'  additional mixture centres placed along the segment between pairs
             of posterior samples, x = a + u (b - a). This fills in between
             samples along whatever direction the posterior actually
             correlates, which is the same idea as an affine-invariant stretch
             move. A bare line proposal has no Lebesgue density in 6-D (a
             segment is measure zero), so the segment points are used as
             mixture *centres* and given the same kernel as 'kde'; the density
             is then exactly computable.

This is deliberately NOT games2. games2 fit a KDE to each tile's own prior
pool, i.e. to a thin curved sliver of the prior, and a full-rank Gaussian
cannot represent that -- confirmed there and independently here (a leaf's
covariance spectrum spans ~13 orders of magnitude and the cloud is curved).
The posterior, by contrast, is a compact correlated blob, which is the shape a
Gaussian mixture handles well.

Combination by strata, not by a single mixture density
------------------------------------------------------
The discrete-pool proposal puts mass on atoms while the generative proposals
have Lebesgue densities; the two measures are mutually singular, so there is no
common density in which to write one mixture. Rather than paper over that (a
good way to introduce a silent bias), each proposal is kept as its own stratum
with its own self-normalised estimate, and the strata are combined as
    I = sum_s beta_s I_s,      sum_s beta_s = 1,
which is valid for any beta. Taking beta_s proportional to the stratum's ESS is
the inverse-variance choice and makes the combined ESS the sum of the strata's,
so a stratum that is doing badly simply contributes little. That is what
implements "use whichever proposal is working" without any tuning.

Consequences worth stating:
  * Each sample carries the weight of the proposal and round that produced it.
    Weights are never recomputed against a later proposal, which would bias.
  * Out-of-prior draws get zero weight and are dropped BEFORE the likelihood
    call, so they cost draws but not likelihood evaluations, and dropping them
    does not bias a self-normalised estimate.
  * Pool exhaustion no longer ends the run: the exhausted stratum's budget is
    handed to the generative stratum for the remaining rounds.

Validated configuration
------------------------
The constructor defaults below ARE the validated settings (see
example_bns/GAMES6_RESULTS.md): local-covariance kernel (``gen_local_k=160``,
not the global ``GaussianMixtureProposal`` fallback), raw coordinates
(``gen_logit=0`` -- the logit/prior-box map made things worse, kept only as a
documented negative result), no defensive component (``gen_defensive=0``),
reserve mode (``gen_fraction=0`` -- the pool keeps the full budget until it
exhausts, then hands over), and weighted-without-replacement centre
subsampling once accumulated samples exceed ``gen_max_centres``. Only two
knobs need a per-case override away from these defaults: ``gen_max_centres``
and ``gen_min_ess`` at the most extreme SNR (see prior_g6final_snr100.ini).
"""
import os
import logging
import tqdm
import h5py
import numpy
import numpy.random
from scipy.special import logsumexp
from pycbc.io import FieldArray

from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler
from .games import call_likelihood, OutOfSamples


class GaussianMixtureProposal:
    """ Gaussian mixture with a shared full-covariance kernel.

    Centres come from the accumulated posterior samples and, optionally, from
    stretch points between pairs of them. The density is a plain mixture, so it
    is exact and cheap: one logsumexp over centres per query point.
    """

    def __init__(self, centres, cov, rng, stretch_pairs=0,
                 stretch_scale=1.0, weights=None, defensive_frac=0.0,
                 defensive_scale=3.0):
        self.rng = rng
        centres = numpy.atleast_2d(centres)
        if weights is None:
            weights = numpy.full(len(centres), 1.0 / len(centres))
        weights = numpy.asarray(weights, dtype=float)
        if stretch_pairs > 0 and len(centres) > 2:
            a = centres[rng.integers(0, len(centres), stretch_pairs)]
            b = centres[rng.integers(0, len(centres), stretch_pairs)]
            u = rng.uniform(0.0, stretch_scale, size=(stretch_pairs, 1))
            # stretch centres inherit the mean component weight, so adding
            # them reshapes the proposal without reweighting the originals
            extra = numpy.full(stretch_pairs, weights.mean())
            centres = numpy.vstack([centres, a + u * (b - a)])
            weights = numpy.concatenate([weights, extra])
        self.weights = weights / weights.sum()
        self.logw = numpy.log(numpy.maximum(self.weights, 1e-300))
        self.centres = centres
        self.n, self.dim = centres.shape

        # Shared kernel: sample covariance times a Scott-like factor. A full
        # covariance matters here -- the posterior is strongly correlated
        # (mchirp/q, and spin1z/spin2z are tightly anti-correlated), and a
        # diagonal kernel would spray draws across those directions.
        factor = self.n ** (-1.0 / (self.dim + 4))
        self.cov = cov * factor ** 2
        # jitter until positive definite; posterior samples can be degenerate
        # in a direction the sampler has barely explored
        eps = 1e-12 * numpy.trace(self.cov) / self.dim
        for _ in range(8):
            try:
                self.chol = numpy.linalg.cholesky(self.cov)
                break
            except numpy.linalg.LinAlgError:
                self.cov = self.cov + eps * numpy.eye(self.dim)
                eps *= 10
        else:
            raise numpy.linalg.LinAlgError('kernel covariance not PD')
        self.log_det = 2.0 * numpy.log(numpy.diag(self.chol)).sum()
        self.inv = numpy.linalg.inv(self.cov)
        self._lognorm = -0.5 * (self.dim * numpy.log(2 * numpy.pi)
                                + self.log_det)

        # Defensive component: a single broad Gaussian covering the whole
        # sample spread, mixed in with probability defensive_frac. Efficiency
        # here is not limited by mild over-dispersion -- a proposal 10% too
        # wide in 6-D would still be ~91% efficient -- but by rare draws where
        # q is small while the target is not, which produce enormous weights
        # and collapse the ESS. Mixing in a broad component bounds the weight
        # ratio pi/q from above everywhere the broad component has support,
        # which is the standard defensive-mixture remedy. It costs a little
        # efficiency in the bulk and removes the tail catastrophe.
        self.defensive_frac = float(defensive_frac)
        if self.defensive_frac > 0:
            self.def_mean = numpy.average(self.centres, axis=0,
                                          weights=self.weights)
            spread = numpy.cov(self.centres.T, aweights=self.weights)
            self.def_cov = spread * defensive_scale ** 2
            eps2 = 1e-12 * numpy.trace(self.def_cov) / self.dim
            for _ in range(8):
                try:
                    self.def_chol = numpy.linalg.cholesky(self.def_cov)
                    break
                except numpy.linalg.LinAlgError:
                    self.def_cov = self.def_cov + eps2 * numpy.eye(self.dim)
                    eps2 *= 10
            else:
                self.defensive_frac = 0.0
        if self.defensive_frac > 0:
            self.def_inv = numpy.linalg.inv(self.def_cov)
            self.def_lognorm = -0.5 * (
                self.dim * numpy.log(2 * numpy.pi)
                + 2.0 * numpy.log(numpy.diag(self.def_chol)).sum())

    def sample(self, size):
        idx = self.rng.choice(self.n, size=size, p=self.weights)
        z = self.rng.normal(size=(size, self.dim))
        x = self.centres[idx] + z @ self.chol.T
        if self.defensive_frac > 0:
            use = self.rng.random(size) < self.defensive_frac
            k = int(use.sum())
            if k:
                zb = self.rng.normal(size=(k, self.dim))
                x[use] = self.def_mean + zb @ self.def_chol.T
        return x

    def logpdf(self, x):
        """ log density, exact: logsumexp over centres minus log(n). """
        out = numpy.empty(len(x))
        # chunked to bound memory at n_centres * chunk * dim
        chunk = max(1, int(4e7 / max(self.n, 1)))
        for s in range(0, len(x), chunk):
            d = x[s:s+chunk, None, :] - self.centres[None, :, :]
            m = numpy.einsum('ijk,kl,ijl->ij', d, self.inv, d)
            out[s:s+chunk] = (logsumexp(-0.5 * m + self.logw[None, :], axis=1)
                              + self._lognorm)
        if self.defensive_frac > 0:
            zb = x - self.def_mean
            mb = numpy.einsum('ij,jk,ik->i', zb, self.def_inv, zb)
            logb = -0.5 * mb + self.def_lognorm
            out = numpy.logaddexp(
                numpy.log1p(-self.defensive_frac) + out,
                numpy.log(self.defensive_frac) + logb)
        return out


class LocalCovarianceProposal:
    """ Mixture with a SEPARATE covariance per centre, from its k nearest
    neighbours.

    A single global covariance cannot hug a thin *curved* manifold: measured at
    SNR 30, the bandwidth broad enough to cover the target puts 41% of draws
    beyond the posterior's own 99th-percentile nearest-neighbour distance, while
    the bandwidth that keeps draws on-manifold under-disperses and leaves the
    tails empty. Estimating each kernel's covariance from its own neighbourhood
    lets the kernel follow the sheet's local orientation and thickness, which is
    what a global covariance averages away.
    """

    def __init__(self, centres, rng, k=40, scale=1.0, weights=None,
                 ridge=1e-3, defensive_frac=0.0, defensive_scale=3.0):
        from scipy.spatial import cKDTree
        self.rng = rng
        self.centres = numpy.atleast_2d(centres)
        self.n, self.dim = self.centres.shape
        if weights is None:
            weights = numpy.full(self.n, 1.0 / self.n)
        self.weights = numpy.asarray(weights, float)
        self.weights /= self.weights.sum()
        self.logw = numpy.log(numpy.maximum(self.weights, 1e-300))
        self.defensive_frac = float(defensive_frac)

        # Everything is built in globally WHITENED coordinates, where every
        # parameter is O(1). Doing this in physical units is a trap: the six
        # parameters span ~6 orders of magnitude (mchirp variance ~1e-8 against
        # lambda ~1e5), so any isotropic regularisation floor applied there is
        # meaningless -- it swamps the tight parameters completely.
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
            # floor is dimensionless here, so it regularises without distorting
            c = c + ridge * f ** 2 * eye
            try:
                L = numpy.linalg.cholesky(c)
            except numpy.linalg.LinAlgError:
                L = numpy.linalg.cholesky(c + 10 * ridge * f ** 2 * eye)
            self.chol[i] = L
            self.inv[i] = numpy.linalg.inv(L @ L.T)
            self.lognorm[i] = -0.5 * (self.dim * numpy.log(2 * numpy.pi)
                                      + 2.0 * numpy.log(numpy.diag(L)).sum())

        # Defensive component: a single broad Gaussian in the same whitened
        # z-space, covering the full weighted spread of centres, mixed in with
        # probability defensive_frac. Local kernels are tight by design (each
        # is sized to its own k-NN neighbourhood), so a centre far from where
        # a query point lands can leave the local mixture density essentially
        # zero there even though the true posterior is not -- the same
        # tail-catastrophe risk the global-covariance proposal's defensive
        # component guards against, here guarding a per-centre kernel instead
        # of a single global one.
        if self.defensive_frac > 0:
            self.def_mean = numpy.average(cz, axis=0, weights=self.weights)
            spread = numpy.cov(cz.T, aweights=self.weights)
            def_cov = spread * defensive_scale ** 2
            eps2 = 1e-12 * numpy.trace(def_cov) / self.dim
            for _ in range(8):
                try:
                    self.def_chol = numpy.linalg.cholesky(def_cov)
                    break
                except numpy.linalg.LinAlgError:
                    def_cov = def_cov + eps2 * numpy.eye(self.dim)
                    eps2 *= 10
            else:
                self.defensive_frac = 0.0
        if self.defensive_frac > 0:
            self.def_inv = numpy.linalg.inv(def_cov)
            self.def_lognorm = -0.5 * (
                self.dim * numpy.log(2 * numpy.pi)
                + 2.0 * numpy.log(numpy.diag(self.def_chol)).sum())

    def sample(self, size):
        idx = self.rng.choice(self.n, size=size, p=self.weights)
        z0 = self.rng.normal(size=(size, self.dim))
        z = self._cz[idx] + numpy.einsum('nij,nj->ni', self.chol[idx], z0)
        if self.defensive_frac > 0:
            use = self.rng.random(size) < self.defensive_frac
            k = int(use.sum())
            if k:
                zb = self.rng.normal(size=(k, self.dim))
                z[use] = self.def_mean + zb @ self.def_chol.T
        return z @ self._Winv

    def logpdf(self, x):
        z = x @ self._W
        out = numpy.full(len(z), -numpy.inf)
        for i in range(self.n):
            d = z - self._cz[i]
            m = numpy.einsum('ij,jk,ik->i', d, self.inv[i], d)
            out = numpy.logaddexp(out, -0.5 * m + self.lognorm[i]
                                  + self.logw[i])
        if self.defensive_frac > 0:
            zb = z - self.def_mean
            mb = numpy.einsum('ij,jk,ik->i', zb, self.def_inv, zb)
            logb = -0.5 * mb + self.def_lognorm
            out = numpy.logaddexp(
                numpy.log1p(-self.defensive_frac) + out,
                numpy.log(self.defensive_frac) + logb)
        # change of variables from whitened z back to physical x
        return out + self._logdetW


class BoundedProposal:
    """ Wraps a mixture proposal in a logit map onto the prior box.

    The posterior here runs hard up against its prior boundaries -- it is nearly
    uniform in lambda2 across the full [0, 1000] range, and 20% of its mass sits
    within 2% of some edge -- and a Gaussian kernel cannot represent a sharp
    edge. Fitting in the raw parameters leaks 26% of draws outside the prior and,
    worse, leaves q deficient just inside each edge (measured: weights within 2%
    of an edge are 3x the interior), because the kernel mass that should sit
    there fell outside instead. Uneven q means weight variance, which is ESS.

    Mapping each parameter to y = logit((x - a) / (b - a)) removes both effects
    at once: the support becomes exactly the prior box, so nothing leaks and
    nothing is rejected, and a density that is flat against the edge in x is
    unbounded-and-smooth in y, which is a shape a Gaussian mixture can represent.
    """

    def __init__(self, inner, lo, hi):
        self.inner = inner
        self.lo = numpy.asarray(lo, float)
        self.hi = numpy.asarray(hi, float)
        self.span = self.hi - self.lo

    @staticmethod
    def to_y(x, lo, span, eps=1e-9):
        u = numpy.clip((x - lo) / span, eps, 1.0 - eps)
        return numpy.log(u / (1.0 - u))

    # Keep draws strictly inside the prior. sigmoid(y) underflows for large
    # negative y, and lo + span*4e-18 rounds back to exactly lo in float64 for
    # a tight parameter like mchirp (span 2e-3), which pycbc's Uniform treats
    # as outside its bounds -- every such draw would then be rejected. An inset
    # of 1e-10 of the span is far above the float64 spacing at these values and
    # displaces a completely negligible amount of probability.
    INSET = 1e-10

    def to_x(self, y):
        u = 1.0 / (1.0 + numpy.exp(-y))
        u = numpy.clip(u, self.INSET, 1.0 - self.INSET)
        return self.lo + self.span * u

    def _log_jac(self, y):
        """ log|dy/dx| = -log(span) - log u - log(1-u), with u = sigmoid(y). """
        logu = -numpy.log1p(numpy.exp(-y))
        log1mu = -numpy.log1p(numpy.exp(y))
        return (-numpy.log(self.span)[None, :] - logu - log1mu).sum(axis=1)

    def sample(self, size):
        return self.to_x(self.inner.sample(size))

    def logpdf(self, x):
        y = self.to_y(x, self.lo, self.span)
        return self.inner.logpdf(y) + self._log_jac(y)


class GameSampler6(DummySampler):
    """games3 with inexhaustible generative proposals combined by strata.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Flat mapping file (games/games3 format).
    treefile : str, optional
        Hierarchical tree file; tiles without an entry fall back to the flat
        pool, exactly as games3 does.
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
    gen_fraction : float
        Initial share of the round budget given to the generative stratum; the
        share then adapts toward whichever stratum delivers more ESS per call.
        Set to 0 to hold the generative proposal in reserve, so that behaviour
        before exhaustion is identical to games3 and the generative stratum only
        takes over once the pools run dry.
    gen_max_centres : int
        Cap on mixture centres taken from the posterior samples.
    gen_refit : int
        Refit the generative proposal every this many rounds (1 = every round,
        0 = fit once and freeze). Refitting is safe because each sample retains
        the weight computed against the proposal in force when it was drawn.
    gen_bandwidth : float
        Extra scale on the kernel covariance. Below 1 tightens the proposal,
        which raises efficiency when the posterior is much narrower than the
        sample spread suggests.
    stretch_pairs : int
        Number of extra stretch centres (0 disables the stretch component).
    stretch_scale : float
        Upper limit of u in x = a + u (b - a); >1 extrapolates past b.
    """
    name = 'games6'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None, treefile=None,
                 loglr_region=25,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 gen_start_round=1,
                 gen_min_ess=100,
                 gen_fraction=0.0,
                 gen_max_centres=20000,
                 gen_refit=1,
                 gen_bandwidth=1.0,
                 gen_bandwidths=None,
                 gen_local_k=160,
                 gen_logit=0,
                 gen_defensive=0.0,
                 gen_defensive_scale=3.0,
                 stretch_pairs=0,
                 stretch_scale=1.0,
                 tree_start_level=0,
                 coarse_fraction=0.5,
                 min_active_points=100,
                 tile_point_cap=200000,
                 tree_ncall_budget=0,
                 gen_seed_leaves=0,
                 gen_seed_temper=-1.0,
                 gen_seed_range=4.0,
                 **kwargs):
        super().__init__(model, *args)

        self.meta = {}
        self.mapfile = mapfile
        self.treefile = treefile
        self.rounds = int(rounds)
        self.dmap = {}
        self.draw = {}

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.loglr_region = float(loglr_region)
        self.gen_start_round = int(gen_start_round)
        self.gen_min_ess = float(gen_min_ess)
        self.gen_fraction = float(gen_fraction)
        self.gen_max_centres = int(gen_max_centres)
        self.gen_refit = int(gen_refit)
        self.gen_bandwidth = float(gen_bandwidth)
        # One generative stratum per bandwidth. Because strata are combined
        # with beta proportional to ESS, the mixture automatically shifts
        # toward whichever bandwidth is working, with no switching rule to
        # tune: a bandwidth that proposes badly earns a small ESS and so a
        # small beta. Budget is also reallocated toward the better strata.
        if gen_bandwidths is None:
            self.gen_bandwidths = [self.gen_bandwidth]
        elif isinstance(gen_bandwidths, str):
            self.gen_bandwidths = [float(v) for v in
                                   gen_bandwidths.replace(',', ' ').split()]
        else:
            self.gen_bandwidths = [float(v) for v in gen_bandwidths]
        self.gen_local_k = int(gen_local_k)
        self.gen_logit = int(gen_logit)
        self.gen_defensive = float(gen_defensive)
        self.gen_defensive_scale = float(gen_defensive_scale)
        self.stretch_pairs = int(stretch_pairs)
        self.stretch_scale = float(stretch_scale)
        self.tree_start_level = int(tree_start_level)
        # If at least this fraction of start-level nodes clears the cut,
        # skip the descent and use those nodes as leaves -- see the long
        # comment at the call site. 0 restores the always-descend
        # behaviour. Only active when tree_start_level > 0.
        self.coarse_fraction = float(coarse_fraction)
        # Fall back to start-level tiles only if the cut leaves fewer
        # points than this IN TOTAL. Set low on purpose. An earlier value
        # of 2000 fired on P-P injection 19, which had 9 leaves holding
        # 1484 points -- a legitimately concentrated active set at SNR
        # 34, where a tight cut is correct -- and swapping those
        # well-placed fine leaves for coarse tiles cost 1.31% -> 0.07%
        # efficiency, an 18x regression. Saturation alone is survivable
        # (the budget passes to the generative stratum); what is not
        # survivable is a pool with nothing in it, e.g. injection 146's
        # single point.
        self.min_active_points = int(min_active_points)
        # Cap on points kept per start-level tile when a coarse path is
        # taken. Their subtrees can be the entire map; a run draws far
        # less. The tile's TRUE population is still used for its prior
        # mass, so capping changes only what is available to draw.
        self.tile_point_cap = int(tile_point_cap)
        # Hard cap on likelihood calls spent RESOLVING which leaves are
        # active; 0 (default) disables it. DO NOT enable without reading
        # the warning in `_find_active_leaves` -- it biases posteriors.
        self.tree_ncall_budget = int(tree_ncall_budget)
        self.gen_seed_leaves = int(gen_seed_leaves)
        self.gen_seed_temper = float(gen_seed_temper)
        self.gen_seed_range = float(gen_seed_range)
        self.leaf_loglr = {}
        # seeded from the legacy global RNG, which pycbc_inference --seed sets,
        # so generative draws are reproducible alongside the shuffles that the
        # discrete-pool path uses
        self._rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 31 - 1))
        self.ncalls = 0
        self.tree_ncalls = 0
        self.stratum_calls = {'pool': 0, 'gen': 0}

    # ---------------- discrete pool stratum (unchanged from games3) ---------

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
        """ One round of games3's discrete-pool draw. Returns
        (psamp, loglr, logw_unnormalised, bin_id) or raises OutOfSamples.

        `lengths` is how many points a bin physically HOLDS -- the cap on
        how many can be drawn. `prior_mass` is how much PRIOR MASS it
        represents, which is what belongs in the importance weight. They
        are the same thing whenever a bin stores every point assigned to
        it, which is why one array served both roles originally; they
        differ as soon as a bin stores a subsample of a larger
        population, as the coarse-tile paths below do (a start-level tile
        can cover millions of points -- loading them all would need GBs
        and there is no reason to, since a run only draws tens of
        thousands). Keeping them separate is also exactly what per-leaf
        refinement needs: oversample a thin leaf, store the extra points,
        and leave its prior mass at the unbiased value.
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

        psamp = FieldArray(total, dtype=self.dtype)
        pweight = numpy.zeros(total)
        bin_id = numpy.zeros(total, dtype=int)
        j = 0
        for i, (c, w) in enumerate(zip(drawcount, drawweight)):
            bdraw = self.draw_samples_from_bin(node_idx[i], c)
            psamp[j:j+len(bdraw)] = FieldArray.from_records(bdraw,
                                                            dtype=self.dtype)
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        loglr = self._likelihoods(psamp, stratum='pool')
        logw = loglr + numpy.log(prior_mass[bin_id]) - pweight
        self._dump('pool', loglr=loglr, logw=logw, bin_id=bin_id,
                   binlen=lengths[bin_id], pweight=pweight,
                   bin_weight=bin_weight, drawcount=drawcount,
                   x=numpy.column_stack([psamp[p] for p in
                                         self.model.variable_params]))
        return psamp, loglr, logw, bin_id

    # ---------------- temporary diagnostics (env-var gated) -----------------

    def _dump(self, tag, **kwargs):
        """ Write per-round weight forensics if GAMES6_DUMP is set to a path
        prefix. Purely diagnostic: no effect on the sampler when unset.
        """
        prefix = os.environ.get('GAMES6_DUMP')
        if not prefix:
            return
        fn = '%s_%s_r%02i.npz' % (prefix, tag.replace(':', '-'),
                                  getattr(self, '_round', 0))
        kwargs['names'] = numpy.array(list(self.model.variable_params))
        numpy.savez(fn, **kwargs)
        logging.info('dumped %s diagnostics to %s', tag, fn)

    # ---------------- generative stratum -----------------------------------

    def gen_round(self, proposal, budget, stratum='gen'):
        """ One round from the generative proposal. Out-of-prior draws are
        dropped before any likelihood call.
        """
        names = list(self.model.variable_params)
        xall = proposal.sample(budget)
        keep = self._in_prior(xall, names)
        if keep.sum() == 0:
            self._dump(stratum + '-rejected', x=xall,
                       is_logit=int(isinstance(proposal, BoundedProposal)),
                       lo=getattr(proposal, 'lo', numpy.zeros(0)),
                       hi=getattr(proposal, 'hi', numpy.zeros(0)))
            return None
        x = xall[keep]
        logq = proposal.logpdf(x)

        psamp = FieldArray(len(x), dtype=self.dtype)
        for k, p in enumerate(names):
            psamp[p] = x[:, k]
        loglr = self._likelihoods(psamp, stratum=stratum)
        self._dump(stratum, x=x, loglr=loglr, logq=logq, budget=budget,
                   nkept=int(keep.sum()),
                   centres=getattr(proposal, 'centres',
                                   getattr(getattr(proposal, 'inner', None),
                                           'centres', numpy.zeros(0))),
                   cweights=getattr(proposal, 'weights',
                                    getattr(getattr(proposal, 'inner', None),
                                            'weights', numpy.zeros(0))),
                   kcov=getattr(proposal, 'cov',
                                getattr(getattr(proposal, 'inner', None),
                                        'cov', numpy.zeros(0))),
                   is_logit=int(isinstance(proposal, BoundedProposal)),
                   lo=getattr(proposal, 'lo', numpy.zeros(0)),
                   hi=getattr(proposal, 'hi', numpy.zeros(0)))
        # uniform prior => its density is a constant that cancels in the
        # stratum's own self-normalisation
        return psamp, loglr, loglr - logq

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
            self.stratum_calls[stratum] += len(args)
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

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as mapfile:
            self.dtype = mapfile['map']['0'].dtype

        treefile = h5py.File(self.treefile, 'r') if self.treefile else None

        # Where to start evaluating. The likelihood is evaluated at each
        # start node's REPRESENTATIVE and everything more than
        # loglr_region below the best is discarded, so the cut is only as
        # good as that representative proxies the likelihood across its
        # whole tile. A coarse representative is a poor proxy -- one at
        # match m sits roughly (1 - m^2) * rho^2/2 below a perfect one,
        # ~46 nats for m=0.7 at SNR 13.5 -- so a window tuned for a fine
        # bank can throw away subtrees that do hold posterior mass.
        # Coarse levels are still wanted for BUILDING the tree (they make
        # routing a draw cheap), so this lets construction and evaluation
        # use different depths instead of forcing one compromise.
        start_nodes = None
        if treefile is not None and self.tree_start_level > 0:
            start_nodes = []
            for k in sorted(treefile['tree'], key=int):
                start_nodes.extend(self._nodes_at_depth(
                    treefile['tree'][k], self.tree_start_level))
            args = [{p: float(n.attrs[f'param_{p}'])
                     for p in self.model.variable_params}
                    for n in start_nodes]
            logging.info('starting evaluation at tree depth %i: %i nodes '
                         '(vs %i at depth 0)', self.tree_start_level,
                         len(start_nodes), len(treefile['tree']))
        else:
            with h5py.File(self.mapfile, 'r') as mapfile:
                bparams = {p: mapfile['bank'][p][:]
                           for p in self.variable_params}
                num_nodes = len(bparams[list(bparams.keys())[0]])
            args = [{p: bparams[p][i] for p in self.model.variable_params}
                    for i in range(num_nodes)]

        logging.info('Calculating likelihood at nodes')
        node_loglrs = numpy.array(list(tqdm.tqdm(
            self.pool.imap(call_likelihood, args), total=len(args))))
        self.meta['setup_ncalls'] = len(args)
        bound = node_loglrs[~numpy.isnan(node_loglrs)].max() - self.loglr_region
        passed = numpy.where((node_loglrs > bound)
                             & ~numpy.isnan(node_loglrs))[0]

        # SHORT-CIRCUIT THE DESCENT WHEN THE CUT CANNOT DISCRIMINATE.
        #
        # The descent evaluates the likelihood at every child of every
        # node it enters, so when nothing prunes it enumerates the WHOLE
        # tree: measured maxima of 445601 (BNS) and 149959 (BBH) are
        # exactly level2+level3 of those trees, i.e. every internal node
        # plus every leaf. That is 4-15x the sampling budget, and one
        # such run took 3.5 h against ~10 min typical.
        #
        # It happens at LOW SNR, and the threshold follows from the cut
        # itself: a node at match m sits (1-m^2) rho^2/2 nats below the
        # peak, so pruning needs rho > sqrt(2*loglr_region/(1-m^2)) --
        # about 12.6 for the 0.9-match nodes of level 1, and NOTHING
        # prunes below rho ~ sqrt(2*loglr_region). Measured: corr(log
        # tree_ncalls, SNR) = -0.81 (BNS) and -0.70 (BBH); the runs that
        # enumerate everything have mean SNR 6.9 against 17.0 for the
        # rest, and they are 21 of the 29 worst-efficiency runs.
        #
        # The cut is not misbehaving -- at rho ~ 5 the posterior really
        # does span most of the prior, so nearly every leaf really is
        # active. The waste is paying one likelihood call per node to
        # rediscover it. So when the passing FRACTION is high, stop and
        # use the start-level nodes themselves as leaves, each carrying
        # its whole subtree's points and its own already-paid-for rep
        # loglr. As rho -> 0 this is not even an approximation: the
        # weights become equal and the subtree clouds union to the prior,
        # which is exactly the target.
        coarse = (start_nodes is not None and self.coarse_fraction > 0
                  and len(passed) >= self.coarse_fraction * len(args))
        if coarse:
            logging.info('%i of %i start nodes passed (>= %.0f%%): the cut '
                         'cannot discriminate, so using start-level tiles '
                         'directly instead of descending. This is the '
                         'low-SNR regime where the posterior spans most of '
                         'the prior anyway.', len(passed), len(args),
                         100 * self.coarse_fraction)

        logging.info('...resolving tree leaves for %i passed tiles',
                     len(passed))
        active_lengths, key = [], 0
        self.leaf_loglr = {}
        tile_mass = {}
        with h5py.File(self.mapfile, 'r') as mapfile:
            for tile in passed:
                leaves = None
                if coarse:
                    pts, tot = self._subtree_points(start_nodes[tile],
                                                    self.tile_point_cap)
                    leaves = [(pts, node_loglrs[tile], tot)] if len(pts) else []
                elif start_nodes is not None:
                    leaves = self._find_active_leaves(start_nodes[tile],
                                                      bound)
                elif treefile is not None and str(tile) in treefile['tree']:
                    leaves = self._find_active_leaves(
                        treefile['tree'][str(tile)], bound)
                if not leaves and start_nodes is None:
                    mp = mapfile['map'][str(tile)][:]
                    leaves = [(mp, node_loglrs[tile], len(mp))]
                for item in (leaves or []):
                    pts, ll = item[0], item[1]
                    self.dmap[key] = pts
                    self.leaf_loglr[key] = ll
                    active_lengths.append(len(pts))
                    if len(item) > 2:
                        tile_mass[key] = item[2]
                    key += 1
        # NB: the tree file stays open until after the starved-set
        # fallback below, which reads subtrees through `start_nodes`.
        # Closing here left those as handles into a closed file and the
        # fallback died with an h5py "invalid identifier" KeyError.
        active_ids = numpy.arange(key)
        lengths = numpy.array(active_lengths)
        prior_mass = numpy.array([tile_mass.get(i, active_lengths[i])
                                  for i in range(key)], dtype=float)
        logging.info('...resolved to %i active leaves (%i points)',
                     key, int(lengths.sum()) if key else 0)

        # STARVED ACTIVE SET. The opposite failure to the one above: the
        # cut prunes so hard that almost nothing survives, and the pool
        # then saturates immediately (`pool_round` caps drawcount at leaf
        # size), which breaks the w ~ L * n_bin / c_bin cancellation.
        #
        # This is not hypothetical. Over 291 wide-space P-P runs the
        # efficiency residual at SNR > 15, after detrending the SNR
        # dependence, correlates with log(tree_ncalls) at +0.31: the
        # smallest-active-set quartile has median efficiency 13.3%
        # against 31.9% for the largest, at comparable SNR. And one run
        # (idx 146, SNR 24) resolved to a single point and reported
        # ESS = 1 from ncalls = 1 -- a "100% efficient" run that produced
        # one sample.
        #
        # Falling back to the start-level tiles keeps the run alive with
        # a coarser but well-populated pool, which is strictly better
        # than a pool of one. Costs no likelihood calls: the tiles and
        # their reps were already evaluated.
        starved = key and lengths.sum() < self.min_active_points
        if (starved or not key) and start_nodes is not None and not coarse:
            logging.warning('active set starved (%i leaves, %i points < %i): '
                            'the cut pruned nearly everything, which '
                            'saturates the pool on the first round. '
                            'Falling back to start-level tiles.',
                            key, int(lengths.sum()) if key else 0,
                            self.min_active_points)
            # Widen the cut until the pool is usable. Falling back to the
            # start-level tiles is not enough on its own: the tile that
            # holds the signal can itself be nearly empty. Measured on
            # P-P injection 146 (SNR ~24) -- exactly ONE of 1104 tiles
            # cleared the cut and its whole subtree held 8 points, so the
            # run returned ESS 1. The tessellation is by MATCH, so a
            # tile's population follows the local prior volume and some
            # corners really are that thin; a signal that lands in one is
            # the adverse case a P-P test is supposed to catch, and it
            # did (1 run in 291).
            #
            # Widening costs NO likelihood calls -- every start node was
            # already evaluated during setup, so admitting more of them
            # is free. It cannot bias the result either: extra tiles are
            # low-likelihood and importance weighting gives them
            # correspondingly low weight. It trades a little efficiency
            # for not returning a one-sample posterior.
            region = self.loglr_region
            for _ in range(6):
                bound_w = numpy.nanmax(node_loglrs) - region
                sel = numpy.where((node_loglrs > bound_w)
                                  & ~numpy.isnan(node_loglrs))[0]
                active_lengths, key = [], 0
                self.dmap, self.leaf_loglr = {}, {}
                tile_mass = {}
                for tile in sel:
                    pts, tot = self._subtree_points(start_nodes[tile],
                                                    self.tile_point_cap)
                    if not len(pts):
                        continue
                    self.dmap[key] = pts
                    self.leaf_loglr[key] = node_loglrs[tile]
                    active_lengths.append(len(pts))
                    tile_mass[key] = tot
                    key += 1
                if sum(active_lengths) >= self.min_active_points:
                    break
                region *= 2.0
                logging.warning('still starved (%i points from %i tiles); '
                                'widening loglr_region to %.0f',
                                sum(active_lengths), len(sel), region)
            active_ids = numpy.arange(key)
            lengths = numpy.array(active_lengths)
            prior_mass = numpy.array([tile_mass.get(i, active_lengths[i])
                                      for i in range(key)], dtype=float)
            logging.info('...fallback gave %i tiles (%i points)',
                         key, int(lengths.sum()) if key else 0)

        if treefile is not None:
            treefile.close()

        # Round-1 draw probability. Plain leaf size is prior VOLUME only:
        # the representative likelihood computed during the cut is used
        # there as a pass/fail threshold and its magnitude then thrown
        # away, so the first round spreads draws by volume and only the
        # adaptive update below (weight[j] += v) starts concentrating them
        # from round 2. Seeding with exp(loglr_rep) * volume uses the
        # information already paid for. This cannot bias the estimate --
        # logw = loglr + log(lengths) - pweight carries log(drawcount), so
        # the importance weights self-correct for any draw probability --
        # it only changes their variance. The adaptive tuning is unchanged
        # and still runs on top of this.
        if self.gen_seed_leaves:
            ll = numpy.array([self.leaf_loglr.get(k, numpy.nan)
                              for k in range(key)], dtype=float)
            ok = numpy.isfinite(ll)
            if ok.any():
                # TEMPERING. exp(loglr_rep) is the likelihood AT the
                # representative, but a leaf's actual weight is the
                # likelihood INTEGRATED over the leaf. Those agree only
                # while the likelihood barely varies inside a leaf. The
                # variation is (1 - m^2) * rho^2/2, so it grows as rho^2:
                # ~1.8 nats at SNR 13.5 but ~9 at SNR 30 and ~95 at SNR
                # 98, i.e. a round-1 weight ratio of 6:1 rising to
                # 7750:1 and beyond. Untempered, that piles every draw
                # into a handful of leaves and drains them immediately --
                # measured, it cost 20.8% -> 12.7% efficiency at SNR 30,
                # where the discrete pool is still the workhorse, while
                # helping at SNR 50/98 where the pool is nearly dead
                # anyway. The integral is far flatter than the point
                # value, so temper the exponent to approximate it.
                spread = float(ll[ok].max() - ll[ok].min())
                if self.gen_seed_temper > 0:
                    temper = self.gen_seed_temper
                else:
                    # auto: hold the effective dynamic range near
                    # gen_seed_range nats regardless of SNR
                    temper = max(1.0, spread / max(self.gen_seed_range,
                                                   1e-6))
                w = lengths.astype(float).copy()
                w[ok] *= numpy.exp((ll[ok] - ll[ok].max()) / temper)
                if w.sum() > 0:
                    weight = w / w.sum()
                    logging.info('seeded round-1 leaf weights: '
                                 'exp((loglr_rep-max)/%.2f) * volume '
                                 '(loglr spread %.1f nats over %i leaves)',
                                 temper, spread, int(ok.sum()))
                else:
                    weight = lengths / lengths.sum()
            else:
                weight = lengths / lengths.sum()
        else:
            weight = lengths / lengths.sum()
        # per-stratum accumulators
        gen_keys = [f'gen:bw{bw:g}' for bw in self.gen_bandwidths]
        strata = {'pool': {'samp': None, 'loglr': None, 'logw': None}}
        for k in gen_keys:
            strata[k] = {'samp': None, 'loglr': None, 'logw': None}
        self.stratum_calls = {k: 0 for k in ['pool'] + gen_keys}
        proposals = {}
        if self.gen_seed_leaves:
            seed = self._seed_proposal_from_leaves()
            if seed is not None:
                proposals = {gen_keys[0]: seed}
        # per-stratum share of the generative budget, adapted by ESS per call
        share = {k: 1.0 / len(gen_keys) for k in gen_keys}
        pool_dead = False
        frac = self.gen_fraction

        for rnd in range(1, self.rounds + 1):
            self._round = rnd
            gen_budget = 0
            if proposals:
                gen_budget = int(self.target_likelihood_calls
                                 * (1.0 if pool_dead else frac))
            pool_budget = self.target_likelihood_calls - gen_budget

            pool_gain = None
            gen_gains = {}
            if pool_budget > 0 and not pool_dead:
                try:
                    ps, pl, pw, bid = self.pool_round(
                        weight / weight.sum(), active_ids, lengths,
                        pool_budget, prior_mass=prior_mass)
                    before = self._ess(strata['pool']['logw'])
                    strata['pool'] = self._accumulate(strata['pool'],
                                                      ps, pl, pw)
                    after = self._ess(strata['pool']['logw'])
                    pool_gain = (after - before) / max(len(pw), 1)
                    # games3's adaptive leaf reweighting, unchanged.
                    # A per-leaf shrinkage estimator of E[L] was tried here
                    # (rep as weak prior, overwritten by observed draws) on
                    # the theory that L(rep) is a poor one-sample estimator
                    # once the within-leaf likelihood spread grows. It
                    # measured NEUTRAL -- tempering beat it everywhere and
                    # it added nothing on top -- so it was removed to keep
                    # this path simple. See commit 4915a9340 for the code
                    # and the numbers.
                    wn = numpy.exp(pw - logsumexp(pw))
                    for j, v in zip(bid, wn):
                        weight[j] += v
                except OutOfSamples:
                    logging.info('pool exhausted at round %i; handing its '
                                 'budget to the generative proposals', rnd)
                    pool_dead = True
                    if not proposals:
                        proposals = self._fit_proposals(strata, gen_keys)
                    if not proposals:
                        logging.info('no generative proposal available; stop')
                        break
                    gen_budget = self.target_likelihood_calls

            if gen_budget > 0 and proposals:
                live = [k for k in gen_keys if k in proposals]
                tot_share = sum(share[k] for k in live) or 1.0
                for k in live:
                    b = int(gen_budget * share[k] / tot_share)
                    if b <= 0:
                        continue
                    got = self.gen_round(proposals[k], b, stratum=k)  # noqa
                    if got is None:
                        continue
                    gs, gl, gw = got
                    # one stratum per (bandwidth, round): each round used a
                    # different refitted proposal, so its samples are only
                    # commensurate with themselves
                    rk = f'{k}:r{rnd}'
                    strata[rk] = {'samp': gs, 'loglr': gl, 'logw': gw}
                    self.stratum_calls.setdefault(rk, 0)
                    self.stratum_calls[rk] = len(gw)
                    self.stratum_calls[k] = self.stratum_calls.get(k, 0)
                    gen_gains[k] = self._ess(gw) / max(len(gw), 1)
                if gen_gains:
                    tot = sum(max(v, 0) for v in gen_gains.values())
                    if tot > 0:
                        for k in gen_keys:
                            g = max(gen_gains.get(k, 0.0), 0.0) / tot
                            share[k] = 0.5 * share[k] + 0.5 * g
                        ssum = sum(share.values()) or 1.0
                        for k in gen_keys:
                            share[k] = max(share[k] / ssum, 0.02)

            # shift budget toward whichever stratum yields more ESS per call.
            # gen_fraction <= 0 means "hold the generative proposal in reserve":
            # the pool keeps the entire budget, so behaviour before exhaustion
            # is identical to games3, and the generative stratum only runs once
            # the pool dies and hands over its budget.
            best_gen = max(gen_gains.values()) if gen_gains else None
            if self.gen_fraction > 0 and pool_gain is not None \
                    and best_gen is not None:
                gen_gain = best_gen
                tot = max(pool_gain, 0) + max(gen_gain, 0)
                if tot > 0:
                    frac = 0.5 * frac + 0.5 * (max(gen_gain, 0) / tot)
                    frac = min(max(frac, 0.05), 0.95)

            # Refit the generative proposal as posterior information
            # accumulates. This is safe because every sample already carries
            # the logw computed against the proposal that was in force when it
            # was drawn -- weights are never recomputed against a later
            # proposal, which would bias the estimate. Refitting is the whole
            # point of an adaptive scheme: a proposal frozen at the first fit
            # can only ever be as good as the posterior estimate available at
            # that round.
            if rnd >= self.gen_start_round and self.gen_refit > 0 and \
                    (rnd - self.gen_start_round) % self.gen_refit == 0:
                new = self._fit_proposals(strata, gen_keys)
                if new:
                    proposals = new

            ess_p = self._ess(strata['pool']['logw'])
            gkeys = [q for q in strata if q.startswith('gen')]
            gsum = sum(self._ess(strata[q]['logw']) for q in gkeys)
            detail = f'{len(gkeys)} sub-strata'
            logging.info('round %i: ESS pool=%.1f gen[%s] total=%.1f '
                         'ncalls=%i', rnd, ess_p, detail, ess_p + gsum,
                         self.ncalls)

        self._finalise(strata)

    def _accumulate(self, st, samp, loglr, logw):
        if st['samp'] is None:
            return {'samp': samp, 'loglr': loglr, 'logw': logw}
        return {'samp': numpy.concatenate([samp, st['samp']]),
                'loglr': numpy.concatenate([loglr, st['loglr']]),
                'logw': numpy.concatenate([logw, st['logw']])}

    def _seed_proposal_from_leaves(self, max_points=200000):
        """ Fit a generative proposal to the ACTIVE LEAVES, before any
        posterior samples exist.

        The setup cut already evaluates the likelihood at each leaf's
        representative -- we pay for those calls to decide the descent and
        then normally discard them. But a leaf's posterior mass is
        approximately exp(loglr_rep) times its prior volume, and its point
        count is proportional to that volume (the points are uniform prior
        draws), so weighting every point in leaf i by exp(loglr_i) turns
        the union of active leaves into an estimate of the posterior --
        the same estimate that makes the discrete-pool stratum work at
        all. If a KDE fitted to accumulated posterior samples is a good
        proposal later, one fitted to this cloud should be a good proposal
        immediately, without waiting for the pool to deliver gen_min_ess.

        That matters most where games6 is weakest: at high SNR the pool
        exhausts after a single round with ESS ~10-20, so the generative
        stratum barely engages under the accumulate-then-fit rule.

        Built as a synthetic stratum and handed to the existing
        _fit_proposal, so it goes through exactly the same weighted
        centre-subsampling and kernel construction as every other fit.
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
        fa = FieldArray(len(samp), dtype=self.dtype)
        for p in names:
            fa[p] = samp[p]
        stratum = {'samp': fa, 'loglr': logw, 'logw': logw}
        prop = self._fit_proposal({'seed': stratum}, self.gen_bandwidth)
        if prop is not None:
            logging.info('seeded generative proposal from %i active leaves '
                         '(%i points, subsampled to %i)', len(keys), total,
                         len(samp))
        return prop

    def _fit_proposals(self, strata, gen_keys):
        out = {}
        for bw, k in zip(self.gen_bandwidths, gen_keys):
            p = self._fit_proposal(strata, bw)
            if p is not None:
                out[k] = p
        if out:
            logging.info('fitted %i generative proposals: %s', len(out),
                         ', '.join(sorted(out)))
        return out

    def _prior_bounds(self):
        names = list(self.model.variable_params)
        lo = numpy.full(len(names), -numpy.inf)
        hi = numpy.full(len(names), numpy.inf)
        try:
            dists = self.model.prior_distribution.distributions
        except AttributeError:
            return None
        for dist in dists:
            b = getattr(dist, 'bounds', None) or {}
            for k, v in b.items():
                if k in names:
                    i = names.index(k)
                    lo[i], hi[i] = float(v[0]), float(v[1])
        if not numpy.all(numpy.isfinite(lo) & numpy.isfinite(hi)):
            return None
        return lo, hi

    def _fit_proposal(self, strata, bandwidth=None):
        """ Fit the generative proposal to the accumulated posterior samples
        from every stratum, resampled to equal weight.
        """
        names = list(self.model.variable_params)
        # Weight each stratum's contribution by its OWN ESS when pooling for
        # the fit, unconditionally. Normalising every stratum to sum to 1 and
        # concatenating gives each equal total mass regardless of quality, so
        # a weak stratum drags the fit toward wherever it happens to have
        # drawn -- the same defect as games.py's per-round self-normalisation,
        # one function away. This matches how _finalise already combines
        # strata (beta proportional to ESS, the inverse-variance choice).
        xs, ws = [], []
        _diag = []
        for kk, st in strata.items():
            if st['samp'] is None:
                continue
            lw = st['logw']
            wv = numpy.exp(lw - logsumexp(lw))
            e = 1.0 / (wv ** 2).sum()
            _diag.append((kk, len(wv), e))
            if e <= 0:
                continue
            xs.append(numpy.column_stack([st['samp'][p] for p in names]))
            ws.append(wv * e)
        if not xs:
            return None
        x = numpy.vstack(xs)
        w = numpy.concatenate(ws)
        if w.sum() <= 0:
            return None
        w = w / w.sum()
        ess = 1.0 / (w ** 2).sum()
        logging.debug('fit blocks=%s sum_stratum_ess=%.1f pooled_ess=%.1f',
                      ';'.join('%s:n=%i,ess=%.1f' % d for d in _diag),
                      sum(d[2] for d in _diag), ess)
        if ess < self.gen_min_ess:
            logging.debug('accumulated ESS %.1f below gen_min_ess %.1f; '
                          'not fitting generative proposal yet',
                          ess, self.gen_min_ess)
            return None
        # A properly WEIGHTED mixture over the samples: every sample is a
        # centre and carries its own posterior weight as its mixture weight.
        # Three constructions were tried for which raw draws become centres
        # when there are more than gen_max_centres of them. Resampling by
        # weight WITH replacement duplicates a handful of points at low ESS
        # (ESS 19 gives ~19 distinct centres), making the covariance
        # degenerate. Taking the top-N by weight concentrates centres at the
        # peak, which under-disperses the proposal -- fatal for importance
        # sampling, since the target's tails then have almost no proposal
        # density and pi/q explodes there. Uniform random selection (ignoring
        # weight) avoids both, but is blind to information content: as the
        # raw pool keeps growing past the cap round over round, each
        # individual high-weight point's chance of survival keeps shrinking
        # even as accumulated ESS climbs -- measured directly (SNR 30, fixed
        # gen_max_centres): the fit-input accumulated ESS rose from ~2700 to
        # ~7900 across rounds while the resulting round's own efficiency
        # *fell* from 27.5% to 20.0%, the opposite of what more information
        # should buy.
        #
        # Weighted sampling WITHOUT replacement gets both properties at once:
        # numpy's choice(replace=False, p=w) can never pick the same point
        # twice (so it cannot reproduce the with-replacement duplication
        # failure), while still preferring higher-weight points, so the kept
        # centres retain most of the pool's weight regardless of how peaked
        # it is (measured: a synthetic pool with weight ~ U(0,1)^8, i.e. very
        # peaked, keeps 99.6% of total weight when subsampled 40000 -> 20000,
        # versus ~50% for the uniform selection it replaces).
        n = min(self.gen_max_centres, len(x))
        if n < len(x):
            # w is already normalised to sum to 1 above
            order = self._rng.choice(len(x), size=n, replace=False, p=w)
        else:
            order = numpy.arange(len(x))
        centres = x[order]
        cw = w[order]
        if cw.sum() <= 0:
            return None
        logging.debug('centres=%i centre_w_ess=%.1f (of %i drawn from %i)',
                      len(cw), cw.sum() ** 2 / (cw ** 2).sum(), n, len(x))
        # weighted covariance over ALL samples, which uses the full posterior
        # information rather than only the retained centres
        mu = (w[:, None] * x).sum(0)
        xc = x - mu
        cov = (w[:, None] * xc).T @ xc / max(1.0 - (w ** 2).sum(), 1e-12)
        if not numpy.all(numpy.isfinite(cov)):
            return None
        # Kernel choice and coordinate transform are ORTHOGONAL, and must be
        # composed rather than raced. An earlier version returned from the logit
        # branch before the gen_local_k test could run, and since gen_logit
        # defaults to 1 that made the local-covariance kernel dead code -- it was
        # then wrongly reported as "tried and ineffective". Build the kernel
        # first, then optionally wrap it in the logit map.
        bw = self.gen_bandwidth if bandwidth is None else bandwidth
        bnd = self._prior_bounds() if self.gen_logit else None

        fit_x, lo, hi = centres, None, None
        if bnd is not None:
            lo, hi = bnd
            fit_x = BoundedProposal.to_y(centres, lo, hi - lo)
            if not numpy.all(numpy.isfinite(fit_x)):
                fit_x, lo, hi = centres, None, None

        try:
            if self.gen_local_k > 0:
                inner = LocalCovarianceProposal(
                    fit_x, self._rng, k=self.gen_local_k, scale=bw,
                    weights=cw, defensive_frac=self.gen_defensive,
                    defensive_scale=self.gen_defensive_scale)
            else:
                fcov = numpy.cov(fit_x.T, aweights=cw) if lo is not None \
                    else cov
                if not numpy.all(numpy.isfinite(fcov)):
                    return None
                inner = GaussianMixtureProposal(
                    fit_x, fcov * bw ** 2, self._rng,
                    stretch_pairs=self.stretch_pairs,
                    stretch_scale=self.stretch_scale, weights=cw,
                    defensive_frac=self.gen_defensive,
                    defensive_scale=self.gen_defensive_scale)
        except numpy.linalg.LinAlgError:
            logging.info('generative proposal fit failed (singular kernel)')
            return None

        prop = inner if lo is None else BoundedProposal(inner, lo, hi)
        kern = ('local%d' % self.gen_local_k if self.gen_local_k > 0
                else 'global')
        logging.info('fitted proposal: kernel=%s space=%s bw=%.2f centres=%i '
                     'accumulated ESS=%.1f', kern,
                     'logit' if lo is not None else 'raw', bw, inner.n, ess)
        return prop

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

        # Setup and tree-cut calls are deliberately NOT folded into ncalls:
        # both are one-time costs paid before/while resolving which leaves
        # are active, independent of how many posterior samples are then
        # drawn, so they amortize away at the ESS this sampler targets.
        # Reported separately because they are still a useful diagnostic --
        # a tree whose root level is too coarse prunes weakly and descends
        # a large fraction of its nodes, which shows up here.
        self.meta['ncalls'] = self.ncalls
        self.meta['ess'] = total_ess
        self.meta['setup_ncalls'] = self.meta.get('setup_ncalls', 0)
        self.meta['tree_ncalls'] = self.tree_ncalls

        nout = max(int(total_ess), 1)
        outs, ols = [], []
        for (k, v), b in zip(parts, beta):
            take = int(round(nout * b))
            if take <= 0:
                continue
            lw = v['logw']
            p = numpy.exp(lw - logsumexp(lw))
            idx = numpy.random.choice(len(p), size=take, replace=True, p=p)
            outs.append(v['samp'][idx])
            ols.append(v['loglr'][idx])
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

    def _subtree_points(self, group, cap=None):
        """ Points under `group`, subsampled to `cap`, plus the TRUE total.

        The counterpart to `_find_active_leaves` for the low-SNR case:
        that one pays a likelihood call per node to decide what to keep,
        this one keeps what is below a node that already passed. Pure
        HDF5 reads, so the cost is IO rather than likelihood calls.

        The cap is not optional in practice. When the cut fails to prune
        at all, EVERY start-level tile passes, and their subtrees are the
        whole map -- 100M points, ~2.4 GB, for a run that will draw a few
        tens of thousands. (Found the hard way: the uncapped version put
        one process at 4.1 GB RSS and made the runs slower than the
        descent it was meant to replace.)

        Returns (points, total) so the caller can keep prior mass and
        physical size apart: `total` is the leaf's share of the prior,
        `len(points)` is only how many are available to draw.
        """
        chunks, total = [], 0

        def gather(g):
            nonlocal total
            if g.attrs['is_leaf']:
                n = int(g.attrs['n_points'])
                total += n
                if n:
                    pts = numpy.zeros(n, dtype=self.dtype)
                    for p in self.variable_params:
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

    def _find_active_leaves(self, tree_root, loglr_bound):
        """ Prune the subtree, evaluating the likelihood at each child's
        representative and descending only where it clears the bound.

        Those evaluations are real likelihood calls and are counted into
        self.tree_ncalls. They used to be counted NOWHERE, which silently
        inflated the reported efficiency (ESS/ncalls) of every tree-based
        run -- and by a tree-shape-dependent amount, so it also made runs
        on differently-shaped trees non-comparable: a tree with few coarse
        root tiles hides many nodes below each one and pays far more here
        than a tree with many fine root tiles.

        BUDGET (`tree_ncall_budget`, DEFAULT OFF -- it BIASES posteriors).
        The idea was that resolving WHERE to sample should not cost more
        than the sampling itself, capping the descent and taking any
        still-unresolved node wholesale. It targets a real gap:
        `coarse_fraction` is all-or-nothing on how many START nodes pass,
        so P-P injection 5 (SNR 8.2) had under half pass, descended, and
        still paid 146113 calls -- five times its sampling budget -- for
        3.06% efficiency. Capping it cut that to 40194 and lifted
        efficiency to 3.82%.

        It is nonetheless WRONG. Against a dynesty reference on that same
        injection, the ordinary descent matches on every parameter
        (offsets <= 0.087 posterior widths) while the budgeted run is off
        by up to **0.81 widths** on mchirp and 0.58 on lambda1 -- a pure
        location bias, since the widths agree to 0.90-1.10 either way.

        The mechanism is the difference from `coarse_fraction`, which IS
        verified: that path takes EVERY tile wholesale, so every bin has
        the same resolution. A budget instead leaves a MIX of fine leaves
        and wholesale subtrees, and a wholesale bin carries one
        representative loglr for a region many leaves wide. Leaf-seeded
        round-1 weighting (`gen_seed_leaves`) then treats incomparable
        bins as comparable. **Uniform bin resolution appears to be the
        thing that matters; mixing scales is what breaks it.**

        Returns (points, loglr, prior_mass) triples. `prior_mass` differs
        from `len(points)` only for wholesale-taken subtrees, which are
        capped for memory -- see `_subtree_points`.
        """
        results = []
        # OFF by default -- see the warning in the docstring above.
        budget = self.tree_ncall_budget

        def descend(group, ll_here):
            if group.attrs['is_leaf']:
                n = int(group.attrs['n_points'])
                pts = numpy.zeros(n, dtype=self.dtype)
                for p in self.variable_params:
                    pts[p] = group[f'points_{p}'][:]
                # carry the likelihood at this leaf's own representative:
                # it was already paid for to decide the descent, and it is
                # what makes the leaf-weighted cloud an estimate of the
                # posterior (see _seed_proposal_from_leaves)
                results.append((pts, ll_here, n))
                return
            if budget and self.tree_ncalls >= budget:
                pts, tot = self._subtree_points(group, self.tile_point_cap)
                if len(pts):
                    results.append((pts, ll_here, tot))
                return
            for c in group['children']:
                child = group['children'][c]
                params = {p: float(child.attrs[f'param_{p}'])
                          for p in self.variable_params}
                ll = call_likelihood(params)
                self.tree_ncalls += 1
                if numpy.isfinite(ll) and ll > loglr_bound:
                    descend(child, ll)

        descend(tree_root, numpy.nan)
        return results
