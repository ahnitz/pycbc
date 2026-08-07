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
"""
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
                 ridge=1e-3):
        from scipy.spatial import cKDTree
        self.rng = rng
        self.centres = numpy.atleast_2d(centres)
        self.n, self.dim = self.centres.shape
        if weights is None:
            weights = numpy.full(self.n, 1.0 / self.n)
        self.weights = numpy.asarray(weights, float)
        self.weights /= self.weights.sum()
        self.logw = numpy.log(numpy.maximum(self.weights, 1e-300))

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

    def sample(self, size):
        idx = self.rng.choice(self.n, size=size, p=self.weights)
        z0 = self.rng.normal(size=(size, self.dim))
        z = self._cz[idx] + numpy.einsum('nij,nj->ni', self.chol[idx], z0)
        return z @ self._Winv

    def logpdf(self, x):
        z = x @ self._W
        out = numpy.full(len(z), -numpy.inf)
        for i in range(self.n):
            d = z - self._cz[i]
            m = numpy.einsum('ij,jk,ik->i', d, self.inv[i], d)
            out = numpy.logaddexp(out, -0.5 * m + self.lognorm[i]
                                  + self.logw[i])
        # change of variables from whitened z back to physical x
        return out + self._logdetW


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
                 gen_start_round=3,
                 gen_min_ess=200,
                 gen_fraction=0.3,
                 gen_max_centres=3000,
                 gen_refit=1,
                 gen_bandwidth=1.0,
                 gen_bandwidths=None,
                 gen_local_k=0,
                 gen_defensive=0.1,
                 gen_defensive_scale=3.0,
                 stretch_pairs=0,
                 stretch_scale=1.0,
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
        self.gen_defensive = float(gen_defensive)
        self.gen_defensive_scale = float(gen_defensive_scale)
        self.stretch_pairs = int(stretch_pairs)
        self.stretch_scale = float(stretch_scale)
        # seeded from the legacy global RNG, which pycbc_inference --seed sets,
        # so generative draws are reproducible alongside the shuffles that the
        # discrete-pool path uses
        self._rng = numpy.random.default_rng(
            numpy.random.randint(0, 2 ** 31 - 1))
        self.ncalls = 0
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

    def pool_round(self, bin_weight, node_idx, lengths, budget):
        """ One round of games3's discrete-pool draw. Returns
        (psamp, loglr, logw_unnormalised, bin_id) or raises OutOfSamples.
        """
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
        logw = loglr + numpy.log(lengths[bin_id]) - pweight
        return psamp, loglr, logw, bin_id

    # ---------------- generative stratum -----------------------------------

    def gen_round(self, proposal, budget, stratum='gen'):
        """ One round from the generative proposal. Out-of-prior draws are
        dropped before any likelihood call.
        """
        names = list(self.model.variable_params)
        x = proposal.sample(budget)
        keep = self._in_prior(x, names)
        if keep.sum() == 0:
            return None
        x = x[keep]
        logq = proposal.logpdf(x)

        psamp = FieldArray(len(x), dtype=self.dtype)
        for k, p in enumerate(names):
            psamp[p] = x[:, k]
        loglr = self._likelihoods(psamp, stratum=stratum)
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

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as mapfile:
            bparams = {p: mapfile['bank'][p][:] for p in self.variable_params}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            self.dtype = mapfile['map']['0'].dtype

        logging.info('Calculating likelihood at nodes')
        args = [{p: bparams[p][i] for p in self.model.variable_params}
                for i in range(num_nodes)]
        node_loglrs = numpy.array(list(tqdm.tqdm(
            self.pool.imap(call_likelihood, args), total=len(args))))
        self.meta['setup_ncalls'] = len(args)
        bound = node_loglrs[~numpy.isnan(node_loglrs)].max() - self.loglr_region
        passed = numpy.where((node_loglrs > bound)
                             & ~numpy.isnan(node_loglrs))[0]

        logging.info('...resolving tree leaves for %i passed tiles',
                     len(passed))
        active_lengths, key = [], 0
        treefile = h5py.File(self.treefile, 'r') if self.treefile else None
        with h5py.File(self.mapfile, 'r') as mapfile:
            for tile in passed:
                leaves = None
                if treefile is not None and str(tile) in treefile['tree']:
                    leaves = self._find_active_leaves(
                        treefile['tree'][str(tile)], bound)
                if not leaves:
                    leaves = [mapfile['map'][str(tile)][:]]
                for pts in leaves:
                    self.dmap[key] = pts
                    active_lengths.append(len(pts))
                    key += 1
        if treefile is not None:
            treefile.close()
        active_ids = numpy.arange(key)
        lengths = numpy.array(active_lengths)
        logging.info('...resolved to %i active leaves', key)

        weight = lengths / lengths.sum()
        # per-stratum accumulators
        gen_keys = [f'gen:bw{bw:g}' for bw in self.gen_bandwidths]
        strata = {'pool': {'samp': None, 'loglr': None, 'logw': None}}
        for k in gen_keys:
            strata[k] = {'samp': None, 'loglr': None, 'logw': None}
        self.stratum_calls = {k: 0 for k in ['pool'] + gen_keys}
        proposals = {}
        # per-stratum share of the generative budget, adapted by ESS per call
        share = {k: 1.0 / len(gen_keys) for k in gen_keys}
        pool_dead = False
        frac = self.gen_fraction

        for rnd in range(1, self.rounds + 1):
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
                        pool_budget)
                    before = self._ess(strata['pool']['logw'])
                    strata['pool'] = self._accumulate(strata['pool'],
                                                      ps, pl, pw)
                    after = self._ess(strata['pool']['logw'])
                    pool_gain = (after - before) / max(len(pw), 1)
                    # games3's adaptive leaf reweighting, unchanged
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
                    got = self.gen_round(proposals[k], b, stratum=k)
                    if got is None:
                        continue
                    gs, gl, gw = got
                    before = self._ess(strata[k]['logw'])
                    strata[k] = self._accumulate(strata[k], gs, gl, gw)
                    after = self._ess(strata[k]['logw'])
                    gen_gains[k] = (after - before) / max(len(gw), 1)
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
            gsum = sum(self._ess(strata[k]['logw']) for k in gen_keys)
            detail = ' '.join(
                f'{k.split(":")[1]}={self._ess(strata[k]["logw"]):.0f}'
                for k in gen_keys)
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

    def _fit_proposal(self, strata, bandwidth=None):
        """ Fit the generative proposal to the accumulated posterior samples
        from every stratum, resampled to equal weight.
        """
        names = list(self.model.variable_params)
        xs, ws = [], []
        for st in strata.values():
            if st['samp'] is None:
                continue
            xs.append(numpy.column_stack([st['samp'][p] for p in names]))
            lw = st['logw']
            ws.append(numpy.exp(lw - logsumexp(lw)))
        if not xs:
            return None
        x = numpy.vstack(xs)
        w = numpy.concatenate(ws)
        w = w / w.sum()
        ess = 1.0 / (w ** 2).sum()
        if ess < self.gen_min_ess:
            logging.debug('accumulated ESS %.1f below gen_min_ess %.1f; '
                          'not fitting generative proposal yet',
                          ess, self.gen_min_ess)
            return None
        # A properly WEIGHTED mixture over the samples: every sample is a
        # centre and carries its own posterior weight as its mixture weight.
        # Two earlier constructions were both wrong. Resampling by weight with
        # replacement duplicates a handful of points at low ESS (ESS 19 gives
        # ~19 distinct centres), making the covariance degenerate. Taking the
        # top-N by weight instead concentrates centres at the peak, which
        # under-disperses the proposal -- and under-dispersion is fatal for
        # importance sampling, because the target's tails then have almost no
        # proposal density and pi/q explodes there. Weighting the components
        # keeps every sample's contribution proportional to its posterior mass
        # while retaining full tail coverage.
        n = min(self.gen_max_centres, len(x))
        if n < len(x):
            order = self._rng.choice(len(x), size=n, replace=False)
        else:
            order = numpy.arange(len(x))
        centres = x[order]
        cw = w[order]
        if cw.sum() <= 0:
            return None
        # weighted covariance over ALL samples, which uses the full posterior
        # information rather than only the retained centres
        mu = (w[:, None] * x).sum(0)
        xc = x - mu
        cov = (w[:, None] * xc).T @ xc / max(1.0 - (w ** 2).sum(), 1e-12)
        if not numpy.all(numpy.isfinite(cov)):
            return None
        try:
            if self.gen_local_k > 0:
                return LocalCovarianceProposal(
                    centres, self._rng, k=self.gen_local_k,
                    scale=(self.gen_bandwidth if bandwidth is None
                           else bandwidth), weights=cw)
            prop = GaussianMixtureProposal(
                centres,
                cov * (self.gen_bandwidth if bandwidth is None
                       else bandwidth) ** 2, self._rng,
                stretch_pairs=self.stretch_pairs,
                stretch_scale=self.stretch_scale, weights=cw,
                defensive_frac=self.gen_defensive,
                defensive_scale=self.gen_defensive_scale)
        except numpy.linalg.LinAlgError:
            logging.info('generative proposal fit failed (singular kernel)')
            return None
        logging.info('fitted generative proposal: %i centres '
                     '(%i stretch), accumulated ESS %.1f',
                     prop.n, self.stretch_pairs, ess)
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
        logging.info('combined ESS=%.1f ncalls=%i (setup %i)', total_ess,
                     self.ncalls, self.meta.get('setup_ncalls', 0))

        self.meta['ncalls'] = self.ncalls
        self.meta['ess'] = total_ess
        self.meta['setup_ncalls'] = self.meta.get('setup_ncalls', 0)

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
        fsamp = numpy.concatenate(outs)
        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = numpy.concatenate(ols)

    def _find_active_leaves(self, tree_root, loglr_bound):
        results = []

        def descend(group):
            if group.attrs['is_leaf']:
                pts = numpy.zeros(int(group.attrs['n_points']),
                                  dtype=self.dtype)
                for p in self.variable_params:
                    pts[p] = group[f'points_{p}'][:]
                results.append(pts)
                return
            for c in group['children']:
                child = group['children'][c]
                params = {p: float(child.attrs[f'param_{p}'])
                          for p in self.variable_params}
                ll = call_likelihood(params)
                if numpy.isfinite(ll) and ll > loglr_bound:
                    descend(child)

        descend(tree_root)
        return results
