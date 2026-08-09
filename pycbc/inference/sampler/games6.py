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
from .games import call_likelihood, OutOfSamples


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
    gen_local_k : int
        Neighbours used to estimate each centre's own kernel covariance.
    """
    name = 'games6'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None, treefile=None,
                 loglr_region=0,
                 loglr_coverage=1e-4,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 gen_start_round=1,
                 gen_min_ess=100,
                 gen_max_centres=20000,
                 gen_refit=1,
                 gen_bandwidth=1.0,
                 gen_local_k=160,
                 tree_start_level=0,
                 min_active_points=100,
                 tile_point_cap=200000,
                 gen_seed_leaves=0,
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
        self.tree_start_level = int(tree_start_level)
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
        self.gen_seed_range = float(gen_seed_range)
        self.leaf_loglr = {}
        # seeded from the global RNG that pycbc_inference --seed sets, so
        # generative draws are reproducible alongside the pool shuffles
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

        psamp = FieldArray(len(x), dtype=self.dtype)
        for k, p in enumerate(names):
            psamp[p] = x[:, k]
        loglr = self._likelihoods(psamp, stratum='gen')
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

        # Where to start evaluating. The cut compares the likelihood at
        # each start node's representative, so it is only as good as that
        # representative proxies the likelihood across its whole tile. A
        # coarser representative is a worse proxy, which is why the depth
        # the cut starts from is separate from the depths the tree was
        # built with.
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

        # Can the cut discriminate at all? The bound is
        # max - loglr_region, so every node passes exactly when
        #
        #     max(node) - min(node) < loglr_region
        #
        # written as a spread rather than an absolute level because
        # `call_likelihood` returns the log likelihood, not the loglr, so
        # only differences are offset-free. When nothing can be pruned,
        # descending would enumerate the whole tree and still discard
        # leaves on representative-level noise, so use the tiles
        # directly: their subtree clouds union to the prior.
        finite = node_loglrs[numpy.isfinite(node_loglrs)]
        spread = (finite.max() - finite.min()) if len(finite) else 0.0
        coarse = start_nodes is not None and spread < self.loglr_region
        if coarse:
            logging.info('%i of %i start nodes passed: cut cannot '
                         'discriminate, using tiles directly',
                         len(passed), len(args))

        if coarse:
            active_ids, lengths, prior_mass = self._tiles_as_bins(
                passed, start_nodes, node_loglrs)
            key = len(active_ids)
        else:
            logging.info('...resolving tree leaves for %i passed tiles',
                         len(passed))
            active_lengths, active_mass, key = [], [], 0
            self.leaf_loglr = {}
            with h5py.File(self.mapfile, 'r') as mapfile:
                for tile in passed:
                    leaves = None
                    if start_nodes is not None:
                        leaves = self._find_active_leaves(start_nodes[tile],
                                                          bound)
                    elif treefile is not None and str(tile) in treefile['tree']:
                        leaves = self._find_active_leaves(
                            treefile['tree'][str(tile)], bound)
                    if not leaves and start_nodes is None:
                        mp = mapfile['map'][str(tile)][:]
                        leaves = [(mp, node_loglrs[tile], float(len(mp)))]
                    for pts, ll, mass in (leaves or []):
                        self.dmap[key] = pts
                        self.leaf_loglr[key] = ll
                        active_lengths.append(len(pts))
                        active_mass.append(mass)
                        key += 1
            active_ids = numpy.arange(key)
            lengths = numpy.array(active_lengths)
            prior_mass = numpy.array(active_mass, dtype=float)
        # NB: the tree file must stay open until after the fallback
        # below, which reads subtrees through `start_nodes`.
        logging.info('...resolved to %i active bins (%i points)',
                     key, int(lengths.sum()) if key else 0)

        # The cut can also starve the pool rather than fail to prune, if
        # very few tiles pass and they hold very few points. Widen until
        # usable: free, since every start node was already evaluated, and
        # unbiased, since the extra tiles are low-likelihood and weighted
        # accordingly.
        if start_nodes is not None and not coarse \
                and (not key or lengths.sum() < self.min_active_points):
            logging.warning('active set starved (%i leaves, %i points): '
                            'falling back to start-level tiles',
                            key, int(lengths.sum()) if key else 0)
            region = self.loglr_region
            for _ in range(6):
                sel = numpy.where(
                    (node_loglrs > numpy.nanmax(node_loglrs) - region)
                    & ~numpy.isnan(node_loglrs))[0]
                active_ids, lengths, prior_mass = self._tiles_as_bins(
                    sel, start_nodes, node_loglrs)
                key = len(active_ids)
                if lengths.sum() >= self.min_active_points:
                    break
                region *= 2.0
                logging.warning('still starved (%i points); widening '
                                'loglr_region to %.0f', lengths.sum(), region)
            logging.info('...fallback gave %i tiles (%i points)',
                         key, int(lengths.sum()) if key else 0)

        if treefile is not None:
            treefile.close()

        weight = self._round1_weights(key, lengths)
        # per-stratum accumulators
        strata = {'pool': {'samp': None, 'loglr': None, 'logw': None}}
        self.stratum_calls = {'pool': 0}
        proposal = self._seed_proposal_from_leaves() \
            if self.gen_seed_leaves else None
        pool_dead = False

        for rnd in range(1, self.rounds + 1):
            self._round = rnd
            # The generative proposal is held in RESERVE: the pool keeps the
            # entire budget, so behaviour before exhaustion is identical to
            # games3, and the generative stratum only runs once the pool dies
            # and hands its budget over.
            gen_budget = self.target_likelihood_calls if pool_dead else 0
            pool_budget = self.target_likelihood_calls - gen_budget

            if pool_budget > 0 and not pool_dead:
                try:
                    ps, pl, pw, bid = self.pool_round(
                        weight / weight.sum(), active_ids, lengths,
                        pool_budget, prior_mass=prior_mass)
                    strata['pool'] = self._accumulate(strata['pool'],
                                                      ps, pl, pw)
                    # concentrate later rounds on the leaves that are
                    # actually carrying posterior weight
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
                    gen_budget = self.target_likelihood_calls

            if gen_budget > 0 and proposal is not None:
                got = self.gen_round(proposal, gen_budget)
                if got is not None:
                    gs, gl, gw = got
                    # one stratum PER ROUND: each round used a different
                    # refitted proposal, so its samples are only commensurate
                    # with themselves
                    rk = 'gen:r%i' % rnd
                    strata[rk] = {'samp': gs, 'loglr': gl, 'logw': gw}
                    self.stratum_calls[rk] = len(gw)

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
            logging.info('round %i: ESS pool=%.1f gen[%i rounds]=%.1f '
                         'total=%.1f ncalls=%i', rnd, ess_p, len(gkeys),
                         gsum, ess_p + gsum, self.ncalls)

        self._finalise(strata)

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
        ok = numpy.isfinite(ll)
        if not ok.any():
            return volume
        spread = float(ll[ok].max() - ll[ok].min())
        temper = max(1.0, spread / max(self.gen_seed_range, 1e-6))
        w = lengths.astype(float).copy()
        w[ok] *= numpy.exp((ll[ok] - ll[ok].max()) / temper)
        if w.sum() <= 0:
            return volume
        logging.info('seeded round-1 leaf weights: temper %.2f over a '
                     'loglr spread of %.1f nats in %i leaves',
                     temper, spread, int(ok.sum()))
        return w / w.sum()

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
        fa = FieldArray(len(samp), dtype=self.dtype)
        for p in names:
            fa[p] = samp[p]
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
        try:
            prop = LocalCovarianceProposal(
                centres, self._rng, k=self.gen_local_k,
                scale=self.gen_bandwidth, weights=cw)
        except numpy.linalg.LinAlgError:
            logging.info('generative proposal fit failed (singular kernel)')
            return None
        logging.info('fitted proposal: k=%i bw=%.2f centres=%i '
                     'accumulated ESS=%.1f', self.gen_local_k,
                     self.gen_bandwidth, prop.n, ess)
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

    def _find_active_leaves(self, tree_root, loglr_bound):
        """ Prune the subtree, evaluating the likelihood at each child's
        representative and descending only where it clears the bound.

        These are real likelihood calls and are counted into
        `self.tree_ncalls`, which is reported separately from `ncalls`.
        """
        results = []

        def descend(group, ll_here):
            if group.attrs['is_leaf']:
                pts = numpy.zeros(int(group.attrs['n_points']),
                                  dtype=self.dtype)
                for p in self.variable_params:
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
