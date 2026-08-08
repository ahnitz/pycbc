""" games6 reading a SHAPE-COMPRESSED tree, with the points regenerated.

games6 is deliberately untouched: it is the known-best baseline and this
module must not be able to perturb it.

What this changes
-----------------
A leaf is ``{x : match(x, rep) > t}`` and match is locally quadratic, so a
leaf's prior draws are approximately uniform on a metric ellipsoid --
measured radial kurtosis 3.96, against 4.47 for uniform-on-ellipsoid and
3.05 for a Gaussian. That is also why denser maps helped while deeper
trees did not: the stored points are a Monte Carlo rendering of a shape
that can be written down. So the map keeps a shape per leaf instead of a
point list, and the points are REGENERATED at load. 845 MB -> 11 MB, and
per-leaf density stops being fixed by what was written to disk.

Why regenerate rather than propose from the shapes
--------------------------------------------------
An earlier attempt used the ellipsoid mixture directly as a proposal and
failed badly: 0.22% efficiency against games6's 22%. The reason is worth
recording. A radius covering every point is an extremal statistic --
median 5.15 Mahalanobis units, so each ellipsoid is ~1.9e4 times the unit
volume, and uniform sampling in 6-D puts 88% of draws beyond 0.7 r_max,
in a shell holding essentially no real leaf points. Shrinking the radius
to fix that then leaves the shape covering only the bulk, and pairing a
bulk proposal with a shell-only residual pool breaks the stratum
combination: games6 combines strata with beta proportional to ESS, which
assumes each stratum independently estimates the FULL posterior, not
disjoint pieces of it.

Regenerating sidesteps all of it. The shapes only have to be a good
COMPRESSION of the point cloud, not a good proposal, and everything
downstream is the validated games6 machinery. Because the regenerated
count is proportional to the leaf's ORIGINAL count, ``lengths`` keeps its
games6 meaning -- prior volume and draw capacity simultaneously -- so
every weight formula is untouched.

``regen_factor`` > 1 exploits the one thing a shape has over a point
list: no fixed sample budget.
"""
import logging

import h5py
import numpy
from scipy.special import logsumexp

from .games6 import GameSampler6
from .games import call_likelihood, OutOfSamples


class GameSampler7(GameSampler6):
    """ games6 over a shape-compressed tree.

    Parameters as games6, plus:

    regen_factor : float
        Points regenerated per leaf as a multiple of its original count.
        1.0 reproduces the original density.
    """
    name = 'games7'

    def __init__(self, model, *args, regen_factor=1.0, **kwargs):
        super().__init__(model, *args, **kwargs)
        self.regen_factor = float(regen_factor)

    def _load_flat(self):
        f = h5py.File(self.treefile, 'r')
        node = {k: f['node/%s' % k][:] for k in
                ('param', 'is_leaf', 'child_start', 'child_count',
                 'child_index', 'leaf_id', 'roots')}
        sh = {k: f['shape/%s' % k][:] for k in ('mu', 'cov_tril', 'radius')}
        leaf = {'shape_id': f['leaf/shape_id'][:],
                'n_total': f['leaf/n_total'][:]}
        f.close()
        return node, sh, leaf

    def _regenerate(self, mu, cov_tril, radius, n):
        """ n draws uniform on the ellipsoid (mu, cov, radius). """
        d = len(self.variable_params)
        c = numpy.zeros((d, d))
        c[numpy.tril_indices(d)] = cov_tril
        c = c + c.T - numpy.diag(numpy.diag(c))
        eps = 1e-12 * max(abs(numpy.trace(c)) / d, 1e-300)
        for _ in range(8):
            try:
                L = numpy.linalg.cholesky(c)
                break
            except numpy.linalg.LinAlgError:
                c = c + eps * numpy.eye(d)
                eps *= 10
        else:
            return None
        # Scale by sqrt(d+2), NOT by the stored radius.
        #
        # A uniform distribution inside a d-dimensional ellipsoid of
        # Mahalanobis radius R has covariance R^2/(d+2) * Sigma -- it is
        # NOT Sigma. Using the stored radius (median 4.44) therefore
        # over-disperses every regenerated cloud by R/sqrt(d+2) = 1.57,
        # measured identically in all six parameters, which is also what
        # pushed regenerated points outside the prior box.
        #
        # sqrt(d+2) is the radius at which a uniform ellipsoid reproduces
        # the point covariance exactly. The cost is that the cloud then
        # stops at 2.83 Mahalanobis units while the real leaf reaches
        # ~4.44: a uniform ellipsoid cannot match both the covariance and
        # the extent, because the real leaf is not uniform. Matching the
        # covariance is the right choice -- it fixes the bulk, which is
        # what the pool weights depend on.
        z = self._rng.normal(size=(n, d))
        z /= numpy.maximum(numpy.linalg.norm(z, axis=1)[:, None], 1e-300)
        z *= (self._rng.random(n) ** (1.0 / d))[:, None] * numpy.sqrt(d + 2.0)
        return mu + z @ L.T

    def run(self):
        node, sh, leaf = self._load_flat()
        logging.info('shape tree: %i nodes, %i leaves, %i shapes',
                     len(node['param']), int((node['is_leaf'] > 0).sum()),
                     len(sh['mu']))

        def evaluate(ids):
            args = [{p: float(node['param'][i][k])
                     for k, p in enumerate(self.variable_params)}
                    for i in ids]
            return numpy.array([call_likelihood(a) for a in args])

        # Start the cut at tree_start_level, exactly as games6 does.
        # Descending from the roots evaluates the cut against 0.7-match
        # representatives, which are ~46 nats below the peak at SNR 13.5
        # and so barely discriminate: measured, that left 2577 active
        # leaves where games6 at level 2 keeps 139, an 18x larger active
        # region, and the KDE fitted over it started at 0.3% instead of
        # 16.45% efficiency.
        frontier = list(node['roots'])
        for _ in range(max(self.tree_start_level, 0)):
            nxt = []
            for i in frontier:
                if node['is_leaf'][i]:
                    nxt.append(int(i))
                    continue
                st = node['child_start'][i]
                nxt.extend(int(c) for c in node['child_index'][
                    st:st + node['child_count'][i]])
            if not nxt:
                break
            frontier = nxt
        logging.info('cut starts at tree depth %i: %i nodes',
                     self.tree_start_level, len(frontier))
        ll = evaluate(frontier)
        self.tree_ncalls += len(frontier)
        self.meta['setup_ncalls'] = len(frontier)
        bound = numpy.nanmax(ll) - self.loglr_region
        active, active_ll = [], []
        while frontier:
            nxt = []
            for i, l in zip(frontier, ll):
                if not numpy.isfinite(l) or l <= bound:
                    continue
                if node['is_leaf'][i]:
                    active.append(int(i))
                    active_ll.append(float(l))
                else:
                    s = node['child_start'][i]
                    nxt.extend(int(c) for c in node['child_index'][
                        s:s + node['child_count'][i]])
            if not nxt:
                break
            ll = evaluate(nxt)
            self.tree_ncalls += len(nxt)
            frontier = nxt
        if not active:
            raise RuntimeError('no active leaves; loglr_region too tight')

        self.dtype = [(p, numpy.float32) for p in self.variable_params]
        bounds = self._prior_bounds()
        if bounds is None:
            logging.warning('no finite prior box; regenerated points cannot '
                            'be filtered and may fall outside the prior')
        lengths, key, regen = [], 0, 0
        self.leaf_loglr = {}
        for i, l in zip(active, active_ll):
            lid = int(node['leaf_id'][i])
            sid = int(leaf['shape_id'][lid])
            if sid < 0:
                continue
            n = max(int(round(int(leaf['n_total'][lid])
                              * self.regen_factor)), 1)
            # Oversample and REJECT outside the prior. The fitted
            # ellipsoid routinely bulges past the prior box -- e.g. a leaf
            # with lambda1 mean 103.5, sd 88 and radius 4.44 reaches
            # lambda1 ~ -287 -- whereas every original leaf point was a
            # valid prior draw. games6 never hits this because its
            # generative stratum drops out-of-prior draws before the
            # likelihood call; regenerated points go into the POOL, which
            # has no such check, so they must be filtered here or they
            # poison the weights (measured: round-1 ESS 1.1).
            x = None
            for _ in range(6):
                cand = self._regenerate(sh['mu'][sid].astype(float),
                                        sh['cov_tril'][sid].astype(float),
                                        float(sh['radius'][sid]),
                                        int(n * 3))
                if cand is None:
                    break
                inb = numpy.ones(len(cand), dtype=bool)
                if bounds is not None:
                    lo, hi = bounds
                    inb = numpy.all((cand >= lo) & (cand <= hi), axis=1)
                cand = cand[inb]
                x = cand if x is None else numpy.vstack([x, cand])
                if len(x) >= n:
                    x = x[:n]
                    break
            if x is None or len(x) == 0:
                continue
            arr = numpy.zeros(len(x), dtype=self.dtype)
            for j, p in enumerate(self.variable_params):
                arr[p] = x[:, j]
            self.dmap[key] = arr
            self.leaf_loglr[key] = l
            lengths.append(len(arr))
            regen += len(arr)
            key += 1
        lengths = numpy.array(lengths)
        active_ids = numpy.arange(key)
        logging.info('active leaves %i, regenerated %i points '
                     '(regen_factor %.2f), tree cut %i calls',
                     key, regen, self.regen_factor, self.tree_ncalls)

        # ---- KDE route ----
        #
        # The regenerated points feed the GENERATIVE stratum, not the
        # discrete pool. Routing them to the pool cannot work: pool_round
        # treats n_bin as the prior volume of a DISJOINT region holding
        # uniform draws, but leaves are tightly-packed Voronoi cells
        # (median sibling separation 0.59 whitened units) while a
        # covariance-matched ellipsoid has radius sqrt(d+2) = 2.83, so
        # regenerated clouds overlap ~5x and the volume bookkeeping is
        # meaningless. Measured: 0.24% efficiency against games6's 22%.
        #
        # A KDE needs no partition. It is fitted to weighted samples and
        # proposes with its own analytic density, so overlapping source
        # clouds are irrelevant -- only that the cloud reproduces the
        # original distribution, which it now does (sd ratio 1.000 after
        # the sqrt(d+2) fix).
        strata = {'seed': {'samp': None, 'loglr': None, 'logw': None}}
        gen_keys = ['gen:bw%g' % bw for bw in self.gen_bandwidths]
        for k in gen_keys:
            strata[k] = {'samp': None, 'loglr': None, 'logw': None}
        self.stratum_calls = {k: 0 for k in gen_keys}

        prop = self._seed_proposal_from_leaves()
        if prop is None:
            raise RuntimeError('could not fit a proposal from the '
                               'regenerated leaf points')
        proposals = {gen_keys[0]: prop}

        for rnd in range(1, self.rounds + 1):
            self._round = rnd
            budget = self.target_likelihood_calls
            for k in list(proposals):
                got = self.gen_round(proposals[k], budget, stratum=k)
                if got is None:
                    continue
                gs, gl, gw = got
                rk = '%s:r%i' % (k, rnd)
                strata[rk] = {'samp': gs, 'loglr': gl, 'logw': gw}
                self.stratum_calls[rk] = len(gw)
            if self.gen_refit > 0 and rnd % self.gen_refit == 0:
                new = self._fit_proposals(strata, gen_keys)
                if new:
                    proposals = new
            tot = sum(self._ess(v['logw']) for v in strata.values()
                      if v['samp'] is not None)
            logging.info('round %i: total ESS=%.1f ncalls=%i', rnd, tot,
                         self.ncalls)

        self._finalise(strata)
