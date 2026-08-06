""" Direct monte carlo sampling using a hierarchical tree of tiles (see
precalc_mapping/build_tree.py), where validated leaves are sampled from
an analytic ellipsoid (Gaussian in the leaf's metric-informative
subspace) instead of a finite discrete pool.

This is a variant of 'games3' (see games3.py), which itself reuses
'games' (games.py) unchanged for discrete-pool leaves -- a tree leaf's
own points are just as much genuine prior-conditional draws as a flat
tile's, so no reweighting change is needed there. What's different here
is *validated* leaves: instead of always popping from that leaf's finite
point pool (which can run dry under heavy resampling, exactly the
scenario a validated leaf is most useful in -- high SNR concentrating
weight onto a handful of small, well-localized regions), draw fresh from
an EllipsoidLeaf (metric_utils.EllipsoidLeaf): an analytic Gaussian in
the eigendirections the metric fit actually constrains, combined with
bootstrapped (real, existing-point) values in the directions the metric
found genuinely unconstrained -- not just poorly resolved, confirmed
during tree building, so forcing a Gaussian shape there would be wrong,
not just imprecise (this is what made the earlier full-rank KDE
approach, 'games2', fail).

An ellipsoid draw's proposal density differs from the implicit,
no-correction-needed density of a discrete-pool draw (which is prior-
correct by construction, since every pool member already *is* a genuine
prior sample), so it needs an explicit importance-sampling correction:
weight get an extra `-logq_informative(x)` term (the true prior's own
density is a global constant here, since every parameter's prior is
uniform, so it cancels in the final self-normalization and never needs
to be evaluated). Leaves without a valid metric fall back to the
discrete-pool draw exactly as games3 does, and tiles/leaves without a
usable tree entry at all fall back further to the flat mapfile pool --
so this never does worse than games3 for any given tile, only better
where a validated ellipsoid is available.
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

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '..', '..', '..', '..', 'precalc_mapping'))
from metric_utils import EllipsoidLeaf  # noqa: E402


class GameSampler4(DummySampler):
    """Direct importance sampling using a hierarchical tree of tiles,
    with validated leaves sampled from an analytic ellipsoid instead of
    a finite discrete pool.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Path to the pre-generated flat mapping file, used for level-0
        node parameters and as the fallback pool for tiles/leaves with
        no usable tree entry.
    treefile : str
        Path to a hierarchical tree file (see build_tree.py).
    loglr_region: int
        Only use regions from the prior volume tiling that are within
        this loglr difference of the maximum tile, applied at every
        level of the tree descent.
    target_likelihood_calls: int
        Try to use this many likelihood calls in each round of the analysis.
    rounds: int
        The number of iterations to use before terminated.
    """
    name = 'games4'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None, treefile=None,
                 loglr_region=25,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 **kwargs):
        super().__init__(model, *args)

        self.meta = {}
        self.mapfile = mapfile
        self.treefile = treefile
        self.rounds = int(rounds)
        self.dmap = {}
        self.draw = {}
        self.ellipsoid_draws = {}
        # Derived from the legacy global numpy RNG (seeded by
        # pycbc_inference's --seed via numpy.random.seed) so ellipsoid
        # draws are reproducible run-to-run, same as the discrete-pool
        # draws elsewhere in this sampler which use numpy.random.* directly.
        self._rng = numpy.random.default_rng(
            numpy.random.randint(0, 2**31 - 1))

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.loglr_region = float(loglr_region)

    def _in_bounds(self, pt_dict):
        """ True prior support check, looped per-point since some prior
        distributions (e.g. uniform's boundary check) aren't vectorized.
        """
        names = list(pt_dict.keys())
        n = len(pt_dict[names[0]])
        ok = numpy.zeros(n, dtype=bool)
        for k in range(n):
            pset = {p: pt_dict[p][k] for p in names}
            ok[k] = numpy.isfinite(self.model.prior_distribution(**pset))
        return ok

    def draw_samples_from_bin(self, i, size):
        """ Get samples from bin i: an analytic ellipsoid draw if this
        leaf validated during tree building, otherwise the same finite-
        pool draw games.py/games3.py use. Always returns (FieldArray of
        length size, logq array of length size); logq is zero for
        discrete-pool draws (no explicit correction needed there).
        """
        leaf = self.dmap[i]
        if isinstance(leaf, EllipsoidLeaf):
            pt_dict, logq = leaf.sample(size, self._in_bounds)
            self.ellipsoid_draws[i] = self.ellipsoid_draws.get(i, 0) + size
            out = FieldArray(size, dtype=self.dtype)
            for p in self.dtype.names:
                out[p] = pt_dict[p]
            return out, logq

        if i not in self.draw:
            self.draw[i] = numpy.arange(0, len(leaf))

        if size > len(self.draw[i]):
            raise OutOfSamples

        numpy.random.shuffle(self.draw[i])
        selected = self.draw[i][:size]
        self.draw[i] = self.draw[i][size:]

        if size > 0:
            remain = len(self.draw[i])
            logging.info('Drew %i, %i remains in bin %i', size, remain, i)

        return leaf[selected], numpy.zeros(size)

    # An EllipsoidLeaf never runs dry -- its informative-subspace draws
    # are fresh Gaussian samples, and its null-subspace draws bootstrap
    # (with replacement) from the leaf's own finite validation pool.
    # That bootstrap pool's *own* empirical distribution has finite-
    # sample error relative to the true local distribution, and that
    # error does not shrink just because we keep resampling it --
    # letting one small, possibly-unrepresentative pool get resampled
    # thousands of times (observed: 1400+ draws/round from a single
    # ~900-point leaf) lets whatever quirk that pool has dominate the
    # posterior far beyond what its actual sample size supports. Capping
    # total draws at a healthy multiple of the pool size keeps the
    # "never exhausts after one round" fix while bounding how much any
    # one small pool can end up dominating the run.
    MAX_ELLIPSOID_OVERSAMPLE = 20

    def _remaining_capacity(self, key, total):
        """ How many more points bin `key` can still supply *this round*
        (see games3.GameSampler3._remaining_capacity for why `lengths[]`
        alone -- the bin's original, fixed pool size -- is not safe to
        reuse across rounds).
        """
        leaf = self.dmap[key]
        if isinstance(leaf, EllipsoidLeaf):
            cap = self.MAX_ELLIPSOID_OVERSAMPLE * total
            used = self.ellipsoid_draws.get(key, 0)
            return max(0, cap - used)
        if key in self.draw:
            return len(self.draw[key])
        return total

    def sample_round(self, bin_weight, node_idx, lengths):
        """ Sample from the posterior using pre-binned sets of points and
        the weighting factor of each bin. Same as games3.py, plus a
        `logq_samp` correction term subtracted for any ellipsoid draws.
        """
        logging.info("...draw template bins")
        drawcount = (bin_weight * self.target_likelihood_calls).astype(int)

        dorder = bin_weight.argsort()[::-1]
        remainder = 0
        for i in dorder:
            bincount = drawcount[i]
            binlen = self._remaining_capacity(node_idx[i], lengths[i])
            if bincount > binlen:
                drawcount[i] = binlen
                remainder += bincount - binlen
            elif bincount < binlen:
                asize = min(binlen - bincount, remainder)
                drawcount[i] += asize
                remainder -= asize

        with numpy.errstate(divide='ignore', invalid='ignore'):
            drawweight = bin_weight / drawcount
        total_draw = drawcount.sum()

        logging.info('...drawn random points within bins')
        psamp = FieldArray(total_draw, dtype=self.dtype)
        pweight = numpy.zeros(total_draw, dtype=float)
        logq_samp = numpy.zeros(total_draw, dtype=float)
        bin_id = numpy.zeros(total_draw, dtype=int)
        j = 0
        for i, (c, w) in enumerate(zip(drawcount, drawweight)):
            bdraw, logq = self.draw_samples_from_bin(node_idx[i], c)
            psamp[j:j+len(bdraw)] = bdraw
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
            logq_samp[j:j+len(bdraw)] = logq
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        logging.info("Possible unique values %s", lengths.sum())
        logging.info("Templates drawn from %s", len(lengths))
        logging.info("Unique values first draw %s", len(psamp))

        args = []
        for i in range(len(psamp)):
            pset = {p: psamp[p][i] for p in self.model.variable_params}
            args.append(pset)

        loglr_samp = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                    total=len(args)))
        loglr_samp = numpy.array(loglr_samp)

        # true prior density is a global constant here (every parameter's
        # prior is uniform), so it cancels in the self-normalization below
        # and never needs to be evaluated. +log(lengths[bin]) is the
        # Horvitz-Thompson correction for a finite pool of lengths[bin]
        # real prior draws; it applies to discrete-pool bins AND to an
        # ellipsoid's null-subspace coordinates (bootstrapped from that
        # same finite pool of lengths[bin] points), so it's kept
        # unconditionally here. Only the informative-subspace Gaussian
        # draw gets an *additional*, explicit density correction, since
        # that part isn't drawn from the finite pool at all.
        logw3 = (loglr_samp + numpy.log(lengths[bin_id]) - pweight
                - logq_samp)

        import os
        if os.environ.get('GAMES4_DUMP'):
            dump_kwargs = {p: psamp[p] for p in self.model.variable_params}
            dump_kwargs.update(
                loglr_samp=loglr_samp, lengths_bin=lengths[bin_id],
                pweight=pweight, logq_samp=logq_samp, bin_id=bin_id,
                is_ellipsoid=numpy.array(
                    [isinstance(self.dmap[node_idx[b]], EllipsoidLeaf)
                     for b in bin_id]),
                weight2=numpy.exp(logw3 - logsumexp(logw3)))
            numpy.savez(os.environ['GAMES4_DUMP'], **dump_kwargs)
            logging.info('dumped round diagnostics to %s',
                         os.environ['GAMES4_DUMP'])

        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)
        return psamp, loglr_samp, weight2, bin_id

    def _find_active_leaves(self, tree_root, loglr_bound):
        """ Recursively descend the tree for one level-0 tile, returning
        a list of (is_valid, metric_or_None, node_params, points) for
        every leaf reached whose own likelihood exceeds loglr_bound.
        """
        results = []

        def descend(group):
            if group.attrs['is_leaf']:
                pts = numpy.zeros(int(group.attrs['n_points']),
                                  dtype=self.dtype)
                for p in self.variable_params:
                    pts[p] = group[f'points_{p}'][:]
                node_params = {p: float(group.attrs[f'param_{p}'])
                              for p in self.variable_params}
                is_valid = bool(group.attrs.get('valid', False))
                metric = group['metric'][:] if is_valid else None
                results.append((is_valid, metric, node_params, pts))
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

    def run(self):
        """ Produce posterior samples """
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as mapfile:
            bparams = {p: mapfile['bank'][p][:] for p in self.variable_params}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            self.dtype = mapfile['map']['0'].dtype

        logging.info('Calculating likelihood at nodes')
        args = []
        for i in range(num_nodes):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            args.append(pset)

        node_loglrs = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args)))
        node_loglrs = numpy.array(node_loglrs)
        loglr_bound = node_loglrs[~numpy.isnan(node_loglrs)].max()
        loglr_bound -= self.loglr_region

        passed = (node_loglrs > loglr_bound) & ~numpy.isnan(node_loglrs)
        passed = numpy.where(passed)[0]

        logging.info("...resolving tree leaves for %i passed tiles",
                    len(passed))
        active_lengths = []
        n_ellipsoid = 0
        key = 0
        treefile = h5py.File(self.treefile, 'r') if self.treefile else None
        with h5py.File(self.mapfile, 'r') as mapfile:
            for tile in passed:
                leaves = None
                if treefile is not None and str(tile) in treefile['tree']:
                    leaves = self._find_active_leaves(
                        treefile['tree'][str(tile)], loglr_bound)
                if not leaves:
                    logging.info('tile %i: using flat fallback pool', tile)
                    pts = mapfile['map'][str(tile)][:]
                    self.dmap[key] = pts
                    active_lengths.append(len(pts))
                    key += 1
                    continue
                for is_valid, metric, node_params, pts in leaves:
                    if is_valid:
                        pts_dict = {p: pts[p] for p in self.variable_params}
                        self.dmap[key] = EllipsoidLeaf(
                            node_params, list(self.variable_params), metric,
                            pts_dict, self._rng)
                        n_ellipsoid += 1
                    else:
                        self.dmap[key] = pts
                    active_lengths.append(len(pts))
                    key += 1
        if treefile is not None:
            treefile.close()

        active_ids = numpy.arange(key)
        lengths = numpy.array(active_lengths)
        logging.info("...resolved to %i active leaves (%i ellipsoid, avg "
                    "%.0f points each)", key, n_ellipsoid, lengths.mean())

        # Sample from posterior
        psamp = None
        loglr_samp = None
        weight2 = None
        bin_ids = None

        weight = lengths / lengths.sum()

        for i in range(self.rounds):
            try:
                psamp_v, loglr_samp_v, weight2_v, bin_id = \
                    self.sample_round(weight / weight.sum(),
                                      active_ids, lengths)
            except OutOfSamples:
                logging.info("No more samples to draw from")
                break

            for j, v in zip(bin_id, weight2_v):
                weight[j] += v

            if psamp is None:
                psamp = psamp_v
                loglr_samp = loglr_samp_v
                weight2 = weight2_v
                bin_ids = bin_id
            else:
                psamp = numpy.concatenate([psamp_v, psamp])
                loglr_samp = numpy.concatenate([loglr_samp_v, loglr_samp])
                weight2 = numpy.concatenate([weight2_v, weight2])
                bin_ids = numpy.concatenate([bin_id, bin_ids])

            ess = 1.0 / ((weight2/weight2.sum()) ** 2.0).sum()
            logging.info("ESS = %s", ess)

        self.meta['ncalls'] = len(weight2)
        self.meta['ess'] = ess

        weight2 /= weight2.sum()
        draw2 = numpy.random.choice(len(psamp), size=int(ess * 1),
                                    replace=True, p=weight2)
        logging.info("Unique values second draw %s",
                     len(numpy.unique(psamp[draw2])))

        fsamp = FieldArray(len(draw2), dtype=self.dtype)
        for i, v in enumerate(draw2):
            fsamp[i] = psamp[v]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = loglr_samp[draw2]
