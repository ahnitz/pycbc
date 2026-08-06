""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space, with an optional hierarchical
tree (see precalc_mapping/build_tree.py) that recursively subdivides a
level-0 tile's own point cloud into smaller, more localized clusters.

This is a variant of the 'games' sampler (see games.py). It reuses that
sampler's core weighting/draw logic unchanged -- a tree leaf's points are
just as much genuine prior-conditional draws as a flat tile's are, so no
reweighting scheme change is needed, unlike the KDE-based 'games2'. What
changes is *which* points are treated as an active tile: instead of
always using a level-0 tile's full point cloud, this descends the tree
for each level-0 tile that passes the initial likelihood cut, evaluating
likelihood at each level's children to find smaller, better-localized
leaves. A level-0 tile that has no tree entry (not prebuilt), or whose
children/descendants all fail the likelihood cut, falls back to that
tile's full flat point cloud from `mapfile` -- so this never does worse
than 'games' for any given tile, only better where the tree helps.
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


class GameSampler3(DummySampler):
    """Direct importance sampling using a preconstructed parameter space
    mapping file, optionally refined by a hierarchical tree of smaller,
    more localized sub-tiles.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Path to the pre-generated file containing the pre-mapped prior
        volume (same format as used by the 'games' sampler).
    treefile : str, optional
        Path to a hierarchical tree file (see build_tree.py) covering
        some or all of mapfile's level-0 tiles. Tiles without a tree
        entry fall back to their flat mapfile pool.
    loglr_region: int
        Only use regions from the prior volume tiling that are within
        this loglr difference of the maximum tile, applied at every
        level of the tree descent.
    target_likelihood_calls: int
        Try to use this many likelihood calls in each round of the analysis.
    rounds: int
        The number of iterations to use before terminated.
    """
    name = 'games3'

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

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.loglr_region = float(loglr_region)

    # -- unchanged from games.py: a tree leaf's points are just as much
    # genuine prior-conditional draws as a flat tile's, so the same
    # finite-pool draw/weighting logic applies without modification.
    def draw_samples_from_bin(self, i, size):
        """ Get samples from the binned prior space """
        if i not in self.draw:
            self.draw[i] = numpy.arange(0, len(self.dmap[i]))

        if size > len(self.draw[i]):
            raise OutOfSamples

        numpy.random.shuffle(self.draw[i])
        selected = self.draw[i][:size]
        self.draw[i] = self.draw[i][size:]

        if size > 0:
            remain = len(self.draw[i])
            logging.info('Drew %i, %i remains in bin %i', size, remain, i)

        return self.dmap[i][selected]

    def sample_round(self, bin_weight, node_idx, lengths):
        """ Sample from the posterior using pre-binned sets of points and
        the weighting factor of each bin. Identical to games.py, including
        its exhaustion behavior: `lengths[]` is each bin's original, fixed
        pool size (not updated across rounds), so a bin depleted by earlier
        rounds raises OutOfSamples in draw_samples_from_bin and stops the
        run -- deliberately kept exactly as-is here (rather than
        redistributing a depleted, high-weight bin's shortfall into
        lower-weight bins, which was tried and found to dilute efficiency
        by spending likelihood calls in low-value regions) so this sampler
        isolates the tree hierarchy's own effect on efficiency, with
        exhaustion left as the same open problem it is for games.py.
        """
        logging.info("...draw template bins")
        drawcount = (bin_weight * self.target_likelihood_calls).astype(int)

        dorder = bin_weight.argsort()[::-1]
        remainder = 0
        for i in dorder:
            bincount = drawcount[i]
            binlen = lengths[i]
            if bincount > binlen:
                drawcount[i] = binlen
                remainder += bincount - binlen
            elif bincount < binlen:
                asize = min(binlen - bincount, remainder)
                drawcount[i] += asize
                remainder -= asize

        drawweight = bin_weight / drawcount
        total_draw = drawcount.sum()

        logging.info('...drawn random points within bins')
        psamp = FieldArray(total_draw, dtype=self.dtype)
        pweight = numpy.zeros(total_draw, dtype=float)
        bin_id = numpy.zeros(total_draw, dtype=int)
        j = 0
        for i, (c, w) in enumerate(zip(drawcount, drawweight)):
            bdraw = self.draw_samples_from_bin(node_idx[i], c)
            psamp[j:j+len(bdraw)] = FieldArray.from_records(bdraw,
                                                            dtype=self.dtype)
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
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

        logw3 = loglr_samp + numpy.log(lengths[bin_id]) - pweight
        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)
        return psamp, loglr_samp, weight2, bin_id

    def _find_active_leaves(self, tree_root, loglr_bound):
        """ Recursively descend the tree for one level-0 tile, returning
        a list of points_structured_array, one per leaf reached whose own
        likelihood exceeds loglr_bound. The root's own likelihood is
        assumed already verified by the caller (it's a member of
        `passed` from the level-0 evaluation); only its descendants are
        checked here. Returns an empty list if the root has children but
        none of them pass (the caller falls back to the flat mapfile
        pool in that case).
        """
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
                    leaves = [mapfile['map'][str(tile)][:]]
                for pts in leaves:
                    self.dmap[key] = pts
                    active_lengths.append(len(pts))
                    key += 1
        if treefile is not None:
            treefile.close()

        active_ids = numpy.arange(key)
        lengths = numpy.array(active_lengths)
        logging.info("...resolved to %i active leaves (avg %.0f points each)",
                    key, lengths.mean())

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
