""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space.

This is a variant of the 'games' sampler (see games.py) that replaces the
finite pregenerated per-tile sample pool with a KDE fit to that same pool.
The pregenerated samples for a tile are themselves draws from the true
prior restricted to that tile's region (that's how the mapfile is built),
so a KDE fit to them approximates the tile's conditional prior density.
Sampling from the KDE instead of popping from the finite pool lets us draw
an unbounded number of fresh proposal points per tile, which both avoids
running out of samples for heavily-weighted tiles and removes the need to
ship a mapfile large enough to cover the full prior volume up front.
"""
import logging
import time
import tqdm
import h5py
import numpy
import numpy.random
from scipy.special import logsumexp
from scipy.stats import gaussian_kde
from pycbc.io import FieldArray

from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler


def call_likelihood(params):
    """ Accessor to update the global model
    """
    models._global_instance.update(**params)
    return models._global_instance.loglikelihood


class GameSampler2(DummySampler):
    """Direct importance sampling using a preconstructed parameter space
    mapping file, with a KDE fit to each tile's pregenerated samples used
    as an inexhaustible proposal in place of the finite sample pool.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Path to the pre-generated file containing the pre-mapped prior volume
    loglr_region: int
        Only use regions from the prior volume tiling that are within
        this loglr difference of the maximum tile.
    target_likelihood_calls: int
        Try to use this many likelihood calls in each round of the analysis.
    rounds: int
        The number of iterations to use before terminated.
    """
    name = 'games2'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None,
                 loglr_region=25,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 **kwargs):
        super().__init__(model, *args)

        self.meta = {}
        self.mapfile = mapfile
        self.rounds = int(rounds)
        self.dmap = {}
        self.kde = {}
        self.log_accept_frac = {}

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.loglr_region = float(loglr_region)

    def fit_kde(self, i):
        """ Fit a KDE proposal to tile i's pregenerated prior samples, and
        estimate what fraction of its (unbounded) probability mass falls
        within the prior's support.

        We draw from this KDE by rejection (see draw_samples_from_bin), so
        the density of an *accepted* draw is the KDE's density truncated
        to the prior's support and renormalized by that accepted fraction,
        not the KDE's raw density. Sampling code needs this fraction to
        weight tiles correctly relative to one another: tiles whose KDE
        support spills further outside the prior (more rejection) would
        otherwise be systematically over-weighted relative to tiles with
        less spillover.
        """
        names = self.dtype.names
        data = numpy.column_stack(
            [self.dmap[i][p].astype(float) for p in names]).T
        try:
            kde = gaussian_kde(data)
        except numpy.linalg.LinAlgError:
            # a handful of pregenerated points can give a singular
            # covariance; nudge them apart and retry
            jitter = data + numpy.random.normal(scale=1e-8, size=data.shape)
            kde = gaussian_kde(jitter)

        nprobe = 5000
        probe = kde.resample(nprobe).T
        params = {p: probe[:, j] for j, p in enumerate(names)}
        accept_frac = numpy.isfinite(self._logprior(params)).mean()
        accept_frac = max(accept_frac, 1.0 / nprobe)
        return kde, numpy.log(accept_frac)

    def _logprior(self, params):
        """ Evaluate the model's prior log-density for a batch of points.

        Some prior distributions (e.g. uniform's boundary check) aren't
        vectorized, so points are evaluated one at a time. This is cheap
        since prior evaluation involves no waveform generation, unlike
        the likelihood.
        """
        names = list(params.keys())
        n = len(params[names[0]])
        out = numpy.zeros(n)
        for k in range(n):
            pset = {p: params[p][k] for p in names}
            out[k] = self.model.prior_distribution(**pset)
        return out

    def draw_samples_from_bin(self, i, size):
        """ Draw fresh proposal samples for tile i from its KDE.

        KDE draws that fall outside the prior's support are rejected and
        redrawn, since the KDE's support can extend slightly past the
        prior's boundaries.
        """
        if size == 0:
            return FieldArray(0, dtype=self.dtype)

        t0 = time.time()
        names = self.dtype.names
        accepted = numpy.zeros((0, len(names)))
        tries = 0
        while len(accepted) < size and tries < 20:
            ndraw = max(size - len(accepted), 1)
            draws = self.kde[i].resample(ndraw).T
            params = {p: draws[:, j] for j, p in enumerate(names)}
            logp = self._logprior(params)
            good = numpy.isfinite(logp)
            accepted = numpy.vstack([accepted, draws[good]])
            tries += 1

        dt = time.time() - t0
        if tries > 1 or dt > 0.5:
            logging.info('bin %i: drew %i in %i tries, %.2f s '
                         '(kde npoints=%i)', i, size, tries, dt,
                         self.kde[i].n)

        out = FieldArray(size, dtype=self.dtype)
        for j, p in enumerate(names):
            out[p] = accepted[:size, j]
        return out

    def sample_round(self, bin_weight, node_idx):
        """ Sample from the posterior using KDE proposals fit to each
        tile's pregenerated prior samples, and the weighting factor of
        each tile of the prior space.

        bin_weight: Array
            The weighting importance factor of each bin of the prior space
        node_idx: Array
            The set of ids into the prebinned prior volume to use. This should
            map to the given weights.
        """
        logging.info("...draw template bins")
        drawcount = (bin_weight * self.target_likelihood_calls).astype(int)
        total_draw = drawcount.sum()

        logging.info('...drawn random points within bins')
        t_draw0 = time.time()
        psamp = FieldArray(total_draw, dtype=self.dtype)
        logkde = numpy.zeros(total_draw, dtype=float)
        bin_id = numpy.zeros(total_draw, dtype=int)
        j = 0
        t_logpdf = 0.0
        for i, c in enumerate(drawcount):
            node = node_idx[i]
            bdraw = self.draw_samples_from_bin(node, c)
            psamp[j:j+len(bdraw)] = bdraw
            if len(bdraw) > 0:
                t0 = time.time()
                pts = numpy.column_stack(
                    [bdraw[p].astype(float) for p in self.dtype.names]).T
                # effective density of an accepted (rejection-sampled)
                # draw = raw kde density / accept_frac
                logkde[j:j+len(bdraw)] = (self.kde[node].logpdf(pts)
                                          - self.log_accept_frac[node])
                t_logpdf += time.time() - t0
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        logging.info("Templates drawn from %s", len(node_idx))
        logging.info("Values drawn %s", len(psamp))
        logging.info('...draw+logpdf stage took %.2f s (logpdf: %.2f s)',
                     time.time() - t_draw0, t_logpdf)

        # Calculate the likelihood values for the drawn parameter space
        # points
        t0 = time.time()
        args = []
        for i in range(len(psamp)):
            pset = {p: psamp[p][i] for p in self.model.variable_params}
            args.append(pset)

        loglr_samp = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                    total=len(args)))
        loglr_samp = numpy.array(loglr_samp)
        logging.info('...likelihood stage took %.2f s', time.time() - t0)

        # Standard importance-sampling correction for the KDE proposal:
        # weight = likelihood * true_prior / kde_proposal_density. This is
        # a strict generalization of the finite-pool weighting used by
        # the 'games' sampler: in the limit of a perfect-fit KDE it
        # reduces to the same result, since a KDE fit to genuine prior
        # draws integrates to (approximately) the tile's true conditional
        # prior density.
        logprior = self._logprior(
            {p: psamp[p] for p in self.model.variable_params})
        logw3 = loglr_samp + logprior - logkde - numpy.log(drawcount[bin_id])

        import os
        if os.environ.get('GAMES2_DUMP'):
            psamp_arrs = {f'psamp_{p}': numpy.asarray(psamp[p])
                         for p in self.dtype.names}
            numpy.savez(os.environ['GAMES2_DUMP'],
                       loglr_samp=loglr_samp, logprior=logprior,
                       logkde=logkde, bin_id=bin_id,
                       drawcount=drawcount, node_idx=node_idx,
                       **psamp_arrs)
            logging.info('dumped round diagnostics to %s',
                         os.environ['GAMES2_DUMP'])

        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)
        return psamp, loglr_samp, weight2, bin_id

    def run(self):
        """ Produce posterior samples """
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as mapfile:
            bparams = {p: mapfile['bank'][p][:] for p in self.variable_params}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            lengths = numpy.array([len(mapfile['map'][str(x)])
                                   for x in range(num_nodes)])
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

        logging.info('Drawing proposal samples from node regions')
        passed = (node_loglrs > loglr_bound) & ~numpy.isnan(node_loglrs)
        passed = numpy.where(passed)[0]

        logging.info("...reading template bins")
        t0 = time.time()
        with h5py.File(self.mapfile, 'r') as mapfile:
            for i in passed:
                self.dmap[i] = mapfile['map'][str(i)][:]
        logging.info('...read %i tiles in %.2f s', len(passed),
                     time.time() - t0)

        t0 = time.time()
        for i in passed:
            t1 = time.time()
            self.kde[i], self.log_accept_frac[i] = self.fit_kde(i)
            dt = time.time() - t1
            if dt > 0.5:
                logging.info('fit_kde bin %i: %i points, %.2f s',
                             i, len(self.dmap[i]), dt)
        logging.info('...fit %i KDEs in %.2f s', len(passed),
                     time.time() - t0)

        # Sample from posterior
        psamp = None
        loglr_samp = None
        weight2 = None
        bin_ids = None

        weight = lengths[passed] / lengths[passed].sum()

        for i in range(self.rounds):
            psamp_v, loglr_samp_v, weight2_v, bin_id = \
                self.sample_round(weight / weight.sum(), passed)

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

        # Prepare the equally weighted output samples
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
