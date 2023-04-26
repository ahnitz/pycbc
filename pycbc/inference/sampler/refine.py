""" Sampler that uses kde refinement of an existing posterior estimate.
"""

import logging
import numpy
import numpy.random

from scipy.special import logsumexp
from scipy.stats import gaussian_kde
from scipy.stats import entropy as sentropy

from pycbc.inference.io import PosteriorFile
from pycbc.inference import models
from pycbc.pool import choose_pool
from pycbc.inference.io import loadfile
from pycbc.io.record import FieldArray

from .base import setup_output, initial_dist_from_config
from .dummy import DummySampler

class RefineSampler(DummySampler):
    """Sampler for kde drawn refinement of existing posterior estimate

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    num_samples: int
        The number of samples to draw from the kde at the conclusion
    iterative_kde_samples: int
        The number of samples to add to the kde during each iterations
    min_refinement_steps: int
        The minimum number of iterations to take.
    max_refinement_steps: The maximum number of refinment steps to take.
    entropy: float
        The target entropy between iterative kdes
    dlogz: float
        The target evidence difference between iterative kde updates
    kde: kde
        The inital kde to use.
    """
    name = 'refine'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 num_samples=int(1e5),
                 iterative_kde_samples = int(1e3),
                 min_refinement_steps = 5,
                 max_refinement_steps = 20,
                 entropy = 0.001,
                 dlogz = 0.01,
                 kde = None,
                 **kwargs):
        super().__init__(model, *args)

        self.model = model
        self.kde = kde
        self.vparam = model.variable_params
        models._global_instance = model
        self.num_samples = int(num_samples)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.num_samples = int(num_samples)
        self.iterative_kde_samples = int(iterative_kde_samples)
        self.min_refinement_steps = int(min_refinement_steps)
        self.max_refinement_steps = int(max_refinement_steps)
        self.entropy = float(entropy)
        self.dlogz_target = float(dlogz)

    def draw_samples(self, size):
        """Draw new samples within the model priors"""
        ksamples = self.kde.resample(size=size)
        params = {k:ksamples[i,:] for i, k in enumerate(self.vparam)}
        keep = self.model.prior_distribution.contains(params)
        return ksamples[:, keep]
        
    @staticmethod
    def compare_kde(kde1, kde2, size=int(1e4)):
        s = kde1.resample(size=size)
        return sentropy(kde1.pdf(s), kde2.pdf(s))

    def converged(self, step, kde_new, factor):   
        """ Check that kde is converged by comparing to previous iteration
        """  
        if not hasattr(self, 'old_logz'):
            self.old_logz =  numpy.inf
        
        entropy_diff = self.compare_kde(self.kde, kde_new)
        
        # Compare how the logz changes when adding new samples
        # this is guaranteed to decrease as old samples included
        logz = logsumexp(factor) - numpy.log(len(factor))
        dlogz = logz - self.old_logz
        self.old_logz = logz
        
        # compare evidence subsets agree
        choice2 = numpy.random.choice(factor, len(factor) // 2)
        logz2 = logsumexp(choice2) - numpy.log(len(choice2))
        choice3 = numpy.random.choice(factor, len(factor) // 2)
        logz3 = logsumexp(choice3) - numpy.log(len(choice3))
        dlogz2 = logz3 - logz2
        
        logging.info('%s: Checking convergence: dlogz_iter=%.4f, dlogz_half=%.4f, entropy=%.4f',
                     step,  dlogz, dlogz2, entropy_diff)
        if (entropy_diff < self.entropy and step >= self.min_refinement_steps
            and max(abs(dlogz), abs(dlogz2)) < self.dlogz_target):
            return True
        else:
            return False

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """This should initialize the sampler given a config file.
        """
        kwargs = {k: cp.get('sampler', k) for k in cp.options('sampler')}
        obj = cls(model, nprocesses=nprocesses, use_mpi=use_mpi, **kwargs)
        obj.set_start_from_config(cp)
        setup_output(obj, output_file, check_nsamples=False, validate=False)
        return obj

    def set_start_from_config(self, cp):
        """Sets the initial state of the sampler from config file
        """
        if cp.has_option('sampler', 'start-file'):
            start_file = cp.get('sampler', 'start-file')
            logging.info("Using file %s for initial positions", start_file)
            samples = loadfile(start_file, 'r').read_samples(self.vparam)
        else:
            init_prior = initial_dist_from_config(
                cp, self.model.variable_params, self.model.static_params)
            ksamples = init_prior.rvs(size=self.iterative_kde_samples)

        ksamples = numpy.array([samples[v] for v in self.vparam])
        self.kde = gaussian_kde(ksamples)                

    def run(self):
        """ Iterative sample from kde and update based on likelihood values
        """
        total_samples = None
        total_logp = None
        total_logw = None
        total_logl = None

        for r in range(self.max_refinement_steps):
            ksamples = self.draw_samples(self.iterative_kde_samples)

            # Calculate likelihood for each samples
            logp = []
            logl = []
            for i in range(len(ksamples[0])):
                param = {k: ksamples[j][i] for j, k in enumerate(self.vparam)}
                self.model.update(**param)
                logp.append(self.model.logposterior)
                logl.append(self.model.loglikelihood)

            logp = numpy.array(logp)
            logl = numpy.array(logl)
            logw = logp - self.kde.logpdf(ksamples)    
            
            if total_samples is not None:
                total_samples = numpy.concatenate([total_samples, ksamples], axis=1)
                total_logp = numpy.concatenate([total_logp, logp])
                total_logw = numpy.concatenate([total_logw, logw])
                total_logl = numpy.concatenate([total_logl, logl])
            else:
                total_samples = ksamples
                total_logp = logp
                total_logw = logw
                total_logl = logl

            ntotal_logw = total_logw - logsumexp(total_logw)
            kde_new = gaussian_kde(total_samples, weights=numpy.exp(ntotal_logw))
            
            if self.converged(r, kde_new, total_logl + total_logw):
                break
            
            self.kde = kde_new
            
        ksamples = self.draw_samples(self.num_samples)
        self._samples = {k: ksamples[k] for k in self.vparam}

