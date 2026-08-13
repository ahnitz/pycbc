""" IO for the games importance samplers, which produce WEIGHTED draws.

The pycbc convention for a weighted sampler -- set by ``dynesty`` -- is that
the file stores the raw draws together with a ``logwt`` dataset and a
``log_evidence`` attribute, and the RESAMPLING TO EQUAL WEIGHT HAPPENS ON
READ. Every downstream tool goes through ``read_samples``, so plotting and
posterior extraction then work without any of them knowing about weights.

Storing the raw draws is what makes the diagnostics that matter available
after the fact: the effective sample size, the evidence and its error, and
the Pareto tail index of the weights all need the individual weights, and
none of them can be recovered from samples that have already been resampled.
"""
import numpy

from .base_hdf import BaseInferenceFile
from .posterior import read_raw_samples_from_file, write_samples_to_file


class GamesFile(BaseInferenceFile):
    """ Weighted draws from the games samplers. """

    name = 'games_file'

    def read_raw_samples(self, fields, raw_samples=False, seed=0):
        """ Draws resampled to equal weight, or the raw weighted draws.

        Parameters
        ----------
        fields : list of str
            Names of the datasets in the file's ``samples`` group to read.
        raw_samples : bool, optional
            Return the weighted draws as stored, without resampling.
        seed : int, optional
            Seed for the resampling, so a plot is reproducible.

        Returns
        -------
        dict : parameter name -> samples.
        """
        samples = read_raw_samples_from_file(self, fields)
        if raw_samples:
            return samples
        logwt = read_raw_samples_from_file(self, ['logwt'])['logwt']
        # Systematic resampling, as DynestyFile does: one stratified pass
        # rather than independent draws, which removes the extra multinomial
        # variance that plain numpy.random.choice would add on top of the
        # weights' own.
        numpy.random.seed(seed)
        w = numpy.exp(logwt - logwt.max())
        n = len(w)
        pos = (numpy.random.random() + numpy.arange(n)) / n
        cdf = numpy.cumsum(w)
        cdf /= cdf[-1]
        idx = numpy.searchsorted(cdf, pos)
        idx = numpy.clip(idx, 0, n - 1)
        return {p: samples[p][idx] for p in samples}

    def write_samples(self, samples, parameters=None):
        return write_samples_to_file(self, samples, parameters=parameters)

    def write_sampler_metadata(self, sampler):
        sampler.model.write_metadata(self)

    def write_resume_point(self):
        pass

    write_run_start_time = write_run_end_time = write_resume_point
