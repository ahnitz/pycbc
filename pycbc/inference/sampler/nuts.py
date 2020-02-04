from .base import BaseSampler

class NutsSampler(object):
    name = 'nuts'

    def __init__(self, model):
        from nuts import nuts6
        self.model = model
        self.checkpoint_file = None
        self.backup_file = None
        self.checkpoint_valid = None
        self.new_checkpoint = None
        self.model_stat = None
        self.io = None
        self.samples = None
        
        self.sampler = nuts6
        self.Madapt = 
        self.steps = 5000
        self.step_size = 0.2

    def from_config(cls, cp, model, nprocesses=1, use_mpi=False):
        """This should initialize the sampler given a config file.
        """
        pass

    def run(self):
        """This function should run the sampler.

        Any checkpointing should be done internally in this function.
        """
        self.samples, lnprob, epsilon = self.smampler(self.logp, 

    def set_initial_conditions(self, initial_distribution=None,
                               samples_file=None):
        """Sets up the starting point for the sampler.

        Should also set the sampler's random state.
        """
        pass

    def checkpoint(self):
        """The sampler must have a checkpoint method for dumping raw samples
        and stats to the file type defined by ``io``.
        """
        pass

    def finalize(self):
        """Do any finalization to the samples file before exiting."""
        pass

    def setup_output(self, output_file, force=False):
        """Sets up the sampler's checkpoint and output files.
        """
        pass

