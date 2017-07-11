""" Chisq based on sine-gaussian tiles """

import numpy
import logging

from pycbc.waveform.utils import apply_fseries_time_shift
from pycbc.filter import sigma
from pycbc.waveform import sinegauss
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.events import newsnr

class SingleDetSGChisq(SingleDetPowerChisq):
    """Class that handles precomputation and memory management for efficiently
    running the sine-Gaussian chisq
    """
    returns = {'sg_chisq': numpy.float32}
    
    def __init__(self, bank, num_bins=0,
                       snr_threshold=None,
                       chisq_locations=None):
        """ Create sine-Gaussian Chisq Calculator

        Parameters
        ----------
        bank: pycbc.waveform.TemplateBank
            The template bank that will be processed.
        num_bins: str
            The string determining the number of power chisq bins
        snr_threshold: float
            The threshold to calculate the sine-Gaussian chisq
        chisq_locations: list of strs
            List of strings which detail where to place a sine-Gaussian.
            The format is 'region-boolean:q1-offset1,q2-offset2'.
            The offset is relative to the end frequency of the approximant.
            The region is a boolean expresion such as 'mtotal>40' and indicates
            where to apply this set of sine-Gaussians.
        """
        if snr_threshold is not None:
            self.do = True
            self.num_bins = num_bins
            self.snr_threshold = snr_threshold
            self.params = {}
            for descr in chisq_locations:
                region, values = descr.split(":")
                mask = bank.table.parse_boolargs([(1, region), (0, 'else')])[0]
                hashes = bank.table['template_hash'][mask.astype(bool)]
                for h in hashes:
                    self.params[h] = values
        else:
            self.do = False
            
    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Sine-Gaussian Chisq")
        group.add_argument("--sgchisq-snr-threshold", type=float,
            help="Minimum SNR threshold to use SG chisq")
        group.add_argument("--sgchisq-locations", type=str, nargs='+',
            help="The frequencies and quality factors of the sine-gaussians"
                 " to use. The format is 'region-boolean:q1-offset1,q2-offset2'."
                 "The offset is relative to the end frequency of the approximant."
                 "The region is a boolean expresion such as 'mtotal>40' and indicates "
                 "where to apply this set of sine-Gaussians.")

    @classmethod
    def from_cli(cls, args, bank, chisq_bins):
        return cls(bank, chisq_bins,
                   args.sgchisq_snr_threshold,
                   args.sgchisq_locations)

    def values(self, stilde, template, psd, snrv, snr_norm,
                     bchisq, bchisq_dof, indices):
        """ Calculate sine-Gaussian chisq

        Parameters
        ----------
        stilde: pycbc.types.Frequencyseries
            The overwhitened strain
        template: pycbc.types.Frequencyseries
            The waveform template being analyzed
        psd: pycbc.types.Frequencyseries
            The power spectral density of the data
        snrv: numpy.ndarray
            The peak unnormalized complex SNR values
        snr_norm: float
            The normalization factor for the snr
        bchisq: numpy.ndarray
            The Bruce Allen power chisq values for these triggers
        bchisq_dof: numpy.ndarray
            The degrees of freedom of the Bruce chisq
        indics: numpy.ndarray
            The indices of the snr peaks.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index
        """
        if not self.do:
            return None

        if template.params.template_hash not in self.params:
            return numpy.ones(len(snrv))

        # This is implemented slowly, so let's not call it often, OK?
        w = 0.25
        stilde = stilde * psd
        chisq = numpy.ones(len(snrv))
        for i, snrvi in enumerate(snrv):
            #Skip if newsnr too low
            snr = abs(snrvi * snr_norm)
            nsnr = newsnr(snr, bchisq[i] / bchisq_dof[i])
            if nsnr < self.snr_threshold:
                continue

            # Subtract template from data        
            time = float(template.epoch) + dt * indices[i]
            scale = snrv[i] * snr_norm / (template.sigmasq(psd)) ** 0.5
            stilde -= scale * apply_fseries_time_shift(template, time)
            ts = stilde.to_timeseries()

            psd = ts.time_slice(time - w * 16, time + w * 16).psd(w)
            p = ts.time_slice(et - w / 2, et + w / 2).psd(w)
            kmin = int(20 / p.delta_f)
            kmax = int(400 / p.delta_f)
            pwr = (p / psd)[kmin:kmax].sum() * p.delta_f / 380

            chisq[i] = pwr
            logging.info('Found Residucal chisq %s', chisq[i])
        return chisq

