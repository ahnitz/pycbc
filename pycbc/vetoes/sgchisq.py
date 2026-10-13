"""Chisq based on sine-gaussian tiles.
See https://arxiv.org/abs/1709.08974 for a discussion.
"""

import numpy

from pycbc.waveform.utils import apply_fseries_time_shift
from pycbc.filter import sigma
from pycbc.waveform import sinegauss
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.events import ranking

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
            The region is a boolean expression such as 'mtotal>40' indicating
            which templates to apply this set of sine-Gaussians to.
        """
        if snr_threshold is not None and chisq_locations is not None:
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
            help="Frequency offsets and quality factors of the sine-Gaussians"
                 " to use, format 'region-boolean:q1-offset1,q2-offset2'. "
                 "Offset is relative to the end frequency of the approximant."
                 " Region is a boolean expression selecting templates to "
                 "apply the sine-Gaussians to, ex. 'mtotal>40'")

    @classmethod
    def from_cli(cls, args, bank, chisq_bins):
        return cls(bank, chisq_bins,
                   args.sgchisq_snr_threshold,
                   args.sgchisq_locations)

    def values_from_fpeak(self, stilde, thash, fpeak, psd, snrv, snr_norm,
                          bchisq, bchisq_dof, indices, flow=20.0, epoch=0.0,
                          template=None):
        """ Calculate sine-Gaussian chisq directly from template peak frequency
        and template hash, with optional template projection subtraction.

        Parameters
        ----------
        stilde: pycbc.types.FrequencySeries
            The overwhitened strain
        thash: int
            The template hash for looking up sine-Gaussian locations
        fpeak: float
            The estimated peak/termination frequency of the template in Hz
        psd: pycbc.types.FrequencySeries
            The power spectral density of the data
        snrv: numpy.ndarray
            The peak unnormalized complex SNR values
        snr_norm: float
            The normalization factor for the snr
        bchisq: numpy.ndarray
            The Bruce Allen power chisq values for these triggers
        bchisq_dof: numpy.ndarray
            The degrees of freedom of the Bruce chisq
        indices: numpy.ndarray
            The indices of the snr peaks.
        flow: float, optional
            Low frequency cutoff (default 20.0 Hz)
        epoch: float, optional
            Epoch offset for time shifting (default 0.0)
        template: pycbc.types.FrequencySeries, optional
            The waveform template being analyzed, used to subtract expected
            signal power projection from sine-Gaussian tiles.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index
        """
        if not self.do:
            return None

        if thash not in self.params:
            return numpy.ones(len(snrv))
        values = self.params[thash].split(',')

        fstop = len(stilde) * stilde.delta_f * 0.9
        delta_f = stilde.delta_f
        N = (len(stilde) - 1) * 2
        dt = 1.0 / (N * delta_f)
        kmin_base = int(flow / psd.delta_f)

        # If any sine-gaussian tile has an upper frequency near nyquist return 1
        for descr in values:
            q, offset = descr.split('-')
            if (fpeak + float(offset) + 50.0) > fstop:
                return numpy.ones(len(snrv))

        # Template norm sigma_h for proper unit-normalized projection
        if template is not None:
            if hasattr(template, 'sigmasq'):
                hsigma = float(numpy.sqrt(template.sigmasq(psd)))
            else:
                flow_tmpl = getattr(template, 'f_lower', float(template.delta_f))
                hsigma = float(sigma(template, psd=psd, low_frequency_cutoff=flow_tmpl))
        else:
            hsigma = 0.0

        chisq = numpy.ones(len(snrv))
        gtem = [None for _ in values]
        gsigma_list = [None for _ in values]
        gh_proj = [None for _ in values]
        for i, snrvi in enumerate(snrv):
            # Skip if newsnr too low
            snr = abs(snrvi * snr_norm)
            nsnr = ranking.newsnr(snr, bchisq[i] / bchisq_dof[i])
            if nsnr < self.snr_threshold:
                continue

            time = float(epoch) + dt * indices[i]
            # Shift the time of interest to be centered on 0
            stilde_shift = apply_fseries_time_shift(stilde, -time)

            # Only apply the sine-Gaussian in a +-50 Hz range around the
            # central frequency
            qwindow = 50
            chisq[i] = 0
            dof = 0

            # Calculate the sum of SNR^2 for the sine-Gaussians specified
            for idxx, descr in enumerate(values):
                # Get the q and frequency offset from the descriptor
                q, offset = descr.split('-')
                q, offset = float(q), float(offset)
                fcen = fpeak + offset
                t_flow = max(kmin_base * delta_f, fcen - qwindow)
                fhigh = fcen + qwindow

                kmin = int(t_flow / delta_f)
                kmax = int(fhigh / delta_f)

                # Calculate sine-gaussian tile and template projection
                if gtem[idxx] is None:
                    gtem[idxx] = sinegauss.fd_sine_gaussian(
                        1.0, q, fcen, t_flow,
                        len(stilde) * delta_f,
                        delta_f
                    ).astype(numpy.complex64)
                    gsigma_list[idxx] = sigma(gtem[idxx], psd=psd,
                                              low_frequency_cutoff=t_flow,
                                              high_frequency_cutoff=fhigh)
                    if template is not None and hsigma > 0:
                        g_tile = gtem[idxx][kmin:kmax]
                        tmpl_overwhite = template[kmin:kmax] / psd[kmin:kmax]
                        gh_proj[idxx] = (g_tile * tmpl_overwhite).sum() * (
                            4.0 * delta_f / (gsigma_list[idxx] * hsigma)
                        )
                    else:
                        gh_proj[idxx] = 0.0

                gsigma = gsigma_list[idxx]
                # Calculate the SNR of the tile
                gsnr = (gtem[idxx][kmin:kmax] * stilde_shift[kmin:kmax]).sum()
                gsnr *= 4.0 * delta_f / gsigma

                # Subtract expected template projection to prevent loud genuine
                # signals from self-vetoing due to high-frequency leakage
                if gh_proj[idxx] is not None and gh_proj[idxx] != 0.0:
                    z = snrvi * snr_norm
                    gsnr -= z * gh_proj[idxx]

                chisq[i] += abs(gsnr)**2.0
                dof += 2
            if dof == 0:
                chisq[i] = 1
            else:
                chisq[i] /= dof
        return chisq

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
        indices: numpy.ndarray
            The indices of the snr peaks.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index
        """
        if not self.do:
            return None

        thash = getattr(getattr(template, 'params', None), 'template_hash', None)
        if thash is None or thash not in self.params:
            return numpy.ones(len(snrv))

        # Get the chisq bins to use as the frequency reference point
        bins = self.cached_chisq_bins(template, psd)

        # Estimate the maximum frequency up to which the waveform has
        # power by approximating power per frequency
        # as constant over the last 2 chisq bins. We cannot use the final
        # chisq bin edge as it does not have to be where the waveform
        # terminates.
        fstep = (bins[-2] - bins[-3])
        fpeak = (bins[-2] + fstep) * template.delta_f

        return self.values_from_fpeak(
            stilde, thash, fpeak, psd, snrv, snr_norm,
            bchisq, bchisq_dof, indices,
            flow=float(getattr(template, 'f_lower', 20.0)),
            epoch=float(getattr(template, 'epoch', 0.0)),
            template=template
        )

