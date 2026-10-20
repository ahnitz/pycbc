"""
Chisq based on sine-gaussian tiles.
See https://arxiv.org/abs/1709.08974 for a discussion.
"""

from collections import OrderedDict

import numpy

from pycbc.events import ranking
from pycbc.filter import sigma
from pycbc.types import FrequencySeries
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.waveform import sinegauss
from pycbc.waveform.utils import apply_fseries_time_shift


class SingleDetSGChisq(SingleDetPowerChisq):
    """
    Class that handles precomputation and memory management for efficiently
    running the sine-Gaussian chisq
    """

    name = "sg"
    columns = {"sg_chisq": numpy.float32}
    default_values = {"sg_chisq": 1}
    returns = {"sg_chisq": numpy.float32}

    def __init__(self, bank, num_bins=0, snr_threshold=None, chisq_locations=None):
        """
        Create sine-Gaussian Chisq Calculator

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
        self.snr_threshold = snr_threshold
        if chisq_locations is not None:
            self.do = True
            self.num_bins = num_bins
            self.params = {}
            for descr in chisq_locations:
                region, values = descr.split(":")
                mask = bank.table.parse_boolargs([(1, region), (0, "else")])[0]
                hashes = bank.table["template_hash"][mask.astype(bool)]
                for h in hashes:
                    self.params[h] = values
        else:
            self.do = False

    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Sine-Gaussian Chisq")
        group.add_argument(
            "--sgchisq-snr-threshold",
            type=float,
            help="Minimum SNR threshold to use SG chisq",
        )
        group.add_argument(
            "--sgchisq-locations",
            type=str,
            nargs="+",
            help="Frequency offsets and quality factors of the sine-Gaussians"
            " to use, format 'region-boolean:q1-offset1,q2-offset2'. "
            "Offset is relative to the end frequency of the approximant."
            " Region is a boolean expression selecting templates to "
            "apply the sine-Gaussians to, ex. 'mtotal>40'",
        )

    @classmethod
    def from_cli(cls, args, bank=None, chisq_bins=None):
        if chisq_bins is None:
            chisq_bins = str(args.chisq_bins)
        return cls(bank, chisq_bins, args.sgchisq_snr_threshold, args.sgchisq_locations)

    def default_steps(self):
        steps = []
        if self.snr_threshold is not None:
            steps.append(f"newsnr:{float(self.snr_threshold)!r}:lower_inc:skip")
        steps.append(self.name)
        return steps

    def evaluate(self, batch, j):
        stilde = batch.stilde
        template = batch.template(j)
        val = self.values(
            stilde,
            template,
            stilde.psd,
            numpy.array([batch.snr[j]]),
            1,
            indices=numpy.array([batch.time_index[j]]),
            template_whitened=True,
            whitened=getattr(stilde, "_is_whitened", False),
        )
        return {"sg_chisq": val[0]}

    def values(
        self,
        stilde,
        template,
        psd,
        snrv,
        snr_norm,
        bchisq=None,
        bchisq_dof=None,
        indices=None,
        template_whitened=False,
        whitened=True,
    ):
        """
        Calculate sine-Gaussian chisq

        Parameters
        ----------
        stilde: pycbc.types.FrequencySeries
            The data, divided by its PSD (overwhitened) unless whitened is
            False
        template: pycbc.types.FrequencySeries
            The waveform template being analyzed
        psd: pycbc.types.FrequencySeries
            The power spectral density of the data
        snrv: numpy.ndarray
            The peak unnormalized complex SNR values
        snr_norm: float
            The normalization factor for the snr
        bchisq: numpy.ndarray, optional
            The Bruce Allen power chisq values for these triggers
        bchisq_dof: numpy.ndarray, optional
            The degrees of freedom of the Bruce chisq
        indices: numpy.ndarray, optional
            The indices of the snr peaks.
        template_whitened: bool, optional
            Whether template is already divided by sqrt(psd).
        whitened: bool, optional
            Whether stilde is already divided by the PSD. If not, only the
            frequencies the tiles read are divided by it.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index

        """
        if not self.do:
            return None

        chisq = numpy.ones(len(snrv))
        template_hash = template.params.template_hash
        if template_hash not in self.params:
            return chisq

        if bchisq is not None and self.snr_threshold is not None:
            # Skip if newsnr too low
            snr = abs(numpy.array(snrv) * snr_norm)
            rchisq = numpy.array(bchisq) / numpy.array(bchisq_dof)
            nsnr = ranking.newsnr(snr, rchisq)
            above = numpy.flatnonzero(nsnr >= self.snr_threshold)
            if not len(above):
                return chisq
        else:
            above = range(len(snrv))

        tiles = self._cached_tiles(
            stilde,
            psd,
            template_hash,
            self.peak_frequency(template, None if template_whitened else psd),
            template.f_lower,
        )
        if tiles is None:
            return chisq

        # Only the frequencies the tiles cover are read. They are time
        # shifted in double precision: the phase f * t of the shift reaches
        # ~1e5 cycles, which single precision cannot hold accurately
        klo = min(tile[1] for tile in tiles)
        khi = max(tile[2] for tile in tiles)
        band = stilde.numpy()[klo:khi].astype(numpy.complex128)
        if not whitened:
            band /= psd.numpy()[klo:khi]
        band = FrequencySeries(band, delta_f=stilde.delta_f)

        num_samples = (len(stilde) - 1) * 2
        dt = 1.0 / (num_samples * stilde.delta_f)
        for i in above:
            time = float(template.epoch) + dt * indices[i]
            # Shift the data so the trigger time is at 0. The band starts at
            # frequency klo * delta_f, so the shift of the band's first
            # sample is applied separately
            shifted = apply_fseries_time_shift(band, -time).numpy() * numpy.exp(
                2j * numpy.pi * klo * stilde.delta_f * time
            )
            chisq[i] = 0
            for gtem, kmin, kmax, gsigma in tiles:
                gsnr = (
                    gtem.numpy()[kmin:kmax] * shifted[kmin - klo : kmax - klo]
                ).sum()
                gsnr *= 4.0 * gtem.delta_f / gsigma
                chisq[i] += abs(gsnr) ** 2.0
            chisq[i] /= 2 * len(tiles)
        return chisq

    def peak_frequency(self, template, psd=None):
        """
        Estimate the frequency up to which the template has power, from
        which the sine-Gaussian tiles are offset.

        The power per frequency is assumed constant over the last two chisq
        bins, so the estimate is the second-to-last bin edge plus the width
        of the second-to-last bin. The last bin edge is not used, as it need
        not be where the waveform ends.
        """
        bins = self.cached_chisq_bins(template, psd)
        fstep = bins[-2] - bins[-3]
        return float((bins[-2] + fstep) * template.delta_f)

    def _cached_tiles(self, stilde, psd, template_hash, fpeak, f_lower):
        """
        Return the result of _tiles, cached for the 16 templates used
        most recently. Triggers are often evaluated one at a time, in an
        order that moves between templates.
        """
        cache = self.__dict__.setdefault("_tile_cache", OrderedDict())
        key = (template_hash, fpeak, f_lower, len(stilde), stilde.delta_f)
        entry = cache.get(key)
        # The psd is held by the entry, so it cannot be another one reusing
        # its id
        if entry is not None and entry[0] is psd:
            cache.move_to_end(key)
            return entry[1]
        tiles = self._tiles(stilde, psd, self.params[template_hash], fpeak, f_lower)
        cache[key] = (psd, tiles)
        if len(cache) > 16:
            cache.popitem(last=False)
        return tiles

    @staticmethod
    def _tiles(stilde, psd, values, fpeak, f_lower):
        """
        Return, for each sine-Gaussian tile, the tile, the first and last
        frequency index it covers and its norm; or None if any tile would
        reach within 10% of the data's Nyquist frequency.
        """
        # Only apply the sine-Gaussian in a +-50 Hz range around the
        # central frequency
        qwindow = 50
        # 90% of the Nyquist frequency of the data, which avoids issues
        # near Nyquist due to resample filtering
        fstop = len(stilde) * stilde.delta_f * 0.9
        delta_f = stilde.delta_f
        kmin = int(f_lower / psd.delta_f)
        tiles = []
        for descr in values.split(","):
            q, offset = descr.split("-")
            q, offset = float(q), float(offset)
            fcen = fpeak + offset
            flow = max(kmin * delta_f, fcen - qwindow)
            fhigh = fcen + qwindow
            if fhigh > fstop:
                return None
            kmin = int(flow / delta_f)
            kmax = int(fhigh / delta_f)
            gtem = sinegauss.fd_sine_gaussian(
                1.0, q, fcen, flow, len(stilde) * delta_f, delta_f
            ).astype(numpy.complex64)
            gsigma = sigma(
                gtem, psd=psd, low_frequency_cutoff=flow, high_frequency_cutoff=fhigh
            )
            tiles.append((gtem, kmin, kmax, gsigma))
        return tiles
