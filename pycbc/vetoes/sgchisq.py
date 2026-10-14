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
        if snr_threshold is not None and chisq_locations is not None:
            self.do = True
            self.num_bins = num_bins
            self.snr_threshold = snr_threshold
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
    def from_cli(cls, args, bank, chisq_bins):
        return cls(bank, chisq_bins, args.sgchisq_snr_threshold, args.sgchisq_locations)

    def values(
        self, stilde, template, psd, snrv, snr_norm, bchisq, bchisq_dof, indices
    ):
        """
        Calculate sine-Gaussian chisq

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

        chisq = numpy.ones(len(snrv))
        if template.params.template_hash not in self.params:
            return chisq

        # Skip triggers whose newsnr is too low
        snr = abs(numpy.array(snrv) * snr_norm)
        rchisq = numpy.array(bchisq) / numpy.array(bchisq_dof)
        nsnr = ranking.newsnr(snr, rchisq)
        above = nsnr >= self.snr_threshold
        if not above.any():
            return chisq

        num_samples = (len(template) - 1) * 2
        dt = 1.0 / (num_samples * template.delta_f)
        times = float(template.epoch) + dt * numpy.array(indices)[above]
        snrs = numpy.array(snrv)[above] * snr_norm
        chisq[above] = self.values_at(stilde, psd, template, psd, times, snrs=snrs)
        return chisq

    def peak_frequency(self, template, psd):
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
        return (bins[-2] + fstep) * template.delta_f

    def values_at(
        self, stilde, psd, template, template_psd, times, whitened=True, snrs=None
    ):
        """
        Calculate the sine-Gaussian chisq of a template at the given
        times.

        Unlike values, the template need not be on the data's frequency
        grid: it is only used to place the tiles, from its peak frequency
        and low frequency cutoff.

        A signal that matches the template also contributes to any tile that
        overlaps the template, by an amount that grows with the square of
        its SNR, so a loud signal could fail the test. If the triggers' SNRs
        are given, that contribution, the trigger's SNR times the
        normalized overlap <tile|template>, is subtracted from
        each tile's SNR. Subtracting it also removes part of the noise, a
        fraction |<tile|template>|^2 of the tile's noise power, so each
        tile's squared SNR is divided by 1 - |<tile|template>|^2 to keep the
        chisq's mean at 1 in noise.

        Parameters
        ----------
        stilde: pycbc.types.FrequencySeries
            The data, divided by its PSD (overwhitened) unless whitened is
            False
        psd: pycbc.types.FrequencySeries
            The power spectral density of the data
        template: pycbc.types.FrequencySeries
            The template, on any frequency grid, with its f_lower and params
            attributes; params.template_hash selects the tiles to use (see
            chisq_locations).
        template_psd: pycbc.types.FrequencySeries
            The PSD on the template's frequency grid: all ones if the
            template is already whitened. Used to find the template's peak
            frequency (see peak_frequency) and its overlap with the tiles.
        times: numpy.ndarray
            The times of the triggers, in seconds from the start of stilde
        whitened: bool, optional
            Whether stilde is already divided by the PSD. If not, only the
            frequencies the tiles read are divided by it.
        snrs: numpy.ndarray, optional
            The normalized complex SNR of each trigger. If given, the
            template's response is removed from each tile.

        Returns
        -------
        chisq: numpy.ndarray
            The sine-Gaussian chisq at each time. It is 1 if the template has
            no tiles, or if a tile would reach within 10% of the Nyquist
            frequency.

        """
        chisq = numpy.ones(len(times))
        template_hash = template.params.template_hash
        if template_hash not in self.params:
            return chisq
        project = snrs is not None
        tiles = self._cached_tiles(
            stilde,
            psd,
            template_hash,
            self.peak_frequency(template, template_psd),
            template.f_lower,
            template if project else None,
            template_psd,
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
        for i, time in enumerate(times):
            # Shift the data so the trigger time is at 0. The band starts at
            # frequency klo * delta_f, so the shift of the band's first
            # sample is applied separately
            shifted = apply_fseries_time_shift(band, -time).numpy() * numpy.exp(
                2j * numpy.pi * klo * stilde.delta_f * time
            )
            chisq[i] = 0
            for gtem, kmin, kmax, gsigma, overlap in tiles:
                gsnr = (
                    gtem.numpy()[kmin:kmax] * shifted[kmin - klo : kmax - klo]
                ).sum()
                gsnr *= 4.0 * gtem.delta_f / gsigma
                if project:
                    gsnr -= snrs[i] * overlap
                    chisq[i] += abs(gsnr) ** 2 / (1 - abs(overlap) ** 2)
                else:
                    chisq[i] += abs(gsnr) ** 2
            chisq[i] /= 2 * len(tiles)
        return chisq

    @staticmethod
    def _overlaps(tiles, template, template_psd, psd):
        """
        Return the normalized overlap of each tile with the template,
        <tile|template> / (|tile| |template|).

        The overlap is computed on the template's frequency grid. The
        template is whitened by template_psd, and the tile by the data's
        PSD at the same frequencies.
        """
        df = template.delta_f
        # Only positive frequencies below the data's Nyquist frequency
        npos = min(len(template), int(len(psd) * psd.delta_f / df))
        tpsd = numpy.asarray(template_psd)[:npos]
        h = numpy.zeros(npos, dtype=complex)
        good = tpsd > 0
        h[good] = numpy.asarray(template)[:npos][good] / numpy.sqrt(tpsd[good])
        hnorm = numpy.sqrt(4 * df * numpy.sum(abs(h) ** 2))

        def on_template_grid(series, k):
            """
            Return series, given on the data's frequency grid, at the
            template's frequency samples k.
            """
            ratio = df / series.delta_f
            if abs(ratio - round(ratio)) < 1e-9:
                return numpy.asarray(series)[(k * round(ratio))]
            data_f = numpy.arange(len(series)) * series.delta_f
            return numpy.interp(k * df, data_f, numpy.asarray(series))

        out = []
        for gtem, kmin, kmax, _ in tiles:
            fmin, fmax = kmin * gtem.delta_f, kmax * gtem.delta_f
            k = numpy.arange(
                int(numpy.ceil(fmin / df)), min(int(numpy.ceil(fmax / df)), npos)
            )
            g = on_template_grid(gtem, k) / numpy.sqrt(on_template_grid(psd, k))
            gnorm = numpy.sqrt(4 * df * numpy.sum(abs(g) ** 2))
            if gnorm == 0 or hnorm == 0:
                out.append(0.0)
                continue
            out.append(4 * df * numpy.sum(g * h[k]) / (gnorm * hnorm))
        return out

    def _cached_tiles(
        self,
        stilde,
        psd,
        template_hash,
        fpeak,
        f_lower,
        template=None,
        template_psd=None,
    ):
        """
        Return the result of _tiles, with each tile's overlap with the
        template appended (0 if no template is given), cached for the 16
        templates used most recently. Triggers are often evaluated one at a
        time, in an order that moves between templates.
        """
        cache = self.__dict__.setdefault("_tile_cache", OrderedDict())
        key = (
            template_hash,
            fpeak,
            f_lower,
            len(stilde),
            stilde.delta_f,
            template is not None,
        )
        entry = cache.get(key)
        # The psd is held by the entry, so it cannot be another one reusing
        # its id
        if entry is not None and entry[0] is psd:
            cache.move_to_end(key)
            return entry[1]
        tiles = self._tiles(stilde, psd, self.params[template_hash], fpeak, f_lower)
        if tiles is not None:
            overlaps = (
                [0] * len(tiles)
                if template is None
                else self._overlaps(tiles, template, template_psd, psd)
            )
            tiles = [(*tile, o) for tile, o in zip(tiles, overlaps, strict=False)]
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
