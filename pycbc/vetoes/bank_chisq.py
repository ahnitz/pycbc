# Copyright (C) 2013  Ian Harry, Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
import logging
from math import sqrt

import numpy

from pycbc.filter import matched_filter_core, overlap_cplx, sigmasq
from pycbc.types import (
    Array,
    FrequencySeries,
    TimeSeries,
    real_same_precision_as,
    zeros,
)
from pycbc.waveform import FilterBank

logger = logging.getLogger(__name__)


def segment_snrs(filters, stilde, psd, low_frequency_cutoff, whitened=True):
    """
    This functions calculates the snr of each bank veto template against
    the segment

    Parameters
    ----------
    filters: list of FrequencySeries
        The list of bank veto templates filters.
    stilde: FrequencySeries
        The current segment of data.
    psd: FrequencySeries
    low_frequency_cutoff: float
    whitened: bool, optional
        Whether stilde is already divided by the PSD.

    Returns
    -------
    snr (list): List of snr time series.
    norm (list): List of normalizations factors for the snr time series.

    """
    snrs = []
    norms = []

    for bank_template in filters:
        # For every template compute the snr against the stilde segment
        h_norm = (
            bank_template.sigmasq(psd)
            if hasattr(bank_template, "sigmasq")
            else sigmasq(
                bank_template, psd=psd, low_frequency_cutoff=low_frequency_cutoff
            )
        )
        snr, _, norm = matched_filter_core(
            bank_template,
            stilde,
            h_norm=h_norm,
            psd=None if whitened else psd,
            low_frequency_cutoff=low_frequency_cutoff,
        )
        # SNR time series stored here
        snrs.append(snr)
        # Template normalization factor stored here
        norms.append(norm)

    return snrs, norms


def template_overlaps(
    bank_filters, template, psd, low_frequency_cutoff, template_whitened=False
):
    """
    This functions calculates the overlaps between the template and the
    bank veto templates.

    Parameters
    ----------
    bank_filters: List of FrequencySeries
    template: FrequencySeries
    psd: FrequencySeries
    low_frequency_cutoff: float
    template_whitened: bool, optional
        Whether template is already divided by sqrt(psd).

    Returns
    -------
    overlaps: List of complex overlap values.

    """
    overlaps = []
    if not template_whitened:
        t_ow = template / psd
        t_sigmasq = (
            template.sigmasq(psd)
            if hasattr(template, "sigmasq")
            else sigmasq(template, psd=psd, low_frequency_cutoff=low_frequency_cutoff)
        )
        dec = 1
        n_pos = len(template)
    else:
        f_nyq = (len(psd) - 1) * psd.delta_f
        n_pos = min(len(template), round(f_nyq / template.delta_f) + 1)
        dec = max(1, round(template.delta_f / psd.delta_f))
        t_pos = FrequencySeries(template[:n_pos], delta_f=template.delta_f)
        b_psd = FrequencySeries(psd[::dec][:n_pos], delta_f=template.delta_f)
        t_ow = t_pos / b_psd**0.5
        t_sigmasq = sigmasq(t_pos, low_frequency_cutoff=low_frequency_cutoff)

    for bank_template in bank_filters:
        if not template_whitened:
            b_pos = bank_template
            b_sigmasq = (
                bank_template.sigmasq(psd)
                if hasattr(bank_template, "sigmasq")
                else sigmasq(
                    bank_template,
                    psd=psd,
                    low_frequency_cutoff=low_frequency_cutoff,
                )
            )
        else:
            b_pos = FrequencySeries(
                bank_template[::dec][:n_pos], delta_f=template.delta_f
            )
            b_sigmasq = sigmasq(
                b_pos, psd=b_psd, low_frequency_cutoff=low_frequency_cutoff
            )
        overlap = overlap_cplx(
            t_ow,
            b_pos,
            low_frequency_cutoff=low_frequency_cutoff,
            normalized=False,
        )
        norm = sqrt(1 / t_sigmasq / b_sigmasq)
        overlaps.append(overlap * norm)
        if abs(overlaps[-1]) > 0.99:
            errMsg = "Overlap > 0.99 between bank template and filter. "
            errMsg += "This bank template will not be used to calculate "
            errMsg += "bank chisq for this filter template. The expected "
            errMsg += "value will be added to the chisq to account for "
            errMsg += "the removal of this template.\n"
            errMsg += "Masses of filter template: %e %e\n" % (
                getattr(template.params, "mass1", float("nan")),
                getattr(template.params, "mass2", float("nan")),
            )
            errMsg += "Masses of bank filter template: %e %e\n" % (
                getattr(bank_template.params, "mass1", float("nan")),
                getattr(bank_template.params, "mass2", float("nan")),
            )
            errMsg += "Overlap: %e" % (abs(overlaps[-1]))
            logger.info(errMsg)
    return overlaps


def bank_chisq_from_filters(
    tmplt_snr, tmplt_norm, bank_snrs, bank_norms, tmplt_bank_matches, indices=None
):
    """
    This function calculates and returns a TimeSeries object containing the
    bank veto calculated over a segment.

    Parameters
    ----------
    tmplt_snr: TimeSeries
        The SNR time series from filtering the segment against the current
        search template
    tmplt_norm: float
        The normalization factor for the search template
    bank_snrs: list of TimeSeries
        The precomputed list of SNR time series between each of the bank veto
        templates and the segment
    bank_norms: list of floats
        The normalization factors for the list of bank veto templates
        (usually this will be the same for all bank veto templates)
    tmplt_bank_matches: list of floats
        The complex overlap between the search template and each
        of the bank templates
    indices: {None, Array}, optional
        Array of indices into the snr time series. If given, the bank chisq
        will only be calculated at these values.

    Returns
    -------
    bank_chisq: TimeSeries of the bank vetos

    """
    if indices is not None:
        tmplt_snr = Array(tmplt_snr, copy=False)
        bank_snrs_tmp = []
        for bank_snr in bank_snrs:
            bank_snrs_tmp.append(bank_snr.take(indices))
        bank_snrs = bank_snrs_tmp

    # Initialise bank_chisq as 0s everywhere
    bank_chisq = zeros(len(tmplt_snr), dtype=real_same_precision_as(tmplt_snr))

    # Loop over all the bank templates
    for i in range(len(bank_snrs)):
        bank_match = tmplt_bank_matches[i]
        if abs(bank_match) > 0.99:
            # Not much point calculating bank_chisquared if the bank template
            # is very close to the filter template. Can also hit numerical
            # error due to approximations made in this calculation.
            # The value of 2 is the expected addition to the chisq for this
            # template
            bank_chisq += 2.0
            continue
        bank_norm = sqrt((1 - bank_match * bank_match.conj()).real)

        bank_SNR = bank_snrs[i] * (bank_norms[i] / bank_norm)
        tmplt_SNR = tmplt_snr * (bank_match.conj() * tmplt_norm / bank_norm)

        bank_SNR = Array(bank_SNR, copy=False)
        tmplt_SNR = Array(tmplt_SNR, copy=False)

        bank_chisq += (bank_SNR - tmplt_SNR).squared_norm()

    if indices is not None:
        return bank_chisq
    return TimeSeries(
        bank_chisq,
        delta_t=tmplt_snr.delta_t,
        epoch=tmplt_snr.start_time,
        copy=False,
    )


class SingleDetBankVeto:
    """
    This class reads in a template bank file for a bank veto, handles the
    memory management of its filters internally, and calculates the bank
    veto TimeSeries.
    """

    name = "bank"
    columns = {"bank_chisq": numpy.float32, "bank_chisq_dof": int}
    default_values = {"bank_chisq": 0, "bank_chisq_dof": 0}

    def __init__(
        self,
        bank_file,
        flen=None,
        delta_f=None,
        f_low=None,
        cdtype=None,
        approximant=None,
        **kwds,
    ):
        if bank_file is not None:
            self.do = True

            self.column_name = "bank_chisq"
            self.table_dof_name = "bank_chisq_dof"

            self.bank_file = bank_file
            self.cdtype = cdtype
            self.delta_f = delta_f
            self.f_low = f_low
            self.approximant = approximant
            self.kwds = kwds
            self.filters = None
            self.dof = 0
            if flen is not None:
                self._init_filters(flen, delta_f, f_low, cdtype)

            self._overlaps_cache = {}
            self._segment_snrs_cache = {}
        else:
            self.do = False

    def default_steps(self):
        return [self.name]

    def evaluate(self, batch, j):
        stilde = batch.stilde
        template = batch.template(j)
        chi, dof = self.values(
            template,
            stilde.psd,
            stilde,
            numpy.array([batch.snr[j]]),
            1,
            numpy.array([batch.time_index[j]]),
            template_whitened=True,
            whitened=getattr(stilde, "_is_whitened", False),
        )
        if chi is None:
            return None
        return {"bank_chisq": float(chi[0]), "bank_chisq_dof": int(dof[0])}

    def _init_filters(self, flen, delta_f, f_low, cdtype):
        self.cdtype = cdtype
        self.delta_f = delta_f
        self.f_low = f_low
        self.seg_len_freq = flen
        self.seg_len_time = (self.seg_len_freq - 1) * 2

        logger.info("Read in bank veto template bank")
        bank_veto_bank = FilterBank(
            self.bank_file,
            self.seg_len_freq,
            self.delta_f,
            self.cdtype,
            low_frequency_cutoff=f_low,
            approximant=self.approximant,
            **self.kwds,
        )

        self.filters = list(bank_veto_bank)
        self.dof = len(bank_veto_bank) * 2

    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Bank Veto")
        group.add_argument(
            "--bank-veto-bank-file",
            type=str,
            help="Template bank file used to compute the bank chi-square veto.",
        )

    @classmethod
    def from_cli(cls, args, flen=None, delta_f=None, f_low=None, cdtype=None, **kwargs):
        kwds = {}
        if getattr(args, "order", None) is not None:
            kwds["phase_order"] = args.order
        return cls(
            args.bank_veto_bank_file,
            flen=flen,
            delta_f=delta_f,
            f_low=f_low
            if f_low is not None
            else getattr(args, "low_frequency_cutoff", None),
            cdtype=cdtype,
            approximant=getattr(args, "approximant", None),
            **kwds,
        )

    def cache_segment_snrs(self, stilde, psd, whitened=True):
        key = (id(stilde), id(psd))
        if key not in self._segment_snrs_cache:
            logger.debug("Precalculate the bank veto template snrs")
            data = segment_snrs(
                self.filters, stilde, psd, self.f_low, whitened=whitened
            )
            self._segment_snrs_cache[key] = data
        return self._segment_snrs_cache[key]

    def cache_overlaps(self, template, psd, template_whitened=False):
        key = (id(template.params), id(psd))
        if key not in self._overlaps_cache:
            logger.debug("...Calculate bank veto overlaps")
            o = template_overlaps(
                self.filters,
                template,
                psd,
                self.f_low,
                template_whitened=template_whitened,
            )
            self._overlaps_cache[key] = o
        return self._overlaps_cache[key]

    def values(
        self,
        template,
        psd,
        stilde,
        snrv,
        norm,
        indices,
        template_whitened=False,
        whitened=True,
    ):
        """
        Returns
        -------
        bank_chisq_from_filters: TimeSeries of bank veto values - if indices
        is None then evaluated at all time samples, if not then only at
        requested sample indices

        bank_chisq_dof: int, approx number of statistical degrees of freedom

        """
        if not self.do:
            return None, None
        logger.debug("...Doing bank veto")
        if self.filters is None:
            self._init_filters(
                len(stilde),
                stilde.delta_f,
                self.f_low if self.f_low is not None else template.f_lower,
                stilde.dtype,
            )
        overlaps = self.cache_overlaps(
            template, psd, template_whitened=template_whitened
        )
        bank_veto_snrs, bank_veto_norms = self.cache_segment_snrs(
            stilde, psd, whitened=whitened
        )
        chisq = bank_chisq_from_filters(
            snrv, norm, bank_veto_snrs, bank_veto_norms, overlaps, indices
        )
        dof = numpy.repeat(self.dof, len(chisq))
        return chisq, dof


class SingleDetSkyMaxBankVeto(SingleDetBankVeto):
    """Stub for precessing bank veto if anyone ever wants to code it up."""

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

    def values(self, *args, **kwargs):
        if self.do:
            err_msg = "Precessing single detector sky-max bank veto has not "
            err_msg += "been written. If you want to use it, why not help "
            err_msg += "write it?"
            raise NotImplementedError(err_msg)
        return None, None
