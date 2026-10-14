"""Tests of the sine-Gaussian chisq with the template's response removed."""

import unittest

import numpy

from pycbc import DYN_RANGE_FAC
from pycbc.filter import matched_filter_core, sigma
from pycbc.io import WaveformArray
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.vetoes.sgchisq import SingleDetSGChisq
from pycbc.waveform import get_fd_waveform

# Tiles centred 0, 5 and 10 Hz above the peak frequency (151 Hz) of a
# 60+60 solar mass template, where each overlaps it: |<tile|template>| is
# 0.19 to 0.27
TILES = "20-0,20-5,20-10"


class Bank:
    # A bank of the one template, to select the templates the tiles are
    # applied to
    table = WaveformArray.from_kwargs(template_hash=[1], mass1=[60.0], mass2=[60.0])


class TestSGChisqProjection(unittest.TestCase):
    def setUp(self):
        self.sr, duration = 2048, 16
        self.n = self.sr * duration
        df = 1.0 / duration
        psd = aLIGOZeroDetHighPower(self.n // 2 + 1, df, 15)
        # Scaled by DYN_RANGE_FAC^2, as pycbc_inspiral does, so that the
        # single precision PSD does not underflow
        self.psd64 = psd * DYN_RANGE_FAC**2
        self.psd = self.psd64.astype(numpy.float32)
        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD", mass1=60, mass2=60, f_lower=20, delta_f=df
        )
        hp.resize(self.n // 2 + 1)
        self.template = hp.astype(numpy.complex64)
        # As a bank template: its params select its tiles, and its f_lower
        # is where its power starts
        self.template.params = Bank.table[0]
        self.template.f_lower = 20.0
        # The chisq bins place the template's peak frequency
        self.sg = SingleDetSGChisq(
            Bank, "16", snr_threshold=0, chisq_locations=["mtotal>0:" + TILES]
        )

    def sg_chisq(self, seed, amp=0, times=(8.0,), project=True):
        """
        Return the sine-Gaussian chisq at the given times (in seconds,
        on the sample grid) of 16 s of Gaussian noise. If amp is nonzero, a
        signal matching the template is added at 8 s with SNR amp. If
        project is true, the template's response is removed from the
        tiles.
        """
        data = noise_from_psd(self.n, 1.0 / self.sr, self.psd64, seed=seed)
        data = data.to_frequencyseries().astype(numpy.complex64)
        if amp:
            norm = amp / sigma(self.template, psd=self.psd, low_frequency_cutoff=20)
            data = data + (self.template * norm).cyclic_time_shift(8.0)
        snr, _, snr_norm = matched_filter_core(
            self.template, data, psd=self.psd, low_frequency_cutoff=20
        )
        index = numpy.rint(numpy.array(times) * self.sr).astype(int)
        kw = {}
        if project:
            kw = dict(snrs=snr.numpy()[index] * snr_norm)
        # The data divided by the PSD, where the PSD is defined (the model
        # is zero below 15 Hz and at Nyquist)
        whitened = data.copy()
        defined = self.psd.numpy() > 0
        whitened.data[defined] /= self.psd.numpy()[defined]
        return self.sg.values_at(
            whitened, self.psd, self.template, self.psd, numpy.array(times), **kw
        )

    def test_template_response_is_removed(self):
        # With the template's response removed, adding a signal that
        # matches the template should leave the chisq as it is for the noise
        # alone. Without the removal, the SNR 200 signal adds about 1200 to
        # the chisq; with it, the change must be under 1e-3 of that (about
        # 1e-5 is measured).
        noise = self.sg_chisq(0)[0]
        signal = self.sg_chisq(0, amp=200)[0]
        unremoved = (
            self.sg_chisq(0, amp=200, project=False)[0]
            - self.sg_chisq(0, project=False)[0]
        )
        self.assertGreater(unremoved, 100)
        self.assertLess(abs(signal - noise), 1e-3 * unremoved)

    def test_noise_is_normalized(self):
        # Removing the template's response also removes a fraction
        # |<tile|template>|^2 of each tile's noise power, which the chisq
        # corrects for, so its mean in noise should stay 1. The mean of
        # these 20000 values has a standard error of about 0.007 (measured),
        # so it should lie within 0.025 of 1. Without the correction it
        # would be lower by the mean |<tile|template>|^2, 0.054.
        rng = numpy.random.default_rng(0)
        vals = numpy.concatenate(
            [
                self.sg_chisq(
                    s, times=rng.integers(2 * self.sr, 14 * self.sr, 50) / self.sr
                )
                for s in range(400)
            ]
        )
        self.assertLess(abs(numpy.mean(vals) - 1), 0.025)


if __name__ == "__main__":
    unittest.main()
