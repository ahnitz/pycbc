"""Unit tests for pycbc.vetoes.chain.SignalConsistencyChain."""

import unittest

import numpy

from pycbc.events import ranking
from pycbc.types import FrequencySeries
from pycbc.vetoes.chain import SignalConsistencyChain
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq


class StubPower(SingleDetPowerChisq):
    """
    A stand-in for SingleDetPowerChisq whose chisq is 10 |SNR|, with 6
    dof (4 bins). It counts how often it is called.
    """

    def __init__(self):
        super().__init__(num_bins=4)
        self.calls = 0

    def values(self, corr, snrv, snr_norm, psd, indices, template):
        self.calls += 1
        return numpy.array([10 * abs(snrv[0])]), numpy.array([2 * 4 - 2])


# The stand-in tests read no data: a placeholder template and correlation
PLACEHOLDER = FrequencySeries(numpy.ones(4), delta_f=1.0)


def run(
    chain,
    snr,
    times,
    tids=None,
    template=lambda j: PLACEHOLDER,
    correlation=lambda j, reach: (None, 0),
    **kwargs,
):
    """
    Run the chain on triggers of the given SNRs and times (of one
    template unless tids is given).
    """
    tids = numpy.zeros(len(snr), dtype=int) if tids is None else tids
    return chain.values(
        numpy.array(snr, dtype=complex),
        numpy.array(times),
        numpy.array(tids),
        template,
        correlation,
        **kwargs,
    )


class TestChain(unittest.TestCase):
    def test_loudest_first_and_window(self):
        power = StubPower()
        chain = SignalConsistencyChain(power, steps=["power"], window=10)
        keep, vals = run(chain, [5, 9, 6, 7], [0, 5, 100, 200])
        self.assertTrue(keep.all())
        # The loudest trigger (SNR 9 at t=5) is computed first, and its value
        # is also given to the trigger at t=0, within the window of 10
        numpy.testing.assert_allclose(vals["chisq"], [90, 90, 60, 70])
        self.assertEqual(power.calls, 3)

    def test_skip(self):
        # Triggers below SNR 6 fail the 'skip' threshold: they are kept, but
        # their chisq is not computed and keeps its default (0, dof -100)
        chain = SignalConsistencyChain(
            StubPower(), steps=["snr:6:lower_inc:skip", "power"], window=0
        )
        _, vals = run(chain, [5, 9, 7], [0, 50, 100])
        numpy.testing.assert_allclose(vals["chisq"], [0, 90, 70])
        numpy.testing.assert_array_equal(vals["chisq_dof"], [-100, 6, 6])

    def test_iteration_cap(self):
        # Each test computes at most max_iterations values, loudest trigger
        # first
        chain = SignalConsistencyChain(
            StubPower(), steps=["power"], window=0, max_iterations=1
        )
        _, vals = run(chain, [5, 9, 7], [0, 50, 100])
        numpy.testing.assert_allclose(vals["chisq"], [0, 90, 0])

    def test_window_is_per_template(self):
        # Two triggers within the window, but of different templates: each
        # is computed, rather than one sharing the other's value
        chain = SignalConsistencyChain(StubPower(), steps=["power"], window=10)
        _, vals = run(chain, [9, 5], [0, 5], tids=[0, 1])
        numpy.testing.assert_allclose(vals["chisq"], [90, 50])

    def test_remove_and_skip(self):
        power = StubPower()
        chain = SignalConsistencyChain(
            power, steps=["snr:4:lower:remove", "snr:6:lower:skip", "power"]
        )
        keep, vals = run(chain, [5, 9, 3], [0, 50, 100])
        # The trigger at SNR 3 is removed from the output; the one at SNR 5
        # is kept, but its chisq is not computed
        numpy.testing.assert_array_equal(keep, [True, True, False])
        numpy.testing.assert_allclose(vals["chisq"], [0, 90, 0])
        self.assertEqual(power.calls, 1)

    def test_threshold_sees_only_triggers_still_tested(self):
        # The final threshold only looks at the trigger still being tested
        # (SNR 9), so the triggers skipped earlier stay in the output
        chain = SignalConsistencyChain(
            StubPower(), steps=["snr:6:lower:skip", "power", "snr:4:lower:remove"]
        )
        keep, _ = run(chain, [5, 9, 3], [0, 50, 100])
        numpy.testing.assert_array_equal(keep, [True, True, True])

    def test_cut_on_ranking_statistic(self):
        # The stub's chisq is 10 |SNR| with 6 dof, so the reduced chisq is
        # 10 |SNR| / 6, and the NewSNRs of triggers at SNR 2.5 and 9 are 1.37
        # and 2.61. A cut at NewSNR 2.4 removes only the first, where a cut
        # on SNR alone would remove neither.
        chain = SignalConsistencyChain(
            StubPower(), steps=["power", "newsnr:2.4:lower:remove"]
        )
        keep, _ = run(chain, [2.5, 9.0], [0, 50])
        numpy.testing.assert_array_equal(keep, [False, True])

    def test_order_is_validated(self):
        # A threshold on NewSNR before the power chisq it needs, a test whose
        # estimator is not configured, and a threshold without an action
        # are all errors
        with self.assertRaises(ValueError):
            SignalConsistencyChain(
                StubPower(), steps=["newsnr:5:lower:remove", "power"]
            )
        with self.assertRaises(ValueError):
            SignalConsistencyChain(StubPower(), steps=["snr:5:lower", "power"])
        with self.assertRaises(ValueError):
            SignalConsistencyChain(StubPower(), steps=["power", "sg"])

    def test_new_ranking_statistic_is_usable(self):
        # A statistic added to pycbc.events.ranking, with the columns it
        # needs, can be used in a threshold without changing the chain
        def get_double_snr(trigs, **kwargs):
            return 2 * numpy.array(trigs["snr"][:])

        ranking.sngls_ranking_function_dict["double_snr"] = get_double_snr
        ranking.reqd_datasets["double_snr"] = ["snr"]
        try:
            chain = SignalConsistencyChain(
                StubPower(), steps=["double_snr:10:lower:remove", "power"]
            )
            keep, _ = run(chain, [4, 6], [0, 50])
            numpy.testing.assert_array_equal(keep, [False, True])
        finally:
            del ranking.sngls_ranking_function_dict["double_snr"]
            del ranking.reqd_datasets["double_snr"]

    def test_defaults_override(self):
        # defaults replaces the power chisq's own default (chisq 0, dof
        # -100) for the triggers it does not compute
        chain = SignalConsistencyChain(
            StubPower(),
            steps=["snr:6:lower_inc:skip", "power"],
            defaults={"chisq": 126, "chisq_dof": 126},
        )
        _, vals = run(chain, [5, 9], [0, 50])
        numpy.testing.assert_allclose(vals["chisq"], [126, 90])
        numpy.testing.assert_array_equal(vals["chisq_dof"], [126, 6])

    def test_skipped_trigger_lets_next_stand_in(self):
        class SkipLoudest(StubPower):
            def values(self, corr, snrv, snr_norm, psd, indices, template):
                if abs(snrv[0]) > 8:
                    return None, None
                return StubPower.values(
                    self, corr, snrv, snr_norm, psd, indices, template
                )

        chain = SignalConsistencyChain(SkipLoudest(), steps=["power"], window=10)
        _, vals = run(chain, [9, 7], [0, 5])
        # No value can be computed for the loudest trigger, so the next
        # loudest is computed instead and its value is given to both
        numpy.testing.assert_allclose(vals["chisq"], [70, 70])

    def test_cut_on_supplied_chisq(self):
        # Here the chisq is not computed by the chain but given in the
        # columns argument, as for a chain run on triggers whose power chisq
        # is already known. With dof 6 (4 bins), the reduced chisqs are 1 and
        # 3, so a cut at 2 removes the second.
        chain = SignalConsistencyChain(
            StubPower(),
            steps=["traditional:2:upper:remove", "power"],
            input_columns=["chisq", "chisq_dof"],
        )
        columns = {"chisq": numpy.array([6.0, 18.0]), "chisq_dof": numpy.array([6, 6])}
        keep, _ = run(chain, [5, 9], [0, 50], columns=columns)
        numpy.testing.assert_array_equal(keep, [True, False])

    def test_later_stages_need_a_value(self):
        class StubSG(SingleDetSGChisq):
            do = True
            calls = 0

            def __init__(self):
                pass

            def values(self, *args, **kwargs):
                self.calls += 1
                return numpy.array([2.0])

        class Data:
            # The sine-Gaussian stage computes the trigger's time from the
            # data's length and frequency step
            delta_f = 1.0
            psd = None

            def __len__(self):
                return 9

        sg = StubSG()
        chain = SignalConsistencyChain(
            StubPower(), sg_chisq=sg, steps=["power", "sg"], window=0, max_iterations=1
        )
        _, vals = run(chain, [5, 9, 7], [0, 50, 100], stilde=Data())
        # Only the trigger given a power chisq (SNR 9, as the cap allows one)
        # has its sine-Gaussian chisq computed; the others keep its default
        # of 1
        self.assertEqual(sg.calls, 1)
        numpy.testing.assert_allclose(vals["sg_chisq"], [1, 2, 1])

    def test_auto_chisq_matches_single_det_auto_chisq(self):
        # For a full-length SNR series, the chain's auto chisq equals
        # SingleDetAutoChisq.values for the same template and data
        from pycbc.filter import matched_filter, sigma
        from pycbc.types import FrequencySeries, TimeSeries
        from pycbc.vetoes.autochisq import SingleDetAutoChisq
        from pycbc.waveform import get_fd_waveform

        rng = numpy.random.default_rng(1)
        n, sr = 8192, 1024.0
        df = sr / n
        hp, _ = get_fd_waveform(
            approximant="TaylorF2", mass1=10, mass2=10, f_lower=30, delta_f=df
        )
        hp.resize(n // 2 + 1)
        psd = FrequencySeries(numpy.ones(n // 2 + 1), delta_f=df)
        # White noise for a one-sided PSD of ones has variance sr / 2
        noise = TimeSeries(rng.normal(size=n) * (sr / 2) ** 0.5, delta_t=1 / sr)
        signal = (
            hp.cyclic_time_shift(2.0).to_timeseries()
            * 20
            / sigma(hp, psd=psd, low_frequency_cutoff=30)
        )
        data = TimeSeries(
            noise.numpy() + signal.numpy(), delta_t=1 / sr
        ).to_frequencyseries()
        snr = matched_filter(hp, data, psd=psd, low_frequency_cutoff=30)
        index = int(numpy.argmax(abs(snr)))
        ref = SingleDetAutoChisq(2, 10, twophase=True)
        expect, dof = ref.values(
            snr, numpy.array([index]), hp, psd, 1.0, low_frequency_cutoff=30
        )

        # The SNR series is given as the correlation, normalized so that
        # ifft(corr) * len(corr) is the SNR, with the template (whitened, as
        # the PSD is ones) on the same full-length frequency grid
        corr = numpy.fft.fft(snr.numpy()) / n
        tmp = numpy.zeros(n, dtype=complex)
        k = numpy.arange(n // 2 + 1)
        # Zero below 30 Hz, the cutoff SingleDetAutoChisq is given
        tmp[k] = hp.numpy() * (k * df >= 30)
        auto = SingleDetAutoChisq(2, 10, twophase=True)
        chain = SignalConsistencyChain(autochisq=auto, steps=["auto"])
        _, vals = run(
            chain,
            [abs(snr[index])],
            [index],
            template=lambda j: tmp,
            correlation=lambda j, reach: (corr, index),
        )
        # cont_chisq is stored in single precision
        numpy.testing.assert_allclose(vals["cont_chisq"], expect, rtol=1e-6)
        self.assertEqual(auto.dof, dof)

    def test_default_chain_from_options(self):
        # The tests' SNR thresholds become the default chain's skip steps
        import argparse

        parser = argparse.ArgumentParser()
        SignalConsistencyChain.insert_option_group(parser)
        args = parser.parse_args(["--chisq-bins", "16", "--chisq-snr-threshold", "6"])
        chain = SignalConsistencyChain.from_cli(args, None, 1024.0)
        self.assertEqual(chain.description, "snr:6.0:lower:skip power")

    def test_bank_chisq_matches_single_det_bank_veto(self):
        import types

        from pycbc.filter import matched_filter_core, sigma
        from pycbc.types import FrequencySeries, TimeSeries
        from pycbc.vetoes.bank_chisq import SingleDetBankVeto
        from pycbc.waveform import get_fd_waveform

        rng = numpy.random.default_rng(2)
        n, sr = 8192, 1024.0
        df = sr / n
        hp, _ = get_fd_waveform(
            approximant="TaylorF2", mass1=10, mass2=10, f_lower=30, delta_f=df
        )
        hp.resize(n // 2 + 1)
        hp.f_lower = 30.0
        hp.params = types.SimpleNamespace(mass1=10.0, mass2=10.0)

        bank_tmplt, _ = get_fd_waveform(
            approximant="TaylorF2", mass1=12, mass2=12, f_lower=30, delta_f=df
        )
        bank_tmplt.resize(n // 2 + 1)
        bank_tmplt.f_lower = 30.0
        bank_tmplt.params = types.SimpleNamespace(mass1=12.0, mass2=12.0)

        psd = FrequencySeries(numpy.ones(n // 2 + 1), delta_f=df)
        noise = TimeSeries(rng.normal(size=n) * (sr / 2) ** 0.5, delta_t=1 / sr)
        signal = (
            hp.cyclic_time_shift(2.0).to_timeseries()
            * 20
            / sigma(hp, psd=psd, low_frequency_cutoff=30)
        )
        stilde = TimeSeries(
            noise.numpy() + signal.numpy(), delta_t=1 / sr
        ).to_frequencyseries()
        stilde.psd = psd

        snrv, _, norm = matched_filter_core(
            hp, stilde, psd=psd, low_frequency_cutoff=30
        )
        index = int(numpy.argmax(abs(snrv)))

        bank_veto = SingleDetBankVeto(None)
        bank_veto.do = True
        bank_veto.f_low = 30.0
        bank_veto.filters = [bank_tmplt]
        bank_veto.dof = 2
        bank_veto._overlaps_cache = {}
        bank_veto._segment_snrs_cache = {}

        expect, expect_dof = bank_veto.values(
            hp,
            psd,
            stilde,
            numpy.array([snrv[index]]),
            norm,
            numpy.array([index]),
            whitened=False,
        )

        chain = SignalConsistencyChain(
            bank_chisq=bank_veto, steps=["bank", "bank:50:upper:remove"]
        )
        keep, vals = run(
            chain,
            [snrv[index] * norm],
            [index],
            template=lambda j: hp,
            stilde=stilde,
        )
        self.assertTrue(keep.all())
        numpy.testing.assert_allclose(vals["bank_chisq"], expect, rtol=1e-5)
        numpy.testing.assert_array_equal(vals["bank_chisq_dof"], expect_dof)


if __name__ == "__main__":
    unittest.main()
