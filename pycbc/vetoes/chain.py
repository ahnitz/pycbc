# Copyright (C) 2026 Alex Nitz
#
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
"""
Run signal-consistency tests on a batch of triggers as an ordered chain.

A chain is a list of steps, run in order. There are two kinds of step:

* A test: 'power', 'auto' or 'sg', computing the power, auto or
  sine-Gaussian chisq with the existing classes in pycbc.vetoes.
* A threshold, written PARAMETER:VALUE:LIMIT:ACTION, with PARAMETER and
  LIMIT as for --trigger-cuts. PARAMETER may be any single-detector
  ranking statistic in pycbc.events.ranking, a reduced chisq choice such
  as 'traditional', or a trigger column such as 'snr'. ACTION says what
  happens to the triggers that fail: with 'remove' they are dropped from
  the output at once, so nothing more is computed for them; with 'skip'
  they are kept in the output, but no later test is computed for them.
  A threshold only looks at the triggers still being tested: those not
  removed or skipped by an earlier threshold.

A test is computed for the loudest trigger first. Its value is then also
given to the other triggers of the same template within --chisq-window of
it, and the next loudest trigger without a value is computed, up to
--chisq-iterations times per test. If the window is no longer than the
clustering window, the trigger that clustering keeps is therefore always
one whose value was computed. A test is only run on triggers that every
earlier test gave a value; the others keep the test's default value.

The tests' own SNR thresholds are steps of the chain too: from_cli builds
the default chain, 'snr:T:lower:skip power auto newsnr:S:lower_inc:skip
sg', from --chisq-snr-threshold T and --sgchisq-snr-threshold S, as
pycbc_inspiral applies them.

SignalConsistencyChain.values is given a batch of triggers from one
segment of data, with two functions: one returning a trigger's template,
the other its template's correlation with the data around it.
"""

import logging
import time
from typing import ClassVar

import numpy as np

from pycbc.events import cuts, ranking
from pycbc.io.hdf import get_chisq_from_file_choice
from pycbc.types import FrequencySeries
from pycbc.vetoes.autochisq import SingleDetAutoChisq, autochisq_from_precomputed
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq

logger = logging.getLogger("pycbc.vetoes.chain")


class _Batch:
    """
    The arguments of SignalConsistencyChain.values, and what the tests
    keep for the batch.
    """

    def __init__(
        self, snr, time_index, template_id, template, correlation, stilde, columns
    ):
        self.snr = np.asarray(snr)
        self.time_index = np.asarray(time_index)
        self.template_id = np.asarray(template_id)
        self.template = template
        self.correlation = correlation
        self.stilde = stilde
        self.columns = columns or {}
        self.autocorrelations = {}


def _unit_psd(template):
    """
    A PSD of ones on the template's frequency grid: the templates are
    already whitened.
    """
    key = (len(template), template.delta_f)
    if key not in _unit_psds:
        _unit_psds[key] = FrequencySeries(
            np.ones(len(template), dtype=np.float32), delta_f=template.delta_f
        )
    return _unit_psds[key]


_unit_psds = {}


# Columns needed by each reduced chisq choice of pycbc.io.hdf
_chisq_choice_columns = {
    "traditional": ["chisq", "chisq_dof"],
    "cont": ["cont_chisq", "cont_chisq_dof"],
    "sg": ["sg_chisq"],
    "max_cont_trad": ["chisq", "chisq_dof", "cont_chisq", "cont_chisq_dof"],
}


class _Stage:
    """
    One test of the chain, computed by estimator. A subclass defines

    name : str
        The step name that selects it.
    columns : dict
        The columns it computes, mapped to their types.
    default(n)
        The value of each column for n triggers it gives no value.
    evaluate(batch, j)
        The value of each column for trigger j, or None if none can be
        computed for it.
    """

    def __init__(self, estimator):
        self.estimator = estimator


class _PowerStage(_Stage):
    name = "power"
    columns: ClassVar[dict] = {"chisq": np.float32, "chisq_dof": int}

    @staticmethod
    def default(n):
        return {
            "chisq": np.zeros(n, dtype=np.float32),
            "chisq_dof": np.zeros(n, dtype=int) - 100,
        }

    def evaluate(self, batch, j):
        corr, index = batch.correlation(j, 0)
        template = batch.template(j)
        chi, dof = self.estimator.values(
            corr,
            np.array([batch.snr[j]]),
            1,
            _unit_psd(template),
            np.array([index]),
            template,
        )
        if chi is None:
            return None
        return {"chisq": chi[0], "chisq_dof": dof[0]}


class _AutoStage(_Stage):
    name = "auto"
    columns: ClassVar[dict] = {"cont_chisq": np.float32}

    @staticmethod
    def default(n):
        return {"cont_chisq": np.zeros(n, dtype=np.float32)}

    def evaluate(self, batch, j):
        est = self.estimator
        corr, index = batch.correlation(j, est.stride * est.num_points)
        if corr is None:
            return None
        corr = np.asarray(corr)
        sn = np.fft.ifft(corr) * len(corr)
        _, achisq = autochisq_from_precomputed(
            sn,
            sn,
            self.autocorrelation(batch, j),
            np.array([index]),
            stride=est.stride,
            num_points=est.num_points,
            oneside=est.one_sided,
            twophase=est.two_phase,
            maxvalued=est.take_maximum_value,
        )
        return {"cont_chisq": achisq[0]}

    @staticmethod
    def autocorrelation(batch, j):
        """
        Return the autocorrelation of trigger j's template, normalized to
        1 at zero lag: the inverse FFT of the whitened template's power over
        the positive frequencies, as SingleDetAutoChisq computes it.
        """
        tid = batch.template_id[j]
        if tid not in batch.autocorrelations:
            tmplt = batch.template(j)
            autocor = getattr(tmplt, "autocorrelation", None)
            if autocor is None:
                power = np.abs(np.asarray(tmplt)) ** 2
                power[len(power) // 2 + 1 :] = 0
                autocor = np.fft.ifft(power)
                autocor = autocor / autocor[0]
                if hasattr(tmplt, "__dict__"):
                    tmplt.autocorrelation = autocor
            batch.autocorrelations[tid] = autocor
        return batch.autocorrelations[tid]


class _SGStage(_Stage):
    name = "sg"
    columns: ClassVar[dict] = {"sg_chisq": np.float32}

    @staticmethod
    def default(n):
        return {"sg_chisq": np.ones(n, dtype=np.float32)}

    def evaluate(self, batch, j):
        stilde = batch.stilde
        template = batch.template(j)
        # The trigger's time in seconds from the start of stilde
        seconds = batch.time_index[j] / (2 * (len(stilde) - 1) * stilde.delta_f)
        val = self.estimator.values_at(
            stilde,
            stilde.psd,
            template,
            _unit_psd(template),
            np.array([seconds]),
            whitened=getattr(stilde, "_is_whitened", False),
            snrs=np.array([batch.snr[j]]),
        )
        return {"sg_chisq": val[0]}


class SignalConsistencyChain:
    """
    Run the power, auto and sine-Gaussian chisq tests, with thresholds
    between them, in a configurable order (see the module documentation).

    Parameters
    ----------
    power_chisq : pycbc.vetoes.SingleDetPowerChisq, optional
        Computes the power chisq. Give it no snr_threshold: which triggers
        the test is computed for is set by the threshold steps.
    autochisq : pycbc.vetoes.SingleDetAutoChisq, optional
        Computes the auto chisq.
    sg_chisq : pycbc.vetoes.SingleDetSGChisq, optional
        Computes the sine-Gaussian chisq.
    steps : list, optional
        The steps of the chain, in order: test names ('power', 'auto',
        'sg') and threshold strings PARAMETER:VALUE:LIMIT:ACTION. By
        default, every test that is configured, in that order, with no
        thresholds.
    window : int
        A computed value is also given to the triggers of the same template
        within this many samples.
    max_iterations : int, optional
        The maximum number of values each test computes per batch.
    input_columns : list of str
        The columns, such as psd_var_val, that values is given in its
        columns argument, which threshold steps may use.
    defaults : dict, optional
        The value to give a column for triggers whose test is not
        computed, in place of the test's own default.

    """

    _ALL_STAGES = (_PowerStage, _AutoStage, _SGStage)
    columns: ClassVar[dict] = {
        col: dtype
        for stage in _ALL_STAGES
        for col, dtype in stage.columns.items()
    }

    def __init__(
        self,
        power_chisq=None,
        autochisq=None,
        sg_chisq=None,
        steps=None,
        window=0,
        max_iterations=None,
        input_columns=(),
        defaults=None,
    ):
        stages = {}
        for stage, est in (
            (_PowerStage, power_chisq),
            (_AutoStage, autochisq),
            (_SGStage, sg_chisq),
        ):
            if est is not None and est.do:
                stages[stage.name] = stage(est)

        if steps is None:
            steps = list(stages)
        self.description = " ".join(steps)

        self.window = window
        self.max_iterations = max_iterations
        self.defaults = defaults or {}
        self.stage_times = {"power": 0.0, "auto": 0.0, "sg": 0.0}

        available = {"snr"} | set(input_columns)
        self.steps = []
        for step in steps:
            if step in ("power", "auto", "sg"):
                if step not in stages:
                    raise ValueError(  # noqa: TRY003
                        f"Stage '{step}' is in the chain but its options are not set"
                    )
                stage = stages[step]
                self.steps.append(stage)
                available |= set(stage.columns)
                if step == "auto":
                    available.add("cont_chisq_dof")
                    self._auto_dof = stage.estimator.dof
            else:
                param, func, value, action = self._parse_threshold(step)
                missing = set(self._required(param)) - available
                if missing:
                    raise ValueError(  # noqa: TRY003
                        f"'{step}' needs {sorted(missing)}, "
                        "which no earlier step provides"
                    )
                self.steps.append((param, func, value, action, step))
        logger.info("Signal consistency chain: %s", self.description)

    @staticmethod
    def _parse_threshold(step):
        fields = step.split(":")
        if len(fields) != 4 or fields[3].lower() not in ("remove", "skip"):
            raise ValueError(  # noqa: TRY003
                "A threshold step is "
                "PARAMETER:VALUE:LIMIT:ACTION, with ACTION "
                "'remove' (drop the triggers that fail from the "
                "output) or 'skip' (keep them, but compute no "
                f"later test for them); got '{step}'"
            )
        action = fields.pop().lower()
        choices = (
            list(ranking.sngls_ranking_function_dict)
            + list(_chisq_choice_columns)
            + ["snr", "psd_var_val"]
        )
        (param, func), value = next(
            iter(cuts.convert_inputstr(":".join(fields), choices).items())
        )
        return param, func, value, action

    @staticmethod
    def _required(param):
        if param in ranking.sngls_ranking_function_dict:
            return ranking.reqd_datasets.get(param, [])
        return _chisq_choice_columns.get(param, [param])

    def values(
        self,
        snr,
        time_index,
        template_id,
        template,
        correlation,
        stilde=None,
        columns=None,
    ):
        """
        Run the chain on a batch of triggers from one segment of data.

        Parameters
        ----------
        snr : complex np.ndarray
            The normalized SNR of each trigger.
        time_index : np.ndarray
            The time of each trigger, in samples from the start of stilde.
        template_id : np.ndarray
            The template of each trigger. Values are only shared between
            triggers of the same template.
        template : function
            template(j) returns the template of trigger j, divided by the
            square root of the PSD (whitened), on the frequency grid of its
            correlation, with its f_lower and params attributes.
        correlation : function
            correlation(j, reach) returns (corr, index): corr is the
            frequency-domain correlation of trigger j's template with a
            block of data around it, normalized so that the SNR time series
            over the block is ifft(corr) * len(corr), and index is the
            trigger's sample in the block. The SNR series must be exact
            (unaffected by the block's edges) within reach samples of the
            trigger: the auto chisq reads it that far either side.
        stilde : pycbc.types.FrequencySeries, optional
            The segment of data, not divided by its PSD, with its psd
            attribute. Needed only for the sine-Gaussian chisq.
        columns : dict, optional
            Further columns of each trigger, such as psd_var_val, that
            threshold steps may use.

        Returns
        -------
        keep : np.ndarray of bool
            Which triggers are kept in the output: those no 'remove'
            threshold dropped.
        values : dict of np.ndarray
            Each column the tests compute, with one entry per trigger of the
            batch.

        """
        batch = _Batch(
            snr, time_index, template_id, template, correlation, stilde, columns
        )
        n = len(batch.snr)
        keep = np.ones(n, dtype=bool)
        # The triggers still being tested: thresholds and tests narrow it
        alive = np.ones(n, dtype=bool)
        vals = {}
        for stage in self._ALL_STAGES:
            for col, arr in stage.default(n).items():
                if col in self.defaults:
                    arr[:] = self.defaults[col]
                vals[col] = arr
        for step in self.steps:
            if isinstance(step, tuple):
                idx = np.flatnonzero(alive)
                failed = idx[~self._threshold(step, batch, vals, idx)]
                alive[failed] = False
                if step[3] == "remove":
                    keep[failed] = False
                logger.debug(
                    "%s: %d of %d triggers %s",
                    step[4],
                    len(failed),
                    len(idx),
                    "removed" if step[3] == "remove" else "skipped",
                )
            else:
                self._run_stage(step, batch, vals, alive)
        return keep, vals

    def _threshold(self, step, batch, vals, idx):
        param, func, value = step[:3]
        trigs = {"snr": np.abs(batch.snr[idx])}
        for col, arr in vals.items():
            trigs[col] = arr[idx]
        for col, arr in batch.columns.items():
            trigs[col] = np.asarray(arr)[idx]
        if "chisq_dof" in trigs:
            # EventManager's convention (2p - 2); ranking expects p bins
            trigs["chisq_dof"] = trigs["chisq_dof"] / 2 + 1
        if hasattr(self, "_auto_dof"):
            trigs["cont_chisq_dof"] = np.repeat(float(self._auto_dof), len(idx))

        if param in ranking.sngls_ranking_function_dict:
            stat = ranking.get_sngls_ranking_from_trigs(trigs, param)
        elif param in _chisq_choice_columns:
            stat = get_chisq_from_file_choice(trigs, param)
        else:
            stat = trigs[param]
        return func(stat, value)

    def _run_stage(self, stage, batch, vals, alive):
        if not np.any(alive):
            return
        start = time.perf_counter()
        priority = np.abs(batch.snr)
        times = batch.time_index
        tids = batch.template_id
        assigned = ~alive
        computed = 0
        for j in np.argsort(priority)[::-1]:
            if not alive[j]:
                continue
            if self.max_iterations is not None and computed >= self.max_iterations:
                break
            if assigned[j]:
                continue
            out = stage.evaluate(batch, j)
            if out is None:
                # Nothing could be computed for this trigger: the next one
                # in its window may stand in for it
                continue
            share = (
                (np.abs(times - times[j]) <= self.window)
                & (tids == tids[j])
                & ~assigned
            )
            for col, val in out.items():
                vals[col][share] = val
            assigned |= share
            computed += 1
        # Later stages only test triggers this one gave a value to
        alive &= assigned
        elapsed = time.perf_counter() - start
        self.stage_times[stage.name] = self.stage_times.get(stage.name, 0.0) + elapsed
        logger.debug(
            "%s chisq: %d evaluations in %.3f s",
            stage.name,
            computed,
            elapsed,
        )

    @staticmethod
    def insert_option_group(parser):
        """
        Add the options of the chain and of each of its tests to an
        argument parser.
        """
        group = parser.add_argument_group("Signal consistency tests")
        group.add_argument(
            "--chisq-window",
            type=float,
            default=1.0,
            help="A test computed for a trigger is also given "
            "to the triggers of the same template within "
            "this many seconds of it.",
        )
        group.add_argument(
            "--chisq-iterations",
            type=int,
            default=10000,
            help="The maximum number of values each test "
            "computes per batch, loudest trigger first. "
            "Triggers left without one keep the test's "
            "default value.",
        )
        group.add_argument(
            "--signal-consistency-chain",
            nargs="+",
            metavar="STEP",
            help="The tests to run and the thresholds between"
            " them, in order. Each STEP is a test "
            "('power', 'auto', 'sg') or a threshold "
            "PARAMETER:VALUE:LIMIT:ACTION, with PARAMETER"
            " and LIMIT as for --trigger-cuts and ACTION "
            "'remove' (drop the triggers that fail from "
            "the output) or 'skip' (keep them, but "
            "compute no later test for them). Ex. "
            "'snr:6:lower:skip power "
            "newsnr:5:lower:remove sg'. By default: "
            "'snr:T:lower:skip power auto "
            "newsnr:S:lower_inc:skip sg', with T the "
            "--chisq-snr-threshold and S the "
            "--sgchisq-snr-threshold, leaving out the "
            "tests that are not configured.",
        )
        SingleDetPowerChisq.insert_option_group(parser)
        SingleDetAutoChisq.insert_option_group(parser)
        SingleDetSGChisq.insert_option_group(parser)

    @classmethod
    def from_cli(
        cls,
        args,
        bank,
        sample_rate,
        power_chisq=None,
        sg_chisq=None,
        defaults=None,
    ):
        """
        Create the chain from the options added by insert_option_group.

        Parameters
        ----------
        args : argparse.Namespace
            The parsed options.
        bank : pycbc.waveform.bank.TemplateBank
            The bank, to select the templates the sine-Gaussian chisq is
            applied to.
        sample_rate : float
            The sample rate of the SNR time series, to convert
            --chisq-window to samples.
        power_chisq : SingleDetPowerChisq, optional
            Custom power-chisq estimator (defaults to SingleDetPowerChisq).
        sg_chisq : SingleDetSGChisq, optional
            Pre-built sine-Gaussian chisq estimator.
        defaults : dict, optional
            Default values for columns of uncomputed triggers.

        """
        if args.autochi_reverse_template:
            raise ValueError("The chain does not support --autochi-reverse-template")  # noqa: TRY003
        # The tests' SNR thresholds are steps of the chain, not of the
        # estimators
        power = (
            SingleDetPowerChisq(args.chisq_bins)
            if power_chisq is None
            else power_chisq
        )
        auto = SingleDetAutoChisq.from_cli(args)
        sg = (
            SingleDetSGChisq.from_cli(args, bank, str(args.chisq_bins))
            if sg_chisq is None
            else sg_chisq
        )
        extra = ["psd_var_val"] if getattr(args, "psdvar_segment", None) else []
        steps = args.signal_consistency_chain
        if steps is None:
            steps = []
            if power.do:
                if args.chisq_snr_threshold:
                    steps.append(f"snr:{args.chisq_snr_threshold!r}:lower:skip")
                steps.append("power")
            if auto.do:
                steps.append("auto")
            if sg.do:
                steps += [f"newsnr:{sg.snr_threshold!r}:lower_inc:skip", "sg"]
        return cls(
            power,
            auto,
            sg,
            steps=steps,
            window=int(args.chisq_window * sample_rate),
            max_iterations=args.chisq_iterations,
            input_columns=extra,
            defaults=defaults,
        )
