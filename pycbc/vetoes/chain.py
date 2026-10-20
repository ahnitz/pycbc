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

* A test: 'power', 'bank', 'auto' or 'sg', computing the power, bank,
  auto or sine-Gaussian chisq with the existing classes in pycbc.vetoes.
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
the default chain from each configured test's default_steps().

SignalConsistencyChain.values is given a batch of triggers from one
segment of data, with two functions: one returning a trigger's template,
the other its template's correlation with the data around it.
"""

import logging
import time
import types

import numpy as np

from pycbc.events import cuts, ranking
from pycbc.io.hdf import chisq_choices, get_chisq_from_file_choice
from pycbc.vetoes.autochisq import SingleDetAutoChisq
from pycbc.vetoes.bank_chisq import SingleDetBankVeto
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq

logger = logging.getLogger("pycbc.vetoes.chain")


class _ColumnProbe(dict):
    """Records which keys get_chisq_from_file_choice reads for a choice."""

    def __missing__(self, key):
        self[key] = np.full(1, 2.0)
        return self[key]


class SignalConsistencyChain:
    """
    Run signal-consistency tests, with thresholds between them, in a
    configurable order (see the module documentation).

    Parameters
    ----------
    *estimators
        Configured test estimators (or passed by keyword as power_chisq,
        bank_chisq, autochisq, sg_chisq).
    steps : list, optional
        The steps of the chain, in order: test names and threshold strings
        PARAMETER:VALUE:LIMIT:ACTION. By default, every test that is
        configured, in order, with no thresholds.
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

    TESTS = (
        SingleDetPowerChisq,
        SingleDetBankVeto,
        SingleDetAutoChisq,
        SingleDetSGChisq,
    )
    columns = {col: dtype for cls in TESTS for col, dtype in cls.columns.items()}

    def __init__(
        self,
        *estimators,
        steps=None,
        window=0,
        max_iterations=None,
        input_columns=(),
        defaults=None,
        **named_estimators,
    ):
        self.stages = {
            est.name: est
            for est in (*estimators, *named_estimators.values())
            if est is not None and est.do
        }

        if steps is None:
            steps = list(self.stages)
        self.description = " ".join(steps)

        self.window = window
        self.max_iterations = max_iterations
        self.defaults = defaults or {}
        self.stage_times = {cls.name: 0.0 for cls in self.TESTS}

        available = {"snr"} | set(input_columns)
        self.steps = []
        for step in steps:
            if ":" not in step:
                if step not in self.stages:
                    raise ValueError(
                        f"Stage '{step}' is in the chain but its options are not set"
                    )
                stage = self.stages[step]
                self.steps.append(stage)
                available |= set(stage.columns)
            else:
                param, func, value, action = self._parse_threshold(step)
                missing = set(self._required(param)) - available
                if missing:
                    raise ValueError(
                        f"'{step}' needs {sorted(missing)}, "
                        "which no earlier step provides"
                    )
                self.steps.append((param, func, value, action, step))
        logger.info("Signal consistency chain: %s", self.description)

    @staticmethod
    def _parse_threshold(step):
        fields = step.split(":")
        if len(fields) != 4 or fields[3].lower() not in ("remove", "skip"):
            raise ValueError(
                "A threshold step is "
                "PARAMETER:VALUE:LIMIT:ACTION, with ACTION "
                "'remove' (drop the triggers that fail from the "
                "output) or 'skip' (keep them, but compute no "
                f"later test for them); got '{step}'"
            )
        action = fields.pop().lower()
        choices = (
            list(ranking.sngls_ranking_function_dict)
            + list(chisq_choices)
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
        if param in chisq_choices:
            probe = _ColumnProbe()
            get_chisq_from_file_choice(probe, param)
            return list(probe)
        return [param]

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
            The segment of data, with its psd attribute.
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
        batch = types.SimpleNamespace(
            snr=np.asarray(snr),
            time_index=np.asarray(time_index),
            template_id=np.asarray(template_id),
            template=template,
            correlation=correlation,
            stilde=stilde,
            columns=columns or {},
        )
        n = len(batch.snr)
        keep = np.ones(n, dtype=bool)
        # The triggers still being tested: thresholds and tests narrow it
        alive = np.ones(n, dtype=bool)
        vals = {
            col: np.full(
                n,
                self.defaults.get(col, src.default_values[col]),
                dtype=dtype,
            )
            for src in (*self.TESTS, *self.stages.values())
            for col, dtype in src.columns.items()
        }
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

        if param in ranking.sngls_ranking_function_dict:
            stat = ranking.get_sngls_ranking_from_trigs(trigs, param)
        elif param in chisq_choices:
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

    @classmethod
    def insert_option_group(cls, parser):
        """
        Add the options of the chain and of each of its tests to an
        argument parser.
        """
        names = ", ".join(f"'{t.name}'" for t in cls.TESTS)
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
            f" them, in order. Each STEP is a test ({names})"
            " or a threshold PARAMETER:VALUE:LIMIT:ACTION, "
            "with PARAMETER and LIMIT as for --trigger-cuts "
            "and ACTION 'remove' (drop the triggers that "
            "fail from the output) or 'skip' (keep them, "
            "but compute no later test for them). Ex. "
            "'snr:6:lower:skip power "
            "newsnr:5:lower:remove sg'.",
        )
        for test_cls in cls.TESTS:
            test_cls.insert_option_group(parser)

    @classmethod
    def from_cli(cls, args, bank, sample_rate, defaults=None):
        """
        Create the chain from the options added by insert_option_group.

        Parameters
        ----------
        args : argparse.Namespace
            The parsed options.
        bank : pycbc.waveform.bank.TemplateBank
            The template bank.
        sample_rate : float
            The sample rate of the SNR time series, to convert
            --chisq-window to samples.
        defaults : dict, optional
            Default values for columns of uncomputed triggers.

        """
        tests = [test_cls.from_cli(args, bank=bank) for test_cls in cls.TESTS]
        steps = args.signal_consistency_chain
        if steps is None:
            steps = [step for test in tests if test.do for step in test.default_steps()]
        # The tests' SNR thresholds are steps of the chain, not of the
        # estimators
        for test in tests:
            if hasattr(test, "snr_threshold"):
                test.snr_threshold = None
        extra = ["psd_var_val"] if getattr(args, "psdvar_segment", None) else []
        return cls(
            *tests,
            steps=steps,
            window=int(args.chisq_window * sample_rate),
            max_iterations=args.chisq_iterations,
            input_columns=extra,
            defaults=defaults,
        )
