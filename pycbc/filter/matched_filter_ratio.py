import collections
import logging
import multiprocessing as mp
import time

import numpy as np

from pycbc.types import complex64, float32

# matchedfilter back end is required for ratio/FIR filtering.
try:
    import matchedfilter as _mf
    from matchedfilter import search as mf_search
    from matchedfilter.time_domain import reference_profile
except ImportError:
    _mf = None
    mf_search = None
    reference_profile = None


def _require_matchedfilter():
    if _mf is None:
        raise ImportError(
            "matchedfilter is required for ratio/FIR filtering but is not installed"
        )


def decimation_factor(tap_sample_rate, data_sample_rate):
    """Integer ratio of the bank's tap sample rate to data's; raises if not exact."""
    exact = float(tap_sample_rate) / float(data_sample_rate)
    factor = int(np.round(exact))
    if abs(exact - factor) > 1e-5 or factor < 1:
        raise ValueError(
            f"The bank sample rate ({tap_sample_rate} Hz) must be an exact integer "
            f"multiple of the data sample rate ({data_sample_rate} Hz); "
            f"the ratio is {exact:.4f}."
        )
    return factor


def analyzed_samples(stilde):
    """Samples searched for primary triggers: the windows if any, else analyze."""
    wins = getattr(stilde, "analyze_windows", None)
    if wins is not None:
        return sum(b - a for a, b in wins)
    return stilde.analyze.stop - stilde.analyze.start


def segments_with_asym_followup(segments, instruments):
    """Return {seg_idx: bool} indicating whether all active detectors overlap in
    each segment.
    """
    if len(instruments) <= 1:
        return collections.defaultdict(lambda: False)
    num_segs = max(len(segments[ifo]) for ifo in instruments)
    can_asym = {}
    for s_num in range(num_segs):
        active = [ifo for ifo in instruments if s_num < len(segments[ifo])]
        if len(active) <= 1:
            can_asym[s_num] = False
            continue
        bounds = []
        for ifo in active:
            st = segments[ifo][s_num]
            t0 = float(st.start_time)
            dt = float(st.delta_t)
            a0 = st.analyze.start if getattr(st, "analyze", None) is not None else 0
            a1 = (
                st.analyze.stop
                if getattr(st, "analyze", None) is not None
                else (len(st) - 1) * 2
            )
            bounds.append((t0 + a0 * dt, t0 + a1 * dt))
        can_asym[s_num] = max(b[0] for b in bounds) < min(b[1] for b in bounds)
    return can_asym


def sample_asym_bins(
    refilter,
    series_length,
    local_id,
    t_center,
    sample_rate=2048.0,
    bin_width=0.1,
    num_bins=100,
    valid_slice=None,
):
    """The peak SNR in each bin of an asymmetric window around t_center, as
    (template ids, samples, complex SNRs).
    """
    bs = max(1, int(round(bin_width * float(sample_rate or 2048.0))))
    win_start_samp = (t_center - bs // 2) - (num_bins // 2) * bs
    win_end_samp = win_start_samp + num_bins * bs

    v_start = 0 if valid_slice is None else valid_slice.start
    v_stop = series_length if valid_slice is None else valid_slice.stop

    asym_start = max(win_start_samp, v_start)
    asym_stop = min(win_end_samp, v_stop)

    if asym_start >= asym_stop:
        return (
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.complex64),
        )

    res = refilter(
        windows=slice(asym_start, asym_stop),
        binsize=bs,
        threshold=0.0,
        template_index=local_id,
    )
    return (
        res.template_indices.astype(np.int32, copy=False),
        res.sample_indices.astype(np.int64, copy=False),
        res.snr.astype(np.complex64, copy=False),
    )


def compute_followup_bin_budget(
    primary_snr,
    num_bins=1000,
    bin_width=0.03,
    min_bins=5,
    scale_factor=1.2,
    target_network_snr=8.5,
    secondary_snr_floor=3.0,
    f_eff=150.0,
    max_bins=None,
):
    """Compute symmetric followup bin budget K(rho_prim) from recovery difficulty."""
    rho_p = float(abs(primary_snr))
    rho_sec_req = np.sqrt(
        max(
            float(secondary_snr_floor) ** 2,
            float(target_network_snr) ** 2 - rho_p**2,
        )
    )

    window_duration = float(num_bins) * float(bin_width)
    n_trials = window_duration * float(f_eff)

    p_exceed = np.exp(-0.5 * rho_sec_req**2)
    mu_noise = n_trials * p_exceed

    k_analytic = np.ceil(mu_noise + 2.0 * np.sqrt(mu_noise) + 1.0)
    k_scaled = int(np.ceil(float(scale_factor) * k_analytic))

    upper_cap = int(num_bins)
    if max_bins is not None and int(max_bins) > 0:
        upper_cap = min(upper_cap, int(max_bins))

    lower_floor = max(1, int(min_bins))
    upper_cap = max(upper_cap, lower_floor)

    return int(np.clip(k_scaled, lower_floor, upper_cap))


def secondary_asym_keep(
    asym_snrs,
    primary_snr,
    secondary_snr_threshold=None,
    network_snr_threshold=None,
    max_followup_bins=0,
):
    """Positions of the secondary triggers (SNRs asym_snrs) kept."""
    if len(asym_snrs) == 0:
        return np.empty(0, dtype=np.int64)

    rho_sec = np.abs(asym_snrs)
    has_sec = (secondary_snr_threshold is not None) and (secondary_snr_threshold > 0.0)
    has_net = (network_snr_threshold is not None) and (network_snr_threshold > 0.0)

    if has_sec or has_net:
        keep_mask = np.zeros(len(asym_snrs), dtype=bool)
        if has_sec:
            keep_mask |= rho_sec >= secondary_snr_threshold
        if has_net:
            keep_mask |= (
                np.sqrt(abs(primary_snr) ** 2 + rho_sec**2) >= network_snr_threshold
            )
        passing_indices = np.flatnonzero(keep_mask)
    else:
        passing_indices = np.arange(len(asym_snrs))

    if (
        max_followup_bins is not None
        and max_followup_bins > 0
        and len(passing_indices) > max_followup_bins
    ):
        order = np.argsort(-rho_sec[passing_indices], kind="stable")
        passing_indices = np.sort(passing_indices[order[:max_followup_bins]])
    return passing_indices


_SHARED_WORKER_CONTEXT = {}


class TopTemplateSearch:
    """One top template's search over every segment and detector: the hierarchy below it
    (the top as root, one inner filter per middle group, each group's fine bank as its
    leaf) run by matchedfilter's HierarchicalSearch, with this object's methods as the
    hooks.
    """

    def __init__(self, task_args, worker_context=None):
        _require_matchedfilter()
        from pycbc.filter.matchedfilter import get_cutoff_indices
        from pycbc.vetoes.chain import SignalConsistencyChain
        from pycbc.vetoes.chisq_fir import (
            ChisqFilters,
            RatioPowerChisq,
            chisq_middle_spectra,
        )

        top_idx, rank, n_top, fine_cnt = task_args
        ctx = _SHARED_WORKER_CONTEXT if worker_context is None else worker_context
        self.ctx = ctx
        self.bank = bank = ctx["bank"]
        self.instruments = instruments = ctx["instruments"]
        self.segments = segments = ctx["segments"]
        self.args = args = ctx["args"]
        self.threshold = args.asym_snr_threshold
        self.sample_rate = sample_rate = ctx["sample_rate"]
        self.inj_reject = ctx.get("inj_reject")
        self.sg_chisq = ctx.get("sg_chisq")
        self.top_idx = top_idx = int(top_idx)
        self.timers = collections.defaultdict(float)
        self.triggers = {ifo: collections.defaultdict(list) for ifo in instruments}
        self.asym = dict.fromkeys(instruments, 0)
        self.search = None

        u_data = ctx["upper_group_data"][top_idx]
        middle_ids = u_data["m_all"]
        if len(middle_ids) == 0:
            return
        counts_u = u_data["counts"]
        taps_u = u_data["taps"]
        sigmas_u = u_data["sigmas"]

        self.top_template = top_template = ctx["top_templates"][top_idx]
        self.sigmasq_ratio = {
            int(mid): float((sigmas_u[row, -1] / sigmas_u[row, 0]) ** 2)
            for row, mid in enumerate(middle_ids)
        }
        logging.info(
            "Top template %d/%d (id=%d) with %d middle groups (%d fine templates)",
            rank + 1,
            n_top,
            top_idx,
            len(middle_ids),
            fine_cnt,
        )

        upper_f_by_coarse = chisq_middle_spectra(
            taps_u, counts_u, middle_ids, int(args.chisq_block_length)
        )

        reach = (
            args.autochi_stride * args.autochi_number_points
            if args.autochi_stride > 0
            else 0
        )
        self.groups = {}
        inner_nodes = []
        t0 = time.perf_counter()
        for row, mid in enumerate(middle_ids):
            coarse_idx = int(mid)
            taps, counts, fine_indices, shifts = bank.get_firs(
                coarse_idx, return_shifts=True
            )
            children = []
            if len(fine_indices):
                chisq = ChisqFilters(
                    taps,
                    counts,
                    args.chisq_block_length,
                    reach=reach,
                    tap_sample_rate=bank.sample_rate,
                    data_sample_rate=sample_rate,
                )
                self.groups[coarse_idx] = (fine_indices, chisq, shifts)
                children.append(
                    mf_search.Leaf(
                        taps,
                        counts,
                        lags=shifts,
                        node_id=coarse_idx,
                    )
                )
            inner_nodes.append(
                mf_search.Inner(
                    taps_u[row],
                    count=int(counts_u[row]),
                    children=children,
                    node_id=coarse_idx,
                )
            )
        self.timers["filter_prep"] += time.perf_counter() - t0
        if not self.groups:
            return

        t0 = time.perf_counter()
        self.seg_ifos = seg_ifos = [
            ifo for ifo in instruments if len(segments[ifo]) > 0
        ]
        first_seg = segments[seg_ifos[0]][0]
        n_series = (len(first_seg) - 1) * 2
        kmin, kmax = get_cutoff_indices(
            top_template.f_lower,
            args.high_frequency_cutoff,
            ctx["delta_f"],
            n_series,
        )
        profile = reference_profile(
            top_template,
            first_seg.psd,
            float(top_template.f_lower or 0.0),
            float(args.high_frequency_cutoff or 0.0),
            ctx["delta_f"],
        )
        self.root = mf_search.Root(
            top_template.numpy(),
            (kmin, kmax),
            children=inner_nodes,
            profile=(profile, ctx["delta_f"]),
            node_id=top_idx,
        )

        self.decimation = decimation_factor(bank.sample_rate, sample_rate)
        min_pad = max(4096, int(args.chisq_block_length) // 2 + reach)
        engine_mode = {
            "matchedfilter": "flat",
            "matchedfilter-hierarchical": "hier",
        }.get(args.ratio_filter_engine, args.ratio_filter_engine)
        asym_bs = int(round(args.asym_bin_width * float(sample_rate)))
        asym_span = (args.asym_num_bins // 2 + 1) * asym_bs + int(
            0.1 * float(sample_rate)
        )
        self.search = mf_search.HierarchicalSearch(
            [self.root],
            threshold=self.threshold,
            leaf_options=dict(
                tap_sample_rate=int(bank.sample_rate),
                data_sample_rate=int(sample_rate),
                engine=engine_mode.lower(),
                false_dismissal=float(args.ratio_filter_false_dismissal),
                first_stage_snr=float(args.ratio_filter_first_stage_snr),
                coarse_band_hz=float(args.ratio_filter_band),
                device=args.ratio_filter_device or None,
                binsize=int(round(float(args.cluster_window) * int(sample_rate)))
                if args.cluster_window
                else None,
            ),
            context=asym_span,
            min_pad=min_pad,
        )
        self.timers["filter_prep"] += time.perf_counter() - t0
        self.top_pad = self.search.pad(self.root)
        guard_sec = float(args.cluster_window or 0.0) + float(args.chisq_window or 0.0)
        self.guard_samples = (
            int(np.ceil(guard_sec * float(sample_rate))) + 2 * int(self.top_pad)
        )

        self.power_chisq = RatioPowerChisq(
            args.chisq_bins,
            sample_rate,
            flow=ctx["flow"],
            f_final=min(
                float(args.high_frequency_cutoff or 800.0),
                float(
                    bank.table["f_final"].max()
                    if hasattr(bank, "table") and "f_final" in bank.table
                    else 800.0
                ),
            ),
            top_template=top_template,
            upper_f_by_coarse=upper_f_by_coarse,
        )
        self.chain = SignalConsistencyChain.from_cli(
            args,
            bank,
            sample_rate,
            power_chisq=self.power_chisq,
            sg_chisq=self.sg_chisq,
        )

        self.group_tau = {}
        if self.inj_reject is not None:
            tau0_fine = self.inj_reject["tau0_fine"]
            cw = self.inj_reject["window"]
            for cid, (fine_indices, _, _) in self.groups.items():
                t_g = tau0_fine[fine_indices]
                self.group_tau[cid] = (float(t_g.min()) - cw, float(t_g.max()) + cw)
        self._norms = {}
        self.seg_clock = {
            "t": time.perf_counter(),
            "leaf": 0.0,
            "chi": 0.0,
            "asym": 0.0,
            "top": 0.0,
            "mid": 0.0,
        }

    def norm(self, ifo, s):
        """The top template's sigmasq over the noise of ifo's segment s."""
        key = (ifo, s)
        if key not in self._norms:
            self._norms[key] = self.top_template.sigmasq(self.segments[ifo][s].psd)
        return self._norms[key]

    def root_scale(self, ifo):
        """The top's output scale for each of ifo's segments."""

        def scale(_root, s):
            st = self.segments[ifo][s]
            norm = (4.0 * st.delta_f) / np.sqrt(self.norm(ifo, s))
            return 0.5 * (norm * st.delta_t) / self.decimation

        return scale

    def windows(self, leaf, s_num, ifo, extra):
        """A group's search windows in a segment, widened by the injection guard band
        and extra samples: None for the whole analysed range, empty to skip the group.
        """
        cid = getattr(leaf, "node_id", leaf)
        st = self.segments[ifo][s_num]
        a0, a1 = st.analyze.start, st.analyze.stop
        margin = self.guard_samples + extra
        if self.inj_reject is not None:
            injs = getattr(st, "analyze_injections", None)
            if injs is not None:
                lo, hi = self.group_tau[cid]
                tau0_ifo = self.inj_reject["tau0_inj"][ifo]
                return [
                    (max(a0, a - margin), min(a1, b + margin))
                    for a, b, row in injs
                    if lo <= tau0_ifo[row] <= hi
                ]
        wins = getattr(st, "analyze_windows", None)
        if wins is None:
            return None
        return [(max(a0, a - margin), min(a1, b + margin)) for a, b in wins]

    def on_peaks(self, ctx):
        """One middle group's fine peaks in a segment: veto and record them in each
        detector, then (in segments with a follow-up) follow the ones passing the
        primary gates up in the other detectors.
        """
        self.power_chisq.begin_group()
        primaries = {
            ifo: self._record(
                ctx, ifo, res.template_indices, res.sample_indices, res.snr
            )
            for ifo, res in ctx.peaks.items()
        }
        if self.ctx["can_asym_by_seg"].get(ctx.chunk, False):
            self.follow_up(ctx, primaries)

    def _record(
        self, ctx, ifo, ids_local, t_idxs, snrs, is_asym=False, threshold=None
    ):
        """Filter triggers to segment analyze bounds and SNR threshold, evaluate the
        signal-consistency chain, and append surviving triggers to self.triggers[ifo].
        """
        cid = ctx.leaf.node_id
        fine_indices, filters_f, shifts = self.groups[cid]
        stilde = ctx.tags[ifo]

        if len(ids_local) > 0 and getattr(stilde, "analyze", None) is not None:
            phys_t = t_idxs - (shifts[ids_local] if shifts is not None else 0)
            in_bounds = (phys_t >= stilde.analyze.start) & (
                phys_t < stilde.analyze.stop
            )
            if not np.all(in_bounds):
                ids_local = ids_local[in_bounds]
                t_idxs = t_idxs[in_bounds]
                snrs = snrs[in_bounds]

        thr = self.threshold if threshold is None else threshold
        if len(ids_local) > 0 and not is_asym and thr and thr > 0:
            keep = np.abs(snrs) >= thr
            if not np.all(keep):
                ids_local = ids_local[keep]
                t_idxs = t_idxs[keep]
                snrs = snrs[keep]

        n_trig = len(ids_local)
        if n_trig == 0:
            return (
                ids_local[:0],
                t_idxs[:0],
                snrs[:0],
                np.empty(0, dtype=float32),
                np.empty(0, dtype=np.int32),
            )

        if self.power_chisq.do or self.args.autochi_stride > 0:
            cthresh = (
                (self.args.chisq_snr_threshold or 0.0)
                if self.power_chisq.do
                else 0.0
            )
            above = np.abs(snrs) >= cthresh
            if np.any(above):
                t0 = time.perf_counter()
                filters_f.prefetch(ids_local[above])
                self.chain.stage_times["power"] += time.perf_counter() - t0

        local_shifts = shifts[ids_local] if shifts is not None else 0
        phys_t = t_idxs - local_shifts
        gids = fine_indices[ids_local]
        ref_snr = ctx.series(ifo)
        try:
            kept, vals = self.chain.values(
                snrs,
                phys_t,
                gids,
                lambda j: self.power_chisq.template(
                    stilde, cid, int(ids_local[j]), int(gids[j]), filters_f, self.bank
                ),
                lambda j, reach=0: self.power_chisq.correlation(
                    stilde,
                    cid,
                    int(ids_local[j]),
                    int(t_idxs[j]),
                    ref_snr,
                    filters_f,
                    reach,
                ),
                stilde=stilde,
            )
        finally:
            self.power_chisq.ref_f_cache.clear()

        if not np.all(kept):
            ids_local = ids_local[kept]
            t_idxs = t_idxs[kept]
            snrs = snrs[kept]
            phys_t = phys_t[kept]
            gids = gids[kept]
            vals = {k: v[kept] for k, v in vals.items()}
            n_trig = len(ids_local)
            if n_trig == 0:
                return (
                    ids_local[:0],
                    t_idxs[:0],
                    snrs[:0],
                    np.empty(0, dtype=float32),
                    np.empty(0, dtype=np.int32),
                )

        time_offset = (
            (-stilde.analyze.start + stilde.cumulative_index)
            if getattr(stilde, "analyze", None) is not None
            else 0
        )
        chi = (
            vals["chisq"].astype(float32, copy=False)
            if "chisq" in vals
            else np.zeros(n_trig, dtype=float32)
        )
        chi_dof = (
            vals["chisq_dof"].astype(np.int32, copy=False)
            if "chisq_dof" in vals
            else np.full(n_trig, -100, dtype=np.int32)
        )
        sg_vals = (
            vals["sg_chisq"].astype(float32, copy=False)
            if "sg_chisq" in vals
            else np.ones(n_trig, dtype=float32)
        )
        h_norm = self.norm(ifo, ctx.chunk) * self.sigmasq_ratio[cid]
        batch = {
            "template_id": gids,
            "time_index": (phys_t + time_offset).astype(np.int64),
            "snr": snrs.astype(complex64),
            "sigmasq": np.full(n_trig, h_norm, dtype=np.float32),
            "chisq": chi,
            "chisq_dof": chi_dof,
            "sg_chisq": sg_vals,
            "is_asym": np.full(n_trig, is_asym, dtype=bool),
        }
        if "cont_chisq" in vals:
            batch["cont_chisq"] = vals["cont_chisq"].astype(float32, copy=False)
        for col, arr in batch.items():
            self.triggers[ifo][col].append(arr)
        return ids_local, t_idxs, snrs, chi, chi_dof

    def follow_up(self, ctx, primaries):
        """Follow up each detector's recorded primaries that pass the NewSNR,
        reduced-chisq, and per-second group-cap gates in every other detector.
        """
        from pycbc.events import ranking

        t0 = time.perf_counter()
        chi0 = sum(self.chain.stage_times.values())
        args = self.args
        num_bins = args.asym_num_bins
        newsnr_min = args.asym_primary_newsnr_threshold
        if newsnr_min is None:
            newsnr_min = args.newsnr_threshold
        if newsnr_min is None:
            newsnr_min = self.threshold
        max_rchisq = args.asym_primary_max_rchisq
        if max_rchisq is not None and max_rchisq <= 0:
            max_rchisq = None
        cap = args.asym_max_primary_per_group_sec

        for ifo_a in ctx.channels:
            ids, t_idxs, p_snrs, p_chi, p_dof = primaries[ifo_a]
            if len(ids) == 0:
                continue
            rchisq = np.where(p_dof > 0, p_chi / np.maximum(1, p_dof), 1.0)
            gate = np.ones(len(p_snrs), dtype=bool)
            if newsnr_min is not None and newsnr_min > 0:
                gate &= ranking.newsnr(np.abs(p_snrs), rchisq) >= newsnr_min
            if max_rchisq is not None:
                gate &= rchisq <= max_rchisq
            followed = np.flatnonzero(gate)
            if len(followed) == 0:
                continue

            stilde_a = ctx.tags[ifo_a]
            t0_a = float(stilde_a.start_time)
            dt_a = float(stilde_a.delta_t)
            if cap is not None and cap > 0 and len(followed) > cap:
                sec_bins = np.floor(t0_a + t_idxs[followed] * dt_a).astype(np.int64)
                throttled = []
                for s in np.unique(sec_bins):
                    s_cands = followed[sec_bins == s]
                    if len(s_cands) <= cap:
                        throttled.extend(s_cands)
                    else:
                        order = np.argsort(-np.abs(p_snrs[s_cands]), kind="stable")
                        throttled.extend(s_cands[order[:cap]])
                followed = np.array(sorted(throttled), dtype=np.int64)

            for ifo_b in ctx.channels:
                if ifo_b == ifo_a:
                    continue
                stilde_b = ctx.tags[ifo_b]
                t0_b = float(stilde_b.start_time)
                dt_b = float(stilde_b.delta_t)
                series_b = ctx.series(ifo_b)

                def refilter(ifo=ifo_b, **kw):
                    return ctx.refilter(ifo, **kw)

                for k in followed:
                    gps_trig = t0_a + t_idxs[k] * dt_a
                    t_center = int(round((gps_trig - t0_b) / dt_b))
                    b_ids, b_tidxs, b_snrs = sample_asym_bins(
                        refilter=refilter,
                        series_length=len(series_b),
                        local_id=ids[k],
                        t_center=t_center,
                        sample_rate=self.sample_rate,
                        bin_width=args.asym_bin_width,
                        num_bins=num_bins,
                        valid_slice=stilde_b.analyze,
                    )
                    self.asym[ifo_b] += 1

                    if args.asym_dynamic_followup_bins:
                        max_bins = compute_followup_bin_budget(
                            primary_snr=p_snrs[k],
                            num_bins=num_bins,
                            bin_width=args.asym_bin_width,
                            min_bins=args.asym_min_followup_bins,
                            scale_factor=args.asym_followup_bin_scale,
                            target_network_snr=args.asym_followup_target_net_snr,
                            secondary_snr_floor=args.asym_followup_sec_snr_floor,
                            max_bins=args.asym_max_followup_bins,
                        )
                    else:
                        max_bins = args.asym_max_followup_bins

                    keep = secondary_asym_keep(
                        b_snrs,
                        primary_snr=p_snrs[k],
                        secondary_snr_threshold=args.asym_secondary_snr_threshold,
                        network_snr_threshold=args.asym_network_snr_threshold,
                        max_followup_bins=max_bins,
                    )
                    self._record(
                        ctx,
                        ifo_b,
                        b_ids[keep],
                        b_tidxs[keep],
                        b_snrs[keep],
                        is_asym=True,
                        threshold=args.asym_secondary_snr_threshold,
                    )

        chi = sum(self.chain.stage_times.values()) - chi0
        self.timers["asym_followup"] += time.perf_counter() - t0 - chi

    def on_chunk(self, ctx):
        """A segment is done: log its time and throughput."""
        now = time.perf_counter()
        clock, timings = self.seg_clock, self.search.timings
        t_seg_dur = max(now - clock["t"], 1e-6)
        ifo0 = ctx.channels[0]
        active_cids = [
            cid
            for cid in self.groups
            if self.windows(cid, ctx.chunk, ifo0, 0) != []
        ]
        n_tmpl_seg = sum(len(self.groups[cid][0]) for cid in active_cids)
        stilde_ref = ctx.tags[ifo0]
        total_analyzed_dur = (
            len(ctx.channels) * analyzed_samples(stilde_ref) * stilde_ref.delta_t
        )
        tperf_seg = total_analyzed_dur / t_seg_dur * n_tmpl_seg
        chi_now = self.chain.stage_times["power"]
        top_now = timings["root"] + self.timers["upper_ref"]
        logging.info(
            "[Top %d] Seg %d/%d: %d batches (%d templates) in %.2f s "
            "(%.0f temp-s/s) [fine=%.2fs chi=%.2fs asym=%.2fs top=%.3fs mid=%.3fs]",
            self.top_idx,
            ctx.chunk + 1,
            self.ctx["num_segs"],
            len(active_cids),
            n_tmpl_seg,
            t_seg_dur,
            tperf_seg,
            timings["leaf"] - clock["leaf"],
            chi_now - clock["chi"],
            self.timers["asym_followup"] - clock["asym"],
            top_now - clock["top"],
            timings["inner"] - clock["mid"],
        )
        clock.update(
            t=now,
            leaf=timings["leaf"],
            chi=chi_now,
            asym=self.timers["asym_followup"],
            top=top_now,
            mid=timings["inner"],
        )

    def run(self):
        """Search every segment; returns (top id, triggers per detector, follow-up
        counts, timers) for the worker pool.
        """
        if self.search is None:
            return (
                self.top_idx,
                {ifo: {} for ifo in self.instruments},
                self.asym,
                dict(self.timers),
            )
        inputs = {
            ifo: mf_search.BlockInput(
                [st.numpy() for st in self.segments[ifo]],
                [(st.analyze.start, st.analyze.stop) for st in self.segments[ifo]],
                self.root_scale(ifo),
                tags=list(self.segments[ifo]),
            )
            for ifo in self.seg_ifos
        }
        self.search.run(
            inputs,
            windows=self.windows,
            on_peaks=self.on_peaks,
            on_chunk=self.on_chunk,
        )
        timings = self.search.timings
        self.timers["upper_ref"] += timings["root"]
        self.timers["upper_fir"] += timings["inner"]
        self.timers["fine_fir"] += timings["leaf"]
        self.timers["filter_prep"] += timings["setup"]
        self.timers["power_chisq"] += self.chain.stage_times["power"]
        self.timers["auto_chisq"] += self.chain.stage_times["auto"]
        self.timers["sg_chisq"] += self.chain.stage_times["sg"]
        self.power_chisq.clear()
        collapsed = {
            ifo: {
                col: np.concatenate(arrs) if arrs else np.empty(0)
                for col, arrs in self.triggers[ifo].items()
            }
            for ifo in self.instruments
        }
        return self.top_idx, collapsed, self.asym, dict(self.timers)


def process_top_template_task(task_args):
    """Filter a single top template across all data segments and detectors."""
    return TopTemplateSearch(task_args).run()


def _validate_bank_and_padding(bank, segments, instruments, bank_file):
    """Verify 2-level real-above/complex-leaf FIR bank and segment padding."""
    if bank.num_levels < 2:
        raise ValueError(
            "pycbc_inspiral_fir requires a hierarchical bank with at least "
            "2 FIR levels."
        )
    for lvl, grp in bank.level_groups.items():
        want_complex = lvl == 1
        for k in bank.level_keys[lvl]:
            if (
                np.issubdtype(grp[k]["taps"].dtype, np.complexfloating)
                != want_complex
            ):
                raise ValueError(
                    "%s: level_%d group %s has %s taps; pycbc_inspiral_fir "
                    "needs complex final (level 1) taps and real taps above "
                    "(rebuild the bank with the current pycbc_fir_bank)"
                    % (bank_file, lvl, k, grp[k]["taps"].dtype)
                )

    max_upper_taps = int(bank.level_groups[2].attrs.get("max_taps", 0)) or max(
        int(c)
        for top_idx in bank.level_indices[2]
        for c in bank.level_groups[2][str(int(top_idx))]["actual_tap_count"][:]
    )
    max_lower_taps = int(bank.level_groups[1].attrs.get("max_taps", 0)) or max(
        int(c)
        for middle_idx in bank.coarse_indices
        for c in bank.level_groups[1][str(int(middle_idx))]["actual_tap_count"][:]
    )
    edge_margin = (max_upper_taps + max_lower_taps) // 2
    for ifo in instruments:
        for stilde in segments[ifo]:
            if (
                stilde.analyze.start < edge_margin
                or (len(stilde) - 1) * 2 - stilde.analyze.stop < edge_margin
            ):
                raise ValueError(
                    "analysis window needs %d samples of padding for the "
                    "composed upper and lower FIRs" % edge_margin
                )


def run_hierarchical_search(
    bank,
    segments,
    instruments,
    args,
    flow,
    delta_f,
    sample_rate,
    sg_chisq=None,
    inj_reject=None,
):
    """Execute the 2-level hierarchical FIR search across all top templates
    (serially or with a fork worker pool) and return (triggers, has_asym_followup,
    timers).
    """
    _validate_bank_and_padding(bank, segments, instruments, args.bank_file)

    triggers = {ifo: collections.defaultdict(list) for ifo in instruments}
    has_asym_followup = dict.fromkeys(instruments, 0)
    timers = collections.defaultdict(float)
    num_segs = max(len(segments[ifo]) for ifo in instruments)
    if sum(len(segments[ifo]) for ifo in instruments) == 0:
        logging.info(
            "Zero valid segments to analyze across all instruments. "
            "Skipping filtering."
        )
        return triggers, has_asym_followup, timers

    coarse_keep = (
        set(bank.coarse_indices)
        if len(bank.coarse_indices) < len(bank.coarse_keys)
        else None
    )
    top_templates = {}
    upper_group_data = {}
    top_tasks = []
    for raw_idx in bank.level_indices[2]:
        top_idx = int(raw_idx)
        grp = bank.level_groups[2][str(top_idx)]
        counts_raw = grp["actual_tap_count"][:]
        order = np.argsort(counts_raw, kind="stable")
        m_all = grp["fine_bank_index"][:][order]
        keep = (
            np.flatnonzero(np.isin(m_all, list(coarse_keep)))
            if coarse_keep is not None
            else slice(None)
        )
        m_kept = m_all[keep]
        if len(m_kept) == 0:
            continue
        top_templates[top_idx] = bank.get_coarse_template(top_idx)
        upper_group_data[top_idx] = {
            "counts": counts_raw[order][keep],
            "m_all": m_kept,
            "taps": grp["taps"][:][order][keep],
            "sigmas": grp["sigmas"][:][order][keep],
        }
        fine_count = sum(bank.group_template_count(int(m)) for m in m_kept)
        top_tasks.append((top_idx, fine_count))

    top_tasks.sort(key=lambda x: x[1], reverse=True)
    n_top = len(top_tasks)
    sorted_top_tasks = [
        (top_idx, rank, n_top, fine_cnt)
        for rank, (top_idx, fine_cnt) in enumerate(top_tasks)
    ]

    _SHARED_WORKER_CONTEXT.clear()
    _SHARED_WORKER_CONTEXT.update(
        {
            "bank": bank,
            "top_templates": top_templates,
            "upper_group_data": upper_group_data,
            "instruments": instruments,
            "segments": segments,
            "num_segs": num_segs,
            "can_asym_by_seg": segments_with_asym_followup(segments, instruments),
            "args": args,
            "delta_f": delta_f,
            "sample_rate": sample_rate,
            "flow": flow,
            "inj_reject": inj_reject,
            "sg_chisq": sg_chisq,
        }
    )

    total_cores = max(1, int(args.cores or 1))
    n_workers = min(total_cores, n_top) if n_top > 0 else 1

    def _collect(results_iter):
        for _, worker_trigs, worker_asym, worker_timers in results_iter:
            for ifo in instruments:
                for col, arr in worker_trigs[ifo].items():
                    if len(arr) > 0:
                        triggers[ifo][col].append(arr)
                has_asym_followup[ifo] += worker_asym[ifo]
            for k, v in worker_timers.items():
                timers[k] += v

    try:
        if n_workers > 1:
            logging.info(
                "Parallelizing %d top templates across %d cores using fork "
                "multiprocessing...",
                n_top,
                n_workers,
            )
            bank.close_filehandler()
            ctx = mp.get_context("fork")
            with ctx.Pool(processes=n_workers) as pool:
                _collect(
                    pool.imap_unordered(process_top_template_task, sorted_top_tasks)
                )
        else:
            logging.info("Filtering %d top templates serially (1 core)...", n_top)
            _collect(map(process_top_template_task, sorted_top_tasks))
    finally:
        _SHARED_WORKER_CONTEXT.clear()

    logging.info(
        "Upper synthesis time: direct-reference %.3f s, FIR %.3f s",
        timers["upper_ref"],
        timers["upper_fir"],
    )
    return triggers, has_asym_followup, timers
