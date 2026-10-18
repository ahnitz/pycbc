import collections
import logging
import time
from typing import NamedTuple, Optional
import numpy as np

from pycbc.types import complex64, float32

_LAZY_CHISQ_EXPORTS = (
    "ChisqFilters",
    "RatioAutoChisq",
    "RatioBlockInputs",
    "RatioPowerChisq",
    "chisq_middle_spectra",
)


def __getattr__(name):
    if name in _LAZY_CHISQ_EXPORTS:
        from pycbc.vetoes import chisq_fir

        return getattr(chisq_fir, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


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


class IfoTriggerData(NamedTuple):
    """Typed single-detector trigger arrays produced by hierarchical FIR search."""

    template_id: np.ndarray
    time_index: np.ndarray
    snr: np.ndarray
    sigmasq: np.ndarray
    chisq: np.ndarray
    chisq_dof: np.ndarray
    sg_chisq: np.ndarray
    is_asym: np.ndarray
    cont_chisq: Optional[np.ndarray] = None
    psd_var_val: Optional[np.ndarray] = None

    def keys(self):
        return [
            f for f in self._fields if getattr(self, f) is not None
        ]

    def __contains__(self, key):
        return key in self._fields and getattr(self, key) is not None

    def __getitem__(self, key):
        if isinstance(key, int):
            return tuple.__getitem__(self, key)
        if key in self._fields:
            val = getattr(self, key)
            if val is not None:
                return val
        raise KeyError(key)

    def to_dict(self):
        return {f: getattr(self, f) for f in self.keys()}

    @classmethod
    def from_dict(cls, data):
        return cls(**{f: data[f] for f in cls._fields if f in data})


def validate_dynamic_snr_options(variance_floor, max_boost, disable=False):
    """Validate dynamic SNR renormalization CLI / runtime options."""
    if disable:
        return
    vfloor = float(variance_floor)
    if not np.isfinite(vfloor) or vfloor <= 0.0:
        raise ValueError(
            f"--dynamic-snr-variance-floor must be finite and positive, got {variance_floor}"
        )
    if max_boost is not None:
        mboost = float(max_boost)
        if not np.isfinite(mboost) or mboost <= 0.0 or mboost > 1e3:
            raise ValueError(
                f"--dynamic-snr-max-boost must be in (0, 1000], got {max_boost}"
            )
    if vfloor < 0.1 and (max_boost is None or float(max_boost) > 10.0):
        raise ValueError(
            f"--dynamic-snr-variance-floor={vfloor} (< 0.1) requires "
            "--dynamic-snr-max-boost <= 10.0 to prevent extreme threshold collapse"
        )


def decimation_factor(tap_sample_rate, data_sample_rate):
    """The integer ratio of the bank's tap sample rate to the data's; raises if it is not one."""
    exact = float(tap_sample_rate) / float(data_sample_rate)
    factor = int(np.round(exact))
    if abs(exact - factor) > 1e-5 or factor < 1:
        raise ValueError(
            f"The bank sample rate ({tap_sample_rate} Hz) must be an exact integer multiple "
            f"of the data sample rate ({data_sample_rate} Hz); the ratio is {exact:.4f}."
        )
    return factor


def leaf_bank_options(
    engine,
    false_dismissal,
    coarse_band_hz,
    first_stage_snr,
    device,
    peak_window,
    tap_sample_rate,
    data_sample_rate,
):
    """Keyword arguments for a fine bank's matchedfilter TimeDomainFilterBank (the search sets
    the threshold). engine is pycbc's --ratio-filter-engine name; peak_window (s) is the
    granularity of the peak search, one peak a template per span this long (pycbc's cluster
    window), so how matchedfilter blocks the series is its own choice.
    """
    _require_matchedfilter()
    decimation_factor(tap_sample_rate, data_sample_rate)
    mode = {
        "matchedfilter": "flat",
        "matchedfilter-hierarchical": "hier",
    }.get(engine, engine)
    return dict(
        tap_sample_rate=int(tap_sample_rate),
        data_sample_rate=int(data_sample_rate),
        engine=mode.lower(),
        false_dismissal=float(false_dismissal),
        first_stage_snr=float(first_stage_snr),
        coarse_band_hz=float(coarse_band_hz),
        device=device or None,
        binsize=int(round(float(peak_window) * int(data_sample_rate)))
        if peak_window
        else None,
    )


def get_dynamic_snr_renorm_factor(
    series,
    dt=None,
    window_duration=8.0,
    hollow_duration=0.5,
    variance_floor=1.0,
    scale=None,
    max_boost_factor=None,
):
    """Compute rolling window hollow dynamic SNR renormalization envelope."""
    factor = DynamicSNRRenormFactor(
        series,
        dt=dt,
        window_duration=window_duration,
        hollow_duration=hollow_duration,
        variance_floor=variance_floor,
        scale=scale,
        max_boost_factor=max_boost_factor,
    )
    return factor[np.arange(len(factor))]


class DynamicSNRRenormFactor:
    """The envelope of get_dynamic_snr_renorm_factor, evaluated where it is read."""

    _BLOCK = 16

    def __init__(
        self,
        series,
        dt=None,
        window_duration=8.0,
        hollow_duration=0.5,
        variance_floor=1.0,
        scale=None,
        max_boost_factor=None,
    ):
        self._n = n = len(series)
        self._ones = n == 0
        if dt is None:
            dt = float(getattr(series, "delta_t", 1.0))
        else:
            dt = float(dt)
        if dt <= 0:
            self._ones = True
        if self._ones:
            self._dtype = np.float32
            return

        if window_duration <= 2.0 * hollow_duration:
            raise ValueError(
                f"window_duration ({window_duration}s) must be greater than 2 * hollow_duration ({2.0 * hollow_duration}s)"
            )

        arr = np.asarray(series)
        if np.iscomplexobj(arr):
            power = arr.real * arr.real
            power += arr.imag * arr.imag
            half = 0.5
        else:
            power = arr * arr
            half = 1.0
        if scale is not None and float(scale) > 0:
            self._factor = half / (float(scale) * float(scale))
        else:
            self._factor = half

        blocks = self._block_prefix(power)
        if not np.isfinite(blocks[-1]):
            np.nan_to_num(power, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
            blocks = self._block_prefix(power)
        self._power, self._blocks = power, blocks

        self._w_outer = int(round((0.5 * window_duration) / dt))
        self._w_inner = int(round(hollow_duration / dt))
        self._floor = max(float(variance_floor), 1e-12)
        self._max_boost = max_boost_factor
        self._dtype = np.float64 if arr.dtype == np.complex128 else np.float32

    _PARTIAL = np.arange(_BLOCK - 1)

    @classmethod
    def _block_prefix(cls, power):
        """blocks[j]: the float64 sum of power over samples [0, j * _BLOCK)."""
        B, n = cls._BLOCK, len(power)
        whole = n // B
        blocks = np.zeros(-(-n // B) + 1, dtype=np.float64)
        blocks[1 : whole + 1] = (
            power[: whole * B].reshape(whole, B).sum(axis=1, dtype=np.float64)
        )
        if len(blocks) > whole + 1:
            blocks[-1] = power[whole * B :].sum(dtype=np.float64)
        return np.cumsum(blocks, out=blocks)

    def _prefix(self, m):
        """Sum of the power over samples [0, m), for a 1D integer array m in [0, n]."""
        B = self._BLOCK
        base = m - m % B
        v = self._power[np.minimum(base[:, None] + self._PARTIAL, self._n - 1)]
        v[self._PARTIAL >= (m - base)[:, None]] = 0
        return self._blocks[m // B] + v.sum(axis=1, dtype=np.float64)

    def __len__(self):
        return self._n

    def __getitem__(self, idx):
        idx = np.asarray(idx)
        if self._ones:
            return np.ones(idx.shape, dtype=self._dtype)
        n = self._n
        l_out = np.clip(idx - self._w_outer, 0, n)
        r_out = np.clip(idx + self._w_outer + 1, 0, n)
        l_in = np.clip(idx - self._w_inner, 0, n)
        r_in = np.clip(idx + self._w_inner + 1, 0, n)

        edges = self._prefix(
            np.concatenate(
                [r_out.ravel(), l_out.ravel(), r_in.ravel(), l_in.ravel()]
            )
        )
        e_ro, e_lo, e_ri, e_li = edges.reshape(4, *idx.shape)
        sum_hollow = self._factor * ((e_ro - e_lo) - (e_ri - e_li))
        cnt_hollow = (r_out - l_out) - (r_in - l_in)

        var_est = sum_hollow / np.maximum(cnt_hollow, 1)
        var_eff = np.where(
            np.isfinite(var_est), np.maximum(var_est, self._floor), self._floor
        )
        renorm_factor = np.where(
            np.isfinite(var_eff), 1.0 / np.sqrt(var_eff), 1.0
        )
        if self._max_boost is not None:
            renorm_factor = np.clip(renorm_factor, 0.0, float(self._max_boost))
        return renorm_factor.astype(self._dtype)


def dynamic_snr_renormalize(
    series,
    dt=None,
    window_duration=8.0,
    hollow_duration=0.5,
    variance_floor=1.0,
    scale=None,
    max_boost_factor=None,
):
    """Dynamically renormalize SNR series using a rolling window with hollow exclusion."""
    n = len(series)
    if n == 0:
        return series

    factor = get_dynamic_snr_renorm_factor(
        series,
        dt=dt,
        window_duration=window_duration,
        hollow_duration=hollow_duration,
        variance_floor=variance_floor,
        scale=scale,
        max_boost_factor=max_boost_factor,
    )
    return series * factor


def process_and_record_triggers(
    ids_local,
    t_idxs,
    snrs,
    stilde,
    h_norm,
    coarse_idx,
    fine_indices,
    filters_f,
    ref_snr,
    s_num,
    ifo,
    triggers,
    args,
    power_chisq,
    timers,
    is_asym=False,
    shifts=None,
    renorm=None,
    primary_snr_threshold=6.0,
    autochisq=None,
    sg_chisq=None,
    bank=None,
):
    """Keep the triggers inside the segment's analysis range (and, for primaries, at or above
    primary_snr_threshold after renormalisation), compute their chisq and append them to
    triggers[ifo]. Returns (positions kept in the inputs, ids, samples, renormalised SNRs,
    chisq, chisq dof) of the recorded triggers.
    """
    n_in = len(ids_local)
    surviving_indices = np.arange(n_in)

    if n_in > 0 and getattr(stilde, "analyze", None) is not None:
        local_shifts_check = shifts[ids_local] if shifts is not None else 0
        phys_t = t_idxs - local_shifts_check
        in_bounds = (phys_t >= stilde.analyze.start) & (
            phys_t < stilde.analyze.stop
        )
        if not np.all(in_bounds):
            ids_local = ids_local[in_bounds]
            t_idxs = t_idxs[in_bounds]
            snrs = snrs[in_bounds]
            surviving_indices = surviving_indices[in_bounds]
            if renorm is not None:
                renorm = renorm[in_bounds]

    snrs_renorm = (
        snrs * renorm if renorm is not None and len(ids_local) > 0 else snrs
    )

    if len(ids_local) > 0 and not is_asym:
        thresh = primary_snr_threshold
        if thresh is not None and thresh > 0:
            keep = np.abs(snrs_renorm) >= thresh
            if not np.all(keep):
                ids_local = ids_local[keep]
                t_idxs = t_idxs[keep]
                snrs = snrs[keep]
                snrs_renorm = snrs_renorm[keep]
                surviving_indices = surviving_indices[keep]

    n_trig = len(ids_local)
    if n_trig > 0:
        auto_on = autochisq is not None and getattr(autochisq, "do", False)
        sg_on = sg_chisq is not None and getattr(sg_chisq, "do", False)
        if (
            power_chisq is not None
            and power_chisq.do
            and np.any(abs(snrs_renorm) >= (power_chisq.snr_threshold or 0.0))
        ):
            t_chi0 = time.perf_counter()
            auto0 = autochisq.elapsed if auto_on else 0.0
            custom_steps = getattr(args, "signal_consistency_chain", None)
            max_rchi = getattr(
                args,
                "autochi_max_rchisq",
                getattr(args, "autochisq_max_rchisq", 5.0),
            )
            chain = power_chisq.chain(
                autochisq=autochisq if auto_on else None,
                max_rchi=max_rchi,
                sg_chisq=sg_chisq if sg_on else None,
                steps=custom_steps,
            )
            sg0 = chain.stage_times.get("sg", 0.0)
            chi, chi_dof, auto, _ = power_chisq.compute(
                ids_local,
                t_idxs,
                snrs_renorm,
                fine_indices,
                stilde,
                ref_snr,
                filters_f,
                coarse_idx,
                args,
                autochisq=autochisq if auto_on else None,
                sg_chisq=sg_chisq if sg_on else None,
                bank=bank,
                shifts=shifts,
            )
            auto_dur = autochisq.elapsed - auto0 if auto_on else 0.0
            sg_dur = chain.stage_times.get("sg", 0.0) - sg0 if sg_on else 0.0
            chi_dur = time.perf_counter() - t_chi0
            timers["power_chisq"] += max(0.0, chi_dur - auto_dur - sg_dur)
            timers["auto_chisq"] += auto_dur
            timers["sg_chisq"] += sg_dur
            sg_vals = (
                getattr(
                    power_chisq,
                    "last_sg_vals",
                    np.ones(n_trig, dtype=float32),
                )
                if sg_on
                else np.zeros(n_trig, dtype=float32)
            )
            chain_kept = getattr(power_chisq, "last_kept", None)
            if chain_kept is not None and not np.all(chain_kept):
                ids_local = ids_local[chain_kept]
                t_idxs = t_idxs[chain_kept]
                snrs = snrs[chain_kept]
                snrs_renorm = snrs_renorm[chain_kept]
                surviving_indices = surviving_indices[chain_kept]
                chi = chi[chain_kept]
                chi_dof = chi_dof[chain_kept]
                if auto is not None:
                    auto = auto[chain_kept]
                sg_vals = sg_vals[chain_kept]
                n_trig = len(ids_local)
                if n_trig == 0:
                    return (
                        np.empty(0, dtype=np.int64),
                        ids_local[:0],
                        t_idxs[:0],
                        snrs_renorm[:0],
                        np.empty(0, dtype=float32),
                        np.empty(0, dtype=np.int32),
                    )
        else:
            num_b = (
                int(getattr(power_chisq, "num_bins", 64) or 64)
                if power_chisq is not None
                else 64
            )
            chi = np.full(n_trig, 2 * num_b - 2, dtype=float32)
            chi_dof = np.full(n_trig, 2 * num_b - 2, dtype=np.int32)
            auto = (
                np.full(n_trig, autochisq.dof, dtype=float32)
                if auto_on
                else None
            )
            sg_vals = (
                np.ones(n_trig, dtype=float32)
                if sg_on
                else np.zeros(n_trig, dtype=float32)
            )

        local_shifts = shifts[ids_local] if shifts is not None else 0
        time_offset = (
            (-stilde.analyze.start + stilde.cumulative_index)
            if getattr(stilde, "analyze", None) is not None
            else 0
        )
        batch = {
            "template_id": fine_indices[ids_local],
            "time_index": (t_idxs + time_offset - local_shifts).astype(
                np.int64
            ),
            "snr": snrs_renorm.astype(complex64),
            "sigmasq": np.full(n_trig, h_norm, dtype=np.float32),
            "chisq": chi.astype(float32),
            "chisq_dof": chi_dof.astype(np.int32),
            "sg_chisq": sg_vals.astype(float32),
            "is_asym": np.full(n_trig, is_asym, dtype=bool),
        }
        if auto_on:
            batch["cont_chisq"] = auto
        ifo_trigs = triggers[ifo]
        for col, arr in batch.items():
            ifo_trigs[col].append(arr)
        return surviving_indices, ids_local, t_idxs, snrs_renorm, chi, chi_dof
    return (
        np.empty(0, dtype=np.int64),
        ids_local[:0],
        t_idxs[:0],
        snrs_renorm[:0],
        np.empty(0, dtype=float32),
        np.empty(0, dtype=np.int32),
    )


def segment_windows(stilde, margin=0):
    """Injection windows of a merged segment, widened by margin samples and
    clipped to the segment's analysis range, or None for a whole segment.
    """
    wins = getattr(stilde, "analyze_windows", None)
    if wins is None:
        return None
    a0, a1 = stilde.analyze.start, stilde.analyze.stop
    return [(max(a0, a - margin), min(a1, b + margin)) for a, b in wins]


def injection_group_windows(stilde, inj_tau, tau_lo, tau_hi, margin=0):
    """Windows of the injections in inj_tau, (start, stop, tau0) per injection,
    whose chirp time lies in [tau_lo, tau_hi], each widened by margin samples
    and clipped to the segment's analysis range.
    """
    a0, a1 = stilde.analyze.start, stilde.analyze.stop
    return [
        (max(a0, a - margin), min(a1, b + margin))
        for a, b, tau in inj_tau
        if tau_lo <= tau <= tau_hi
    ]


def injection_guard_samples(args, sample_rate, block_len):
    """Guard band searched around each injection window, in samples."""
    cw = float(args.cluster_window or 0.0)
    xw = float(args.chisq_window or 0.0)
    return int(np.ceil((cw + xw) * float(sample_rate))) + 2 * int(block_len)


def asym_half_span(args, sample_rate, num_bins):
    """How far an asym follow-up reaches from its primary, in samples, plus
    light travel between detectors.
    """
    bs = int(round(args.asym_bin_width * float(sample_rate)))
    return (num_bins // 2 + 1) * bs + int(0.1 * float(sample_rate))


def analyzed_samples(stilde):
    """Samples searched for primary triggers: the windows if any, else analyze."""
    wins = getattr(stilde, "analyze_windows", None)
    if wins is not None:
        return sum(b - a for a, b in wins)
    return stilde.analyze.stop - stilde.analyze.start


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


def select_primary_asym_candidates(
    a_ids,
    a_tidxs,
    a_snrs,
    sample_rate,
    asym_snr_threshold,
    max_primary_per_group_sec=1,
    stilde=None,
    gps_start_time=None,
):
    """Select candidate primary triggers to instigate secondary followup windowing."""
    if len(a_ids) == 0:
        return np.empty(0, dtype=np.int64)

    cand_mask = abs(a_snrs) >= asym_snr_threshold
    if not np.any(cand_mask):
        return np.empty(0, dtype=np.int64)

    cand_indices = np.flatnonzero(cand_mask)
    if len(cand_indices) <= 1:
        return cand_indices

    if (
        max_primary_per_group_sec is not None
        and max_primary_per_group_sec > 0
        and len(cand_indices) > max_primary_per_group_sec
    ):
        sample_rate_val = float(sample_rate or 2048.0)
        cum_idx = (
            getattr(stilde, "cumulative_index", 0) if stilde is not None else 0
        )
        an_start = (
            stilde.analyze.start
            if (stilde is not None and getattr(stilde, "analyze", None) is not None)
            else 0
        )
        gps_start_val = (
            float(gps_start_time)
            if gps_start_time is not None
            else float(getattr(stilde, "start_time", 0.0))
        )
        gps_times = (
            gps_start_val
            + (a_tidxs[cand_indices] - an_start + cum_idx) / sample_rate_val
        )
        sec_bins = np.floor(gps_times).astype(np.int64)

        unique_secs = np.unique(sec_bins)
        throttled = []
        for s in unique_secs:
            s_cands = cand_indices[sec_bins == s]
            if len(s_cands) <= max_primary_per_group_sec:
                throttled.extend(s_cands)
            else:
                s_order = np.argsort(-abs(a_snrs[s_cands]), kind="stable")
                throttled.extend(s_cands[s_order[:max_primary_per_group_sec]])
        cand_indices = np.array(sorted(throttled), dtype=np.int64)

    return cand_indices


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
    """Compute symmetric followup bin budget K(rho_prim) based on recovery difficulty."""
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
    if upper_cap < lower_floor:
        upper_cap = lower_floor

    return int(np.clip(k_scaled, lower_floor, upper_cap))


def secondary_asym_keep(
    asym_snrs,
    primary_snr,
    secondary_snr_threshold=None,
    network_snr_threshold=None,
    max_followup_bins=0,
):
    """Positions of the secondary triggers (renormalised SNRs asym_snrs) kept."""
    if len(asym_snrs) == 0:
        return np.empty(0, dtype=np.int64)

    rho_sec = np.abs(asym_snrs)
    has_sec = (secondary_snr_threshold is not None) and (
        secondary_snr_threshold > 0.0
    )
    has_net = (network_snr_threshold is not None) and (
        network_snr_threshold > 0.0
    )

    if has_sec or has_net:
        keep_mask = np.zeros(len(asym_snrs), dtype=bool)
        if has_sec:
            keep_mask |= rho_sec >= secondary_snr_threshold
        if has_net:
            keep_mask |= (
                np.sqrt(abs(primary_snr) ** 2 + rho_sec**2)
                >= network_snr_threshold
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
    (the top as root, one inner filter per middle group, each group's fine bank as its leaf)
    run by matchedfilter's HierarchicalSearch, with this object's methods as the hooks.
    """

    def __init__(self, task_args, worker_context=None):
        _require_matchedfilter()
        from pycbc.filter.matchedfilter import get_cutoff_indices
        from pycbc.vetoes.chisq_fir import (
            ChisqFilters,
            RatioAutoChisq,
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
        self.threshold = ctx["primary_snr_threshold"]
        self.sample_rate = sample_rate = ctx["sample_rate"]
        self.inj_reject = ctx.get("inj_reject")
        self.sg_chisq = ctx.get("sg_chisq")
        self.top_idx = top_idx = int(top_idx)
        self.timers = collections.defaultdict(float)
        self.triggers = {
            ifo: collections.defaultdict(list) for ifo in instruments
        }
        self.asym = {ifo: 0 for ifo in instruments}
        self.search = None

        u_data = ctx["upper_group_data"][top_idx]
        counts_raw = u_data["counts"]
        order_raw = np.argsort(counts_raw, kind="stable")
        m_all = u_data["m_all"][order_raw]
        keep_idx = slice(None)
        if ctx["coarse_keep"] is not None:
            row_mask = np.isin(m_all, list(ctx["coarse_keep"]))
            if not row_mask.any():
                return
            keep_idx = np.flatnonzero(row_mask)
        middle_ids = m_all[keep_idx]
        counts_u = counts_raw[order_raw][keep_idx]
        taps_u = u_data["taps"][order_raw][keep_idx]
        sigmas_u = u_data["sigmas"][order_raw][keep_idx]

        self.top_template = top_template = ctx["top_templates"][top_idx]
        self.sigmasq_ratio = {
            int(mid): float((sigmas_u[row, 2] / sigmas_u[row, 0]) ** 2)
            for row, mid in enumerate(middle_ids)
        }
        logging.info(
            f"Top template {rank+1}/{n_top} (id={top_idx}) with {len(middle_ids)} middle groups ({fine_cnt} fine templates)"
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

        self.renorm_on = not args.disable_dynamic_snr_renorm
        self.vfloor = args.dynamic_snr_variance_floor
        self.max_boost = args.dynamic_snr_max_boost
        validate_dynamic_snr_options(
            self.vfloor, self.max_boost, disable=not self.renorm_on
        )
        rf_max = 1.0 / np.sqrt(max(float(self.vfloor), 1e-12))
        if self.max_boost is not None:
            rf_max = min(rf_max, float(self.max_boost))
        self.decimation = decimation_factor(bank.sample_rate, sample_rate)
        min_pad = max(4096, int(args.chisq_block_length) // 2 + reach)
        self.search = mf_search.HierarchicalSearch(
            [self.root],
            threshold=self.threshold,
            envelope_floor=(1.0 / rf_max) if self.renorm_on else 1.0,
            leaf_options=leaf_bank_options(
                args.ratio_filter_engine,
                args.ratio_filter_false_dismissal,
                args.ratio_filter_band,
                args.ratio_filter_first_stage_snr,
                args.ratio_filter_device,
                args.cluster_window,
                bank.sample_rate,
                sample_rate,
            ),
            context=asym_half_span(
                args, sample_rate, ctx["effective_asym_num_bins"]
            ),
            min_pad=min_pad,
        )
        self.timers["filter_prep"] += time.perf_counter() - t0
        self.top_pad = self.search.pad(self.root)

        self.power_chisq = RatioPowerChisq(
            args.chisq_bins,
            args.chisq_snr_threshold,
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
            chisq_window=args.chisq_window,
            chisq_iterations=args.chisq_iterations,
            top_template=top_template,
            upper_f_by_coarse=upper_f_by_coarse,
        )
        self.autochisq = RatioAutoChisq(self.power_chisq, args)

        self.group_tau = {}
        self.rej_pairs = [0, 0]
        if self.inj_reject is not None:
            tau0_fine = self.inj_reject["tau0_fine"]
            cw = self.inj_reject["window"]
            for cid, (fine_indices, _, _) in self.groups.items():
                t_g = tau0_fine[fine_indices]
                self.group_tau[cid] = (
                    float(t_g.min()) - cw,
                    float(t_g.max()) + cw,
                )
        self.kept_by_seg = {}
        self._norms = {}
        self.renorm_factors = {}
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
            self._norms[key] = self.top_template.sigmasq(
                self.segments[ifo][s].psd
            )
        return self._norms[key]

    def root_scale(self, ifo):
        """The top's output scale for each of ifo's segments."""

        def scale(_root, s):
            st = self.segments[ifo][s]
            norm = (4.0 * st.delta_f) / np.sqrt(self.norm(ifo, s))
            return 0.5 * (norm * st.delta_t) / self.decimation

        return scale

    def kept_groups(self, s_num):
        """The groups searched in a segment (None: all) and its injections per detector."""
        if s_num not in self.kept_by_seg:
            kept, seg_inj = None, {}
            if self.inj_reject is not None:
                for ifo in self.seg_ifos:
                    if s_num < len(self.segments[ifo]):
                        injs = getattr(
                            self.segments[ifo][s_num],
                            "analyze_injections",
                            None,
                        )
                        seg_inj[ifo] = (
                            None
                            if injs is None
                            else [
                                (a, b, self.inj_reject["tau0_inj"][ifo][row])
                                for a, b, row in injs
                            ]
                        )
                if all(v is not None for v in seg_inj.values()):
                    seg_tau = [t for v in seg_inj.values() for _, _, t in v]
                    kept = {
                        cid
                        for cid, (lo, hi) in self.group_tau.items()
                        if any(lo <= t <= hi for t in seg_tau)
                    }
                self.rej_pairs[0] += (
                    len(self.groups) if kept is None else len(kept)
                )
                self.rej_pairs[1] += len(self.groups)
            self.kept_by_seg[s_num] = (kept, seg_inj)
        return self.kept_by_seg[s_num]

    def windows(self, leaf, s_num, ifo, extra):
        """A group's search windows in a segment, widened by the injection guard band and
        extra samples: None for the whole analysed range, empty to skip the group.
        """
        kept, seg_inj = self.kept_groups(s_num)
        cid = leaf.node_id
        if kept is not None and cid not in kept:
            return []
        st_g = self.segments[ifo][s_num]
        margin = (
            injection_guard_samples(
                self.args, self.ctx["ifo_sample_rates"][ifo], self.top_pad
            )
            + extra
        )
        if kept is None:
            return segment_windows(st_g, margin)
        lo, hi = self.group_tau[cid]
        return injection_group_windows(st_g, seg_inj[ifo], lo, hi, margin)

    def on_series(self, ctx):
        """The top's output per detector: each renormalisation envelope, built from it, and
        given to the search as the threshold envelope 1 / rf.
        """
        t0 = time.perf_counter()
        envelopes = {}
        for ifo in ctx.channels:
            if not self.renorm_on:
                self.renorm_factors[ifo] = None
                continue
            st = ctx.tags[ifo]
            rf = DynamicSNRRenormFactor(
                ctx.series[ifo],
                dt=float(st.delta_t),
                scale=float(st.delta_t / self.decimation),
                window_duration=self.args.dynamic_snr_window,
                hollow_duration=self.args.dynamic_snr_hollow,
                variance_floor=self.vfloor,
                max_boost_factor=self.max_boost,
            )
            self.renorm_factors[ifo] = rf
            envelopes[ifo] = lambda t, rf=rf: 1.0 / rf[t]
        self.timers["upper_ref"] += time.perf_counter() - t0
        return envelopes

    def on_peaks(self, ctx):
        """One middle group's fine peaks in a segment: renormalise, select, veto and record
        them in each detector, then (in segments with a follow-up) follow the ones passing the
        primary gates up in the other detectors.
        """
        args = self.args
        primaries = {}
        self.power_chisq.begin_group()
        for ifo in ctx.channels:
            res = ctx.peaks[ifo]
            e = ctx.envelope(ifo)
            renorm = (1.0 / e).astype(np.float32) if e is not None else None
            cand = select_primary_asym_candidates(
                a_ids=res.template_indices,
                a_tidxs=res.sample_indices,
                a_snrs=res.snr * renorm if renorm is not None else res.snr,
                sample_rate=self.ctx["ifo_sample_rates"][ifo],
                asym_snr_threshold=self.threshold,
                max_primary_per_group_sec=args.asym_max_primary_per_group_sec,
                stilde=ctx.tags[ifo],
                gps_start_time=self.ctx["ifo_gps_cli_starts"][ifo],
            )
            ids, t_idxs = res.template_indices[cand], res.sample_indices[cand]
            kept, _, _, p_snrs, p_chi, p_dof = self._record(
                ctx,
                ifo,
                ids,
                t_idxs,
                res.snr[cand],
                renorm[cand] if renorm is not None else None,
            )
            primaries[ifo] = (ids[kept], t_idxs[kept], p_snrs, p_chi, p_dof)
        if self.ctx["can_asym_by_seg"].get(ctx.chunk, False):
            self.follow_up(ctx, primaries)

    def _record(
        self, ctx, ifo, ids, t_idxs, snrs, renorm, is_asym=False, threshold=None
    ):
        """process_and_record_triggers for this group's triggers in ifo's segment."""
        cid = ctx.leaf.node_id
        fine_indices, chisq, shifts = self.groups[cid]
        return process_and_record_triggers(
            ids,
            t_idxs,
            snrs,
            ctx.tags[ifo],
            self.norm(ifo, ctx.chunk) * self.sigmasq_ratio[cid],
            cid,
            fine_indices,
            chisq,
            ctx.series(ifo),
            ctx.chunk,
            ifo,
            self.triggers,
            self.args,
            self.power_chisq,
            self.timers,
            is_asym=is_asym,
            shifts=shifts,
            renorm=renorm,
            primary_snr_threshold=self.threshold
            if threshold is None
            else threshold,
            autochisq=self.autochisq,
            sg_chisq=self.sg_chisq,
            bank=self.bank,
        )

    def follow_up(self, ctx, primaries):
        """Follow up each detector's recorded primaries that pass the NewSNR and
        reduced-chisq gates in every other detector.
        """
        from pycbc.events import ranking

        t0 = time.perf_counter()
        chi0 = (
            self.timers["power_chisq"]
            + self.timers["auto_chisq"]
            + self.timers["sg_chisq"]
        )
        args, s_num = self.args, ctx.chunk
        shifts = self.groups[ctx.leaf.node_id][2]
        rates, gps0, seg0 = (
            self.ctx["ifo_sample_rates"],
            self.ctx["ifo_gps_starts"],
            self.ctx["seg_starts"],
        )
        num_bins = self.ctx["effective_asym_num_bins"]
        newsnr_min = args.asym_primary_newsnr_threshold
        if newsnr_min is None:
            newsnr_min = args.newsnr_threshold
        if newsnr_min is None:
            newsnr_min = self.threshold
        max_rchisq = args.asym_primary_max_rchisq
        if max_rchisq is not None and max_rchisq <= 0:
            max_rchisq = None

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

            for ifo_b in ctx.channels:
                if ifo_b == ifo_a:
                    continue
                stilde_b = ctx.tags[ifo_b]
                series_b = ctx.series(ifo_b)
                rf_b = self.renorm_factors.get(ifo_b)
                refilter = lambda ifo=ifo_b, **kw: ctx.refilter(ifo, **kw)
                for k in followed:
                    gps_trig = (
                        gps0[ifo_a]
                        + (seg0[(ifo_a, s_num)] + t_idxs[k]) / rates[ifo_a]
                    )
                    t_center = (
                        int(round((gps_trig - gps0[ifo_b]) * rates[ifo_b]))
                        - seg0[(ifo_b, s_num)]
                    )
                    b_ids, b_tidxs, b_snrs = sample_asym_bins(
                        refilter=refilter,
                        series_length=len(series_b),
                        local_id=ids[k],
                        t_center=t_center,
                        sample_rate=rates[ifo_b],
                        bin_width=args.asym_bin_width,
                        num_bins=num_bins,
                        valid_slice=stilde_b.analyze,
                    )
                    self.asym[ifo_b] += 1

                    r_b = None
                    b_snrs_eval = b_snrs
                    if rf_b is not None and len(b_ids) > 0:
                        r_b = rf_b[
                            np.clip(b_tidxs - shifts[b_ids], 0, len(rf_b) - 1)
                        ]
                        b_snrs_eval = b_snrs * r_b

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
                        b_snrs_eval,
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
                        r_b[keep] if r_b is not None else None,
                        is_asym=True,
                        threshold=args.asym_secondary_snr_threshold,
                    )

        chi = (
            self.timers["power_chisq"]
            + self.timers["auto_chisq"]
            + self.timers["sg_chisq"]
            - chi0
        )
        self.timers["asym_followup"] += time.perf_counter() - t0 - chi

    def on_chunk(self, ctx):
        """A segment is done: log its time and throughput."""
        now = time.perf_counter()
        clock, timings = self.seg_clock, self.search.timings
        t_seg_dur = max(now - clock["t"], 1e-6)
        kept, _ = self.kept_groups(ctx.chunk)
        n_groups = len(self.groups) if kept is None else len(kept)
        n_tmpl_seg = sum(
            len(g[0])
            for cid, g in self.groups.items()
            if kept is None or cid in kept
        )
        stilde_ref = ctx.tags[ctx.channels[0]]
        total_analyzed_dur = (
            len(ctx.channels) * analyzed_samples(stilde_ref) * stilde_ref.delta_t
        )
        tperf_seg = total_analyzed_dur / t_seg_dur * n_tmpl_seg
        s_fine = timings["leaf"] - clock["leaf"]
        s_chi = self.timers["power_chisq"] - clock["chi"]
        s_asym = self.timers["asym_followup"] - clock["asym"]
        top_now = timings["root"] + self.timers["upper_ref"]
        s_top = top_now - clock["top"]
        s_mid = timings["inner"] - clock["mid"]
        logging.info(
            f"[Top {self.top_idx}] Seg {ctx.chunk+1}/{self.ctx['num_segs']}: {n_groups} batches ({n_tmpl_seg} templates) in {t_seg_dur:.2f} s ({tperf_seg:.0f} temp-s/s) [fine={s_fine:.2f}s chi={s_chi:.2f}s asym={s_asym:.2f}s top={s_top:.3f}s mid={s_mid:.3f}s]"
        )
        clock.update(
            t=now,
            leaf=timings["leaf"],
            chi=self.timers["power_chisq"],
            asym=self.timers["asym_followup"],
            top=top_now,
            mid=timings["inner"],
        )
        self.renorm_factors.clear()

    def run(self):
        """Search every segment; returns (top id, triggers per detector, follow-up counts,
        timers) for the worker pool.
        """
        if self.search is None:
            return (
                self.top_idx,
                {ifo: {} for ifo in self.instruments},
                self.asym,
                dict(self.timers),
            )
        inputs = {
            ifo: mf_search_block_input(self.segments[ifo], self.root_scale(ifo))
            for ifo in self.seg_ifos
        }
        self.search.run(
            inputs,
            windows=self.windows,
            on_series=self.on_series,
            on_peaks=self.on_peaks,
            on_chunk=self.on_chunk,
        )
        timings = self.search.timings
        self.timers["upper_ref"] += timings["root"]
        self.timers["upper_fir"] += timings["inner"]
        self.timers["fine_fir"] += timings["leaf"]
        self.timers["filter_prep"] += timings["setup"]
        if self.inj_reject is not None:
            logging.info(
                "[Top %d] injection rejection: %d of %d (group, segment) pairs filtered",
                self.top_idx,
                self.rej_pairs[0],
                self.rej_pairs[1],
            )
        self.power_chisq.clear()
        self.autochisq.clear()
        collapsed = {}
        for ifo in self.instruments:
            collapsed[ifo] = {
                col: np.concatenate(arrs) if arrs else np.empty(0)
                for col, arrs in self.triggers[ifo].items()
            }
        return self.top_idx, collapsed, self.asym, dict(self.timers)


def mf_search_block_input(segs, root_scale):
    """A detector's segments as matchedfilter search input blocks."""
    _require_matchedfilter()
    return mf_search.BlockInput(
        [st.numpy() for st in segs],
        [(st.analyze.start, st.analyze.stop) for st in segs],
        root_scale,
        tags=list(segs),
    )


def process_top_template_task(task_args):
    """Filter a single top template across all data segments and detectors."""
    return TopTemplateSearch(task_args).run()
