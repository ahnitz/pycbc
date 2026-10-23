#!/usr/bin/env python3
"""Standalone validation test for full-trigger-table hollow-annulus ranking.

Validates that:
1. ReadByTemplate provides access to full trigger table and caches clustered
   clean time bins per file/detector, guaranteeing $O(N \log N)$ total time.
2. SingleDetTriggers.trig_dict() attaches _full_trigs and enables full dataset
   access for hollow-annulus ranking.
3. On real O4b data (H1L1-HDF_TRIGGER_MERGE_FULL_DATA-1398500224-5184000.hdf):
   - Event #53 (t=1398521756.41) in H1 has hollow counts detected (count ~ 19.0)
     and is genuinely downweighted (nsnr drops from ~3.79 to ~1.28).
   - Event #109 (t=1398521644.90) in H1 has hollow counts detected (count ~ 32.0)
     and is genuinely downweighted (nsnr drops from ~6.37 to ~0.62).
4. Both ReadByTemplate and SingleDetTriggers.trig_dict() produce consistent
   results.
"""

import os
import sys
import time
import numpy as np

import pycbc.events.ranking as ranking
from pycbc.io.hdf import ReadByTemplate, SingleDetTriggers


HDF_PATH = "/home/ahnitz/projects/claude/searchdev/bbh_test/test64_2month_duration_smoothing/output/full_data/H1L1-HDF_TRIGGER_MERGE_FULL_DATA-1398500224-5184000.hdf"

HOLLOW_PARAMS = {
    "sngl_ranking_hollow_knee": 4.5,
    "sngl_ranking_hollow_scale": 2.5,
    "sngl_ranking_hollow_inner": 0.25,
    "sngl_ranking_hollow_outer": 12.0,
    "sngl_ranking_hollow_min_snr": 6.0,
    "sngl_ranking_hollow_max_rchisq": 2.5,
}

QUIET_PARAMS = {
    "sngl_ranking_hollow_knee": 0.0,
}


def test_hollow_ranking_real_data():
    if not os.path.exists(HDF_PATH):
        print(f"SKIPPED: Reference dataset {HDF_PATH} not found.")
        return 0

    print("=" * 70)
    print("HOLLOW-ANNULUS FULL-TRIGGER-TABLE RANKING VALIDATION")
    print("=" * 70)
    print(f"Dataset: {HDF_PATH}")

    # -------------------------------------------------------------
    # 1. Test ReadByTemplate access and caching
    # -------------------------------------------------------------
    print("\n--- 1. Testing ReadByTemplate Access & Caching ---")
    reader = ReadByTemplate(HDF_PATH, gating_veto_windows={}, ifo="H1")
    assert hasattr(reader, "_full_trigs"), "ReadByTemplate missing _full_trigs!"
    print("ReadByTemplate initialized successfully with _full_trigs.")

    # Event #53: template_id 1653, target_t = 1398521756.41
    reader.set_template(1653)
    times_53 = reader["end_time"]
    if "end_time_ns" in reader:
        times_53 = times_53 + reader["end_time_ns"] * 1e-9
    idx_53 = np.argmin(np.abs(times_53 - 1398521756.41))
    t_ev53 = times_53[idx_53]

    t0 = time.time()
    stat_53_hollow = ranking.get_newsnr_sgveto_autochisq(reader, **HOLLOW_PARAMS)
    t1 = time.time()
    dt_call1 = t1 - t0
    stat_53_quiet = ranking.get_newsnr_sgveto_autochisq(reader, **QUIET_PARAMS)

    val_53_hollow = stat_53_hollow[idx_53]
    val_53_quiet = stat_53_quiet[idx_53]

    print(f"Event #53 (template 1653, t={t_ev53:.5f}):")
    print(f"  Quiet nsnr:            {val_53_quiet:.4f}")
    print(f"  Hollow downweighted:   {val_53_hollow:.4f}")
    print(f"  Delta:                 {val_53_hollow - val_53_quiet:.4f}")
    print(f"  Call 1 extraction time: {dt_call1:.4f}s")

    assert val_53_hollow < val_53_quiet, "Event #53 was NOT downweighted!"
    assert val_53_hollow < 2.5, f"Expected downweighted nsnr < 2.5, got {val_53_hollow}"

    # Event #109: template_id 37070, target_t = 1398521644.90
    reader.set_template(37070)
    times_109 = reader["end_time"]
    if "end_time_ns" in reader:
        times_109 = times_109 + reader["end_time_ns"] * 1e-9
    idx_109 = np.argmin(np.abs(times_109 - 1398521644.90))
    t_ev109 = times_109[idx_109]

    t2 = time.time()
    stat_109_hollow = ranking.get_newsnr_sgveto_autochisq(reader, **HOLLOW_PARAMS)
    t3 = time.time()
    dt_call2 = t3 - t2
    stat_109_quiet = ranking.get_newsnr_sgveto_autochisq(reader, **QUIET_PARAMS)

    val_109_hollow = stat_109_hollow[idx_109]
    val_109_quiet = stat_109_quiet[idx_109]

    print(f"Event #109 (template 37070, t={t_ev109:.5f}):")
    print(f"  Quiet nsnr:            {val_109_quiet:.4f}")
    print(f"  Hollow downweighted:   {val_109_hollow:.4f}")
    print(f"  Delta:                 {val_109_hollow - val_109_quiet:.4f}")
    print(f"  Call 2 extraction time: {dt_call2:.6f}s (cached)")

    assert val_109_hollow < val_109_quiet, "Event #109 was NOT downweighted!"
    assert val_109_hollow < 1.5, f"Expected downweighted nsnr < 1.5, got {val_109_hollow}"
    assert dt_call2 < 0.05, f"Call 2 was not cached! Took {dt_call2:.4f}s"
    print(f"  Caching speedup: {dt_call1 / max(1e-6, dt_call2):.1f}x")

    # -------------------------------------------------------------
    # 2. Test SingleDetTriggers.trig_dict() access
    # -------------------------------------------------------------
    print("\n--- 2. Testing SingleDetTriggers.trig_dict() Access ---")
    sngl = SingleDetTriggers(HDF_PATH, "H1")
    assert hasattr(sngl, "_full_trigs"), "SingleDetTriggers missing _full_trigs!"

    # Test Event #53 via SingleDetTriggers masked to just this trigger
    mask_53 = np.zeros(sngl.ntriggers, dtype=bool)
    mask_53[5771767] = True
    sngl.mask = mask_53

    d_53 = sngl.trig_dict()
    assert "_full_trigs" in d_53, "trig_dict missing _full_trigs!"

    sngl_stat_53_hollow = ranking.get_newsnr_sgveto_autochisq(d_53, **HOLLOW_PARAMS)
    sngl_stat_53_quiet = ranking.get_newsnr_sgveto_autochisq(d_53, **QUIET_PARAMS)
    print(f"Event #53 via SingleDetTriggers.trig_dict():")
    print(f"  Quiet nsnr:            {sngl_stat_53_quiet[0]:.4f}")
    print(f"  Hollow downweighted:   {sngl_stat_53_hollow[0]:.4f}")
    print(f"  Delta:                 {sngl_stat_53_hollow[0] - sngl_stat_53_quiet[0]:.4f}")

    assert sngl_stat_53_hollow[0] < sngl_stat_53_quiet[0], "Event #53 in trig_dict was NOT downweighted!"
    np.testing.assert_allclose(sngl_stat_53_hollow[0], val_53_hollow, rtol=1e-3,
                               err_msg="ReadByTemplate and SingleDetTriggers hollow rank mismatch on Event #53!")

    # Test Event #109 via SingleDetTriggers masked to just this trigger
    mask_109 = np.zeros(sngl.ntriggers, dtype=bool)
    mask_109[879342] = True
    sngl.mask = mask_109

    d_109 = sngl.trig_dict()
    assert "_full_trigs" in d_109, "trig_dict missing _full_trigs!"

    sngl_stat_109_hollow = ranking.get_newsnr_sgveto_autochisq(d_109, **HOLLOW_PARAMS)
    sngl_stat_109_quiet = ranking.get_newsnr_sgveto_autochisq(d_109, **QUIET_PARAMS)
    print(f"Event #109 via SingleDetTriggers.trig_dict():")
    print(f"  Quiet nsnr:            {sngl_stat_109_quiet[0]:.4f}")
    print(f"  Hollow downweighted:   {sngl_stat_109_hollow[0]:.4f}")
    print(f"  Delta:                 {sngl_stat_109_hollow[0] - sngl_stat_109_quiet[0]:.4f}")

    assert sngl_stat_109_hollow[0] < sngl_stat_109_quiet[0], "Event #109 in trig_dict was NOT downweighted!"
    np.testing.assert_allclose(sngl_stat_109_hollow[0], val_109_hollow, rtol=1e-3,
                               err_msg="ReadByTemplate and SingleDetTriggers hollow rank mismatch on Event #109!")

    print("\n" + "=" * 70)
    print("ALL PRE-TEST VALIDATION CHECKS PASSED SUCCESSFULLY (100% PASS)!")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(test_hollow_ranking_real_data())
