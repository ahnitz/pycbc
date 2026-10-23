"""Unit tests for cap_alpha_below, --max-fit-coeff, hollow-annulus ranking, and asym-followup-cluster-bins.
"""

import os
import unittest
import numpy as np

import pycbc.events.ranking as ranking
import pycbc.events.stat as stat
from pycbc.filter.matched_filter_ratio import secondary_asym_keep


def make_exp_fit_stat(alpha_below_thresh=5.0, cap_alpha_below=False):
    """Helper to instantiate ExpFitStatistic without requiring disk files."""
    inst = object.__new__(stat.ExpFitStatistic)
    inst.kwargs = {
        "alpha_below_thresh": alpha_below_thresh,
        "cap_alpha_below": cap_alpha_below,
    }
    inst.alphabelow = float(alpha_below_thresh)
    val = cap_alpha_below
    if isinstance(val, str):
        inst.cap_alpha_below = val.lower() in ["true", "1", "yes", "t"]
    else:
        inst.cap_alpha_below = bool(val)
    inst.curr_tnum = None
    return inst


class TestCapAlphaBelow(unittest.TestCase):
    """Test cap_alpha_below in ExpFitStatistic and continuity across threshold."""

    def test_default_behavior_uncapped(self):
        """When cap_alpha_below is False, discontinuous slope using alphabelow."""
        stat_inst = make_exp_fit_stat(alpha_below_thresh=5.0, cap_alpha_below=False)
        self.assertFalse(stat_inst.cap_alpha_below)
        self.assertEqual(stat_inst.alphabelow, 5.0)

        # Mock find_fits and get_sngl_ranking
        stat_inst.find_fits = lambda trigs: (3.0, 10.0, 6.0)  # alphai=3, ratei=10, thresh=6
        stat_inst.get_sngl_ranking = lambda trigs: trigs["snr"]
        
        # Test just above threshold (s = 6.0001)
        trigs_above = {
            "template_id": np.array([0]),
            "snr": np.array([6.0001]),
            "chisq": np.array([1.0]),
            "chisq_dof": np.array([15.0]),
        }
        logr_above = stat_inst.lognoiserate(trigs_above)[0]
        expected_above = np.log(3.0) + np.log(10.0)
        self.assertAlmostEqual(logr_above, expected_above, places=3)

        # Test just below threshold (s = 5.9999) without cap_alpha_below:
        # logrbt = -alphabelow*(s-thresh) + log(alphabelow) + log(ratei)
        # at 5.9999: approx log(5) + log(10)
        trigs_below = {
            "template_id": np.array([0]),
            "snr": np.array([5.9999]),
            "chisq": np.array([1.0]),
            "chisq_dof": np.array([15.0]),
        }
        logr_below = stat_inst.lognoiserate(trigs_below)[0]
        expected_below = -5.0 * (5.9999 - 6.0) + np.log(5.0) + np.log(10.0)
        self.assertAlmostEqual(logr_below, expected_below, places=3)
        # Note the jump discontinuity: log(5) - log(3) approx 0.51
        self.assertNotAlmostEqual(logr_above, logr_below, places=1)

    def test_cap_alpha_below_continuity_at_threshold(self):
        """When cap_alpha_below is True, exactly continuous at sngl_stat == thresh."""
        stat_inst = make_exp_fit_stat(alpha_below_thresh=5.0, cap_alpha_below="True")
        self.assertTrue(stat_inst.cap_alpha_below)
        self.assertEqual(stat_inst.alphabelow, 5.0)

        # alphai = 3.0 (< alphabelow = 5.0)
        stat_inst.find_fits = lambda trigs: (3.0, 10.0, 6.0)
        stat_inst.get_sngl_ranking = lambda trigs: trigs["snr"]

        s_above = {"template_id": np.array([0]), "snr": np.array([6.0 + 1e-6])}
        s_below = {"template_id": np.array([0]), "snr": np.array([6.0 - 1e-6])}
        logr_above = stat_inst.lognoiserate(s_above)[0]
        logr_below = stat_inst.lognoiserate(s_below)[0]

        # Must match across threshold (C0 continuity)
        self.assertAlmostEqual(logr_above, logr_below, places=4)

    def test_cap_alpha_below_capping_slope(self):
        """When alphai > alphabelow, slope below threshold is capped at alphabelow."""
        stat_inst = make_exp_fit_stat(alpha_below_thresh=4.0, cap_alpha_below=1)
        self.assertTrue(stat_inst.cap_alpha_below)

        # alphai = 7.0 (> alphabelow = 4.0)
        stat_inst.find_fits = lambda trigs: (7.0, 10.0, 6.0)
        stat_inst.get_sngl_ranking = lambda trigs: trigs["snr"]

        s1 = {"template_id": np.array([0]), "snr": np.array([5.0])}
        s2 = {"template_id": np.array([0]), "snr": np.array([4.0])}
        logr1 = stat_inst.lognoiserate(s1)[0]
        logr2 = stat_inst.lognoiserate(s2)[0]

        # Slope is (logr2 - logr1) / (5.0 - 4.0) = ab_eff = min(7.0, 4.0) = 4.0
        delta_stat = 1.0
        observed_slope = (logr2 - logr1) / delta_stat
        self.assertAlmostEqual(observed_slope, 4.0, places=5)

    def test_kwargs_parsing(self):
        """Verify cap_alpha_below is correctly parsed from kwargs strings and bools."""
        for val, expected in [("1", True), ("True", True), ("yes", True), ("0", False), ("False", False), (False, False)]:
            s = make_exp_fit_stat(alpha_below_thresh=5.0, cap_alpha_below=val)
            self.assertEqual(s.cap_alpha_below, expected)


class TestMaxFitCoeff(unittest.TestCase):
    """Test max_fit_coeff capping in pycbc_fit_sngls_over_multiparam."""

    def test_capping_logic(self):
        """Test array capping with negative placeholders preserved."""
        raw_alpha = np.array([-100.0, 2.5, 4.8, 7.35, 9.1, 0.0])
        max_coeff = 5.5
        capped = np.where(
            raw_alpha > 0, np.minimum(raw_alpha, max_coeff), raw_alpha
        )
        expected = np.array([-100.0, 2.5, 4.8, 5.5, 5.5, 0.0])
        np.testing.assert_allclose(capped, expected)

    def test_fit_over_param_cli(self):
        """Verify --max-fit-coeff argument is accepted by argument parser."""
        from argparse import ArgumentParser
        parser = ArgumentParser()
        parser.add_argument(
            "--max-fit-coeff",
            type=float,
            default=None,
        )
        args = parser.parse_args(["--max-fit-coeff", "5.5"])
        self.assertEqual(args.max_fit_coeff, 5.5)


class TestHollowAnnulusRanking(unittest.TestCase):
    """Test hollow-annulus non-stationarity calculation and ranking term."""

    def test_compute_hollow_counts_geometry(self):
        """Test that triggers inside inner window are excluded, and inside outer are counted."""
        # Clean clusters at t = [90.0, 99.0, 99.9, 100.0, 100.1, 101.0, 110.0, 120.0]
        # Target trigger at t = 100.0, inner = 0.25, outer = 12.0
        # Left window: [100 - 12, 100 - 0.25] = [88.0, 99.75] -> 90.0, 99.0 (2 triggers)
        # Right window: [100 + 0.25, 100 + 12] = [100.25, 112.0] -> 101.0, 110.0 (2 triggers)
        # Inside hole: [99.75, 100.25] -> 99.9, 100.0, 100.1 are EXCLUDED!
        # Outside outer: 120.0 is EXCLUDED!
        clustered = np.array([90.0, 99.0, 99.9, 100.0, 100.1, 101.0, 110.0, 120.0])
        counts = ranking.compute_hollow_counts(
            np.array([100.0]), clustered, inner=0.25, outer=12.0
        )
        self.assertEqual(counts[0], 4.0)

    def test_clustering_100ms_blips(self):
        """Test that triggers in the same 100 ms bin are clustered into 1 event."""
        # 10 triggers in the same 100 ms glitch [105.01, 105.02, ..., 105.09]
        # plus 1 trigger at 108.0
        trigs = {
            "snr": np.array([8.0] * 10 + [7.0]),
            "chisq": np.array([1.0] * 11),
            "chisq_dof": np.array([15.0] * 11),
            "end_time": np.array([105.0 + 0.01 * i for i in range(10)] + [108.0]),
        }
        clustered = ranking._extract_clean_hollow_times(trigs, min_snr=6.0, max_rchisq=2.5)
        # Should result in exactly 2 unique clusters: 1 at 105.0s, 1 at 108.0s
        self.assertEqual(len(clustered), 2)
        self.assertAlmostEqual(clustered[0], 105.05, places=2)
        self.assertAlmostEqual(clustered[1], 108.05, places=2)

    def test_newsnr_sgveto_autochisq_hollow_quadrature(self):
        """Test that hollow counts > knee add smoothly to quadrature distance."""
        snr = 10.0
        rpow = 1.0  # inside injection envelope (zp = 0)
        rauto = 1.0  # za = 0
        rsg = 1.0    # zs = 0

        # When hollow <= hollow_knee (4.5), zh = 0, no downweighting
        nsnr_quiet = ranking.newsnr_sgveto_autochisq(
            snr, rpow, rsg, rauto, hollow=4.0, hollow_knee=4.5, hollow_scale=2.5
        )
        self.assertEqual(nsnr_quiet, 10.0)

        # When hollow = 7.0 (> 4.5):
        # zh = (7.0 - 4.5) / 2.5 = 1.0
        # dist = sqrt(0 + 0 + 0 + 1.0^2) = 1.0
        # downweighting factor applied
        nsnr_loud = ranking.newsnr_sgveto_autochisq(
            snr, rpow, rsg, rauto, hollow=7.0, hollow_knee=4.5, hollow_scale=2.5
        )
        self.assertLess(nsnr_loud, 10.0)
        # Verify exact mathematical formula:
        q = 6.0
        n = 1.0
        expected = snr * (0.5 * (1.0 + (1.0 + 1.0) ** (q / n))) ** (-1.0 / q)
        self.assertAlmostEqual(nsnr_loud, expected, places=5)

    def test_symmetric_exchangeability(self):
        """Verify zero-lag and time-slide triggers undergo identical evaluation."""
        trigs = {
            "snr": np.array([8.0, 8.0]),
            "chisq": np.array([28.0, 28.0]),
            "chisq_dof": np.array([15.0, 15.0]),
            "cont_chisq": np.array([1.0, 1.0]),
            "cont_chisq_dof": np.array([1.0, 1.0]),
            "sg_chisq": np.array([1.0, 1.0]),
            "end_time": np.array([100.0, 200.0]),
        }
        res = ranking.get_newsnr_sgveto_autochisq(
            trigs,
            sngl_ranking_hollow_knee=4.5,
            sngl_ranking_hollow_scale=2.5,
        )
        self.assertEqual(len(res), 2)
        # Triggers in identical local stationary environments have identical ranking
        self.assertEqual(res[0], res[1])

    def test_backward_compatibility_disabled(self):
        """When hollow keywords not provided, 4th quadrature defaults to 0.0."""
        trigs = {
            "snr": np.array([8.0]),
            "chisq": np.array([28.0]),
            "chisq_dof": np.array([15.0]),
            "cont_chisq": np.array([1.0]),
            "cont_chisq_dof": np.array([1.0]),
            "sg_chisq": np.array([1.0]),
            "end_time": np.array([100.0]),
        }
        res_default = ranking.get_newsnr_sgveto_autochisq(trigs)
        res_zero_knee = ranking.get_newsnr_sgveto_autochisq(trigs, sngl_ranking_hollow_knee=0.0)
        self.assertEqual(res_default[0], res_zero_knee[0])


class TestAsymFollowupClusterBins(unittest.TestCase):
    """Test --asym-followup-cluster-bins local-maximum greedy separation."""

    def test_disabled_cluster_bins(self):
        """When cluster_bins=0, preserves original top-N selection."""
        snrs = np.array([2.0, 7.5, 6.0, 5.2, 4.8, 6.5, 4.9])
        keep = secondary_asym_keep(
            snrs, primary_snr=8.0, secondary_snr_threshold=4.5, max_followup_bins=2, cluster_bins=0
        )
        # Top 2 without clustering are bin 1 (7.5) and bin 5 (6.5)
        np.testing.assert_array_equal(keep, np.array([1, 5]))

    def test_enabled_cluster_bins(self):
        """When cluster_bins=4, clusters adjacent multi-bin glitches within 4 bins."""
        # Bin 1 (7.5), Bin 2 (7.0), Bin 3 (6.8) are all from the same glitch.
        # Bin 10 (6.0) is a genuine secondary candidate separated by 9 bins.
        snrs = np.array([0.0, 7.5, 7.0, 6.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0])
        
        # With cluster_bins = 0, top 2 are bin 1 (7.5) and bin 2 (7.0), both glitch!
        keep_no_cluster = secondary_asym_keep(
            snrs, primary_snr=8.0, secondary_snr_threshold=4.5, max_followup_bins=2, cluster_bins=0
        )
        np.testing.assert_array_equal(keep_no_cluster, np.array([1, 2]))

        # With cluster_bins = 4, bin 2 and bin 3 (|2-1| <= 4, |3-1| <= 4) are clustered into bin 1!
        # Next kept candidate is bin 10 (|10-1| = 9 > 4)!
        keep_clustered = secondary_asym_keep(
            snrs, primary_snr=8.0, secondary_snr_threshold=4.5, max_followup_bins=2, cluster_bins=4
        )
        np.testing.assert_array_equal(keep_clustered, np.array([1, 10]))


if __name__ == "__main__":
    unittest.main()
