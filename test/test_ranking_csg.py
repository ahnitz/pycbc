"""Unit tests for ranking c_sg tolerance and minifollowup batching."""

import unittest
import numpy
from pycbc.events import ranking
from pycbc.workflow.plotting import PlotExecutable
from pycbc.workflow.configuration import WorkflowConfigParser


class TestRankingCsg(unittest.TestCase):
    def test_newsnr_sgveto_default_csg(self):
        # High SNR test vector (GW241011-like, SNR ~ 33.87)
        # s2 = min(33.87^2, 5000) = 1147.1769
        # tau_s = 3.20 + 0.012 * 1147.1769 = 16.9661
        snr = numpy.array([33.87])
        rchisq = numpy.array([1.0])
        autochisq = numpy.array([1.0])

        # sgchisq within tau_s should not downweight
        sgchisq_in = numpy.array([16.0])
        res_in = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq_in, autochisq=autochisq
        )
        self.assertAlmostEqual(res_in[0], 33.87, places=4)

        # sgchisq well above tau_s should smoothly downweight
        sgchisq_out = numpy.array([35.0])
        res_out = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq_out, autochisq=autochisq
        )
        self.assertLess(res_out[0], 33.87)
        self.assertGreater(res_out[0], 0.0)

    def test_moderate_snr_intercept_protection(self):
        # Moderate SNR trigger (e.g. GW240512-like, SNR ~ 8.13, sg_chisq ~ 3.0)
        # With base intercept 3.20: tau_s = 3.20 + 0.012 * 8.13^2 = 3.99 > 3.0 -> zero downweighting!
        # With intercept 1.50: tau_s was 1.50 + 0.012 * 66.1 = 2.29 < 3.0 -> was heavily downweighted!
        snr = numpy.array([8.13])
        rchisq = numpy.array([1.0])
        autochisq = numpy.array([1.0])
        sgchisq = numpy.array([3.0])
        res = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=autochisq
        )
        self.assertAlmostEqual(res[0], 8.13, places=4)

        # At near-zero SNR, intercept is strictly 3.20
        snr_low = numpy.array([0.001])
        res_in = ranking.newsnr_sgveto_autochisq(
            snr_low, rchisq, sgchisq=numpy.array([3.19]), autochisq=autochisq
        )
        self.assertAlmostEqual(res_in[0], 0.001, places=4)
        res_out = ranking.newsnr_sgveto_autochisq(
            snr_low, rchisq, sgchisq=numpy.array([3.25]), autochisq=autochisq
        )
        self.assertLess(res_out[0], 0.001)

    def test_newsnr_sgveto_custom_csg(self):
        snr = numpy.array([20.0])
        rchisq = numpy.array([1.0])
        autochisq = numpy.array([1.0])
        sgchisq = numpy.array([10.0])

        # With c_sg=0.012: tau_s = 3.20 + 0.012 * 400 = 8.00 -> downweighted
        res_c012 = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=autochisq, c_sg=0.012
        )
        # With c_sg=0.020: tau_s = 3.20 + 0.020 * 400 = 11.20 -> inside
        res_c020 = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=autochisq, c_sg=0.020
        )
        self.assertLess(res_c012[0], 20.0)
        self.assertAlmostEqual(res_c020[0], 20.0, places=4)

    def test_tau_p_strict_cap_at_300(self):
        # At high SNR (e.g. SNR = 30), s2_pa is capped at 300.0
        # tau_p = 1.50 + 0.0018 * 300.0 = 2.04 (NOT 10.5 or 3.12)
        # sig_p = 0.45 + 0.0008 * 300.0 = 0.69
        snr = numpy.array([30.0])
        autochisq = numpy.array([1.0])
        sgchisq = numpy.array([1.0])

        # rchisq = 2.04 is within tau_p (2.04), zero downweighting
        res_in = ranking.newsnr_sgveto_autochisq(
            snr, brchisq=numpy.array([2.04]), sgchisq=sgchisq, autochisq=autochisq
        )
        self.assertAlmostEqual(res_in[0], 30.0, places=4)

        # rchisq = 2.50 (e.g. L1 glitch at SNR ~ 26, rchisq ~ 2.42) is outside tau_p (2.04)
        # If tau_p were 10.5 (or 3.12), rchisq = 2.50 would NOT be downweighted.
        # With s2_pa capped at 300.0, it is strictly downweighted.
        res_out = ranking.newsnr_sgveto_autochisq(
            snr, brchisq=numpy.array([2.50]), sgchisq=sgchisq, autochisq=autochisq
        )
        self.assertLess(res_out[0], 30.0)

        # At extreme SNR = 100 (where s2_sg reaches 5000), tau_p remains strictly capped at 2.04, NOT 10.5
        snr_100 = numpy.array([100.0])
        res_extreme = ranking.newsnr_sgveto_autochisq(
            snr_100, brchisq=numpy.array([4.90]), sgchisq=sgchisq, autochisq=autochisq
        )
        # If tau_p were 10.5 (from 5000 cap), rchisq = 4.90 would NOT be downweighted!
        # With 300 cap, it is strongly downweighted.
        self.assertLess(res_extreme[0], 100.0)

    def test_tau_a_strict_cap_at_300(self):
        # At high SNR (e.g. SNR = 30), s2_pa is capped at 300.0
        # tau_a = 1.85 + 0.0028 * 300.0 = 2.69 (NOT 15.8 or 4.37)
        # sig_a = 0.65 + 0.0012 * 300.0 = 1.01
        snr = numpy.array([30.0])
        rchisq = numpy.array([1.0])
        sgchisq = numpy.array([1.0])

        # autochisq = 2.69 is within tau_a (2.69), zero downweighting
        res_in = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=numpy.array([2.69])
        )
        self.assertAlmostEqual(res_in[0], 30.0, places=4)

        # autochisq = 3.50 is outside tau_a (2.69) -> strictly downweighted
        # If tau_a were 15.8 (or 4.37), autochisq = 3.50 would NOT be downweighted.
        res_out = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=numpy.array([3.50])
        )
        self.assertLess(res_out[0], 30.0)

        # At extreme SNR = 100, tau_a remains capped at 2.69, NOT 15.8
        snr_100 = numpy.array([100.0])
        res_extreme = ranking.newsnr_sgveto_autochisq(
            snr_100, rchisq, sgchisq=sgchisq, autochisq=numpy.array([5.0])
        )
        # If tau_a were 15.8, autochisq = 5.0 would not be downweighted!
        self.assertLess(res_extreme[0], 100.0)

    def test_tau_s_scales_with_5000_cap(self):
        # Extreme SNR trigger: SNR = 100 -> SNR^2 = 10000, capped at s2_sg = 5000.0
        snr = numpy.array([100.0])
        rchisq = numpy.array([1.0])
        autochisq = numpy.array([1.0])
        # At s2_sg = 5000.0:
        # tau_s = 3.20 + 0.012 * 5000.0 = 63.2
        # sig_s = 1.20 + 0.0035 * 5000.0 = 18.7
        # Inside envelope (sgchisq = 63.0 <= 63.2): zero downweighting
        sgchisq_in = numpy.array([63.0])
        res_in = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq_in, autochisq=autochisq
        )
        self.assertAlmostEqual(res_in[0], 100.0, places=3)

        # Outside envelope (sgchisq = 65.0 > 63.2): downweighted
        sgchisq_out = numpy.array([65.0])
        res_out = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq_out, autochisq=autochisq
        )
        self.assertLess(res_out[0], 100.0)

        # If tau_s were capped at 300.0 like tau_p/tau_a:
        # tau_s would be 3.20 + 0.012 * 300 = 6.80, so sgchisq = 63.0 would be severely downweighted.
        # This confirms tau_s scales with 5000.0 cap while tau_p/tau_a are capped at 300.0.

    def test_get_newsnr_sgveto_autochisq_kwargs(self):
        # Mock trigger dict
        trigs = {
            "chisq_dof": numpy.array([4, 4]),
            "chisq": numpy.array([6.0, 6.0]),
            "sg_chisq": numpy.array([10.0, 10.0]),
            "auto_chisq": numpy.array([1.0, 1.0]),
            "snr": numpy.array([25.0, 25.0]),
        }
        res_default = ranking.get_newsnr_sgveto_autochisq(trigs)
        res_custom = ranking.get_newsnr_sgveto_autochisq(
            trigs, sngl_ranking_c_sg=0.025
        )
        # Higher c_sg gives wider tolerance -> higher or equal newsnr
        self.assertGreaterEqual(res_custom[0], res_default[0])


class TestMinifollowupClustering(unittest.TestCase):
    def test_plot_executable_clustering_profile(self):
        cp = WorkflowConfigParser()
        cp.add_section("workflow")
        cp.add_section("executables")
        cp.set("executables", "page_coincinfo", "/bin/true")
        cp.add_section("workflow-minifollowups")
        cp.set("workflow-minifollowups", "batch-size", "15")

        exe = PlotExecutable(cp, "page_coincinfo")
        self.assertTrue(hasattr(exe, "group_jobs"))
        self.assertEqual(exe.group_jobs, "15")
        self.assertEqual(exe.profiles.get(("pegasus", "clusters.size")), "15")


if __name__ == "__main__":
    unittest.main()
