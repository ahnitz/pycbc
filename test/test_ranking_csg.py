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

    def test_s2_cap_at_5000(self):
        # Extreme SNR trigger: SNR = 100 -> SNR^2 = 10000, capped at 5000
        snr = numpy.array([100.0])
        rchisq = numpy.array([1.0])
        autochisq = numpy.array([1.0])
        # At s2 = 5000, tau_s = 3.20 + 0.012 * 5000 = 63.2
        # sig_s = 1.20 + 0.0035 * 5000 = 18.7
        sgchisq = numpy.array([62.0])
        res = ranking.newsnr_sgveto_autochisq(
            snr, rchisq, sgchisq=sgchisq, autochisq=autochisq
        )
        self.assertAlmostEqual(res[0], 100.0, places=3)

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
