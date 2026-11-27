"""Unit test for coincident ranking statistic implementations."""

import os
import tempfile
import unittest

import h5py
import numpy as np
from utils import parse_args_cpu_only, simple_exit

from pycbc.events.stat import ExpFitStatistic, PhaseTDStatistic, statistic_dict

# this test only needs to happen on the CPU
parse_args_cpu_only("coinc stats")


class CoincStatTest(unittest.TestCase):
    def setUp(self):
        # one could loop over all single rankings
        # and detector combinations for a complete test
        self.single_rank = "snr"
        self.detectors = ["H1", "L1"]

        # simulate some single-detector triggers from Gaussian noise
        self.num_trigs = n = 100
        self.single_trigs = {}
        for d in self.detectors:
            self.single_trigs[d] = {
                "snr": np.random.chisquare(2, size=n) ** 0.5,
                "coa_phase": np.random.uniform(0, 2 * np.pi, size=n),
                "end_time": 1295441120 + np.random.uniform(-0.01, 0.01, size=n),
                "sigmasq": np.random.uniform(1, 10, size=n),
            }

        self.time_slide_shift = 0.1
        self.time_slide_vector = np.zeros(n, dtype=int)

        # here we should also prepare some files needed by some stats,
        # like the PTA histograms or fake trigger fits
        self.aux_files = []


# dynamically insert a test case method for each available statistic
for stat_name in statistic_dict:
    # FIXME until we can fake the required files,
    # do not test stats that require them
    if "phasetd" in stat_name or "exp_fit" in stat_name or "ExpFit" in stat_name:
        continue

    def stat_test_method(self, sn=stat_name):
        # instantiate an object for this stat
        stat = statistic_dict[sn](
            self.single_rank, files=self.aux_files, ifos=self.detectors
        )

        # get single-detector statistic at each detector from the fake triggers
        single_ranks = [(d, stat.single(self.single_trigs[d])) for d in self.detectors]

        # pretend the fake triggers are coincident, and rank the coincidences
        coinc_ranks = stat.rank_stat_coinc(
            single_ranks,
            self.time_slide_vector,
            self.time_slide_shift,
            (0, 0),  # does this make sense?
        )

        # basic sanity check on the produced ranks
        self.assertEqual(len(coinc_ranks), self.num_trigs)
        self.assertFalse(np.isnan(coinc_ranks).any())

        # here one could add a statistical test based on known analytic
        # behavior of a particular statistic in Gaussian noise

    setattr(CoincStatTest, "test_" + stat_name, stat_test_method)


class PhaseTDSmoothingTest(unittest.TestCase):
    """Test optional 3D mismatch smoothing on PhaseTDStatistic."""

    def test_phasetd_3d_mismatch_smoothing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            hist_path = os.path.join(tmpdir, "phasetd_newsnr_H1L1.hdf")
            pdtype = np.dtype([("c0", np.int16), ("c1", np.int16), ("c2", np.int16)])
            with h5py.File(hist_path, "w") as hf:
                hf.attrs["stat"] = "phasetd_newsnr_H1L1"
                hf.attrs["ifos"] = np.array(["H1", "L1"], dtype=object)
                hf.attrs["twidth"] = 0.001
                hf.attrs["pwidth"] = 2.0 * np.pi / 8.0
                hf.attrs["swidth"] = 0.2
                hf.attrs["srbmin"] = 0
                hf.attrs["srbmax"] = 6
                hf.attrs["sensitivity_ratios"] = np.array([1.0, 1.0])
                for ifo in ("H1", "L1"):
                    grp = hf.create_group(ifo)
                    # Put peak at c1=0 (phase 0), c2=5 (sdif=1.0) and floor at corners
                    grp["param_bin"] = np.array(
                        [(0, 0, 5), (-3, 0, 0), (3, 7, 6)], dtype=pdtype
                    )
                    grp["weights"] = np.array([100.0, 1.0, 1.0], dtype=np.float32)

            with self.assertRaises(ValueError):
                PhaseTDStatistic(
                    "snr",
                    files=[hist_path],
                    ifos=["H1", "L1"],
                    phasetd=True,
                    phasetd_smooth_t=-0.5,
                )

            stat_none = PhaseTDStatistic(
                "snr", files=[hist_path], ifos=["H1", "L1"], phasetd=True
            )
            stat_zero = PhaseTDStatistic(
                "snr",
                files=[hist_path],
                ifos=["H1", "L1"],
                phasetd=True,
                phasetd_smooth_t=0.0,
                phasetd_smooth_p=0.0,
                phasetd_smooth_s=0.0,
            )
            np.testing.assert_array_equal(
                stat_none.two_det_weights["H1"], stat_zero.two_det_weights["H1"]
            )
            self.assertEqual(stat_none.hist_max, stat_zero.hist_max)

            stat_sm = PhaseTDStatistic(
                "snr",
                files=[hist_path],
                ifos=["H1", "L1"],
                phasetd=True,
                phasetd_smooth_t=1.0,
                phasetd_smooth_p=1.2,
                phasetd_smooth_s=0.5,
            )
            w_sm = stat_sm.two_det_weights["H1"]
            c0_mid = stat_sm.c0_size["H1"] // 2
            c1_mid = stat_sm.c1_size["H1"] // 2
            c2_mid = stat_sm.c2_size["H1"] // 2

            # Empty negative c1/c2 padding bins stay at max_penalty
            self.assertAlmostEqual(
                float(w_sm[c0_mid, c1_mid - 1, c2_mid + 5]), float(stat_sm.max_penalty)
            )
            self.assertAlmostEqual(
                float(w_sm[c0_mid, c1_mid, c2_mid - 1]), float(stat_sm.max_penalty)
            )
            # Periodic phase wrap: c1=7 (next to 2pi) receives smoothed weight from c1=0
            self.assertGreater(
                float(w_sm[c0_mid, c1_mid + 7, c2_mid + 5]), float(stat_sm.max_penalty)
            )
            # Peak decreases after smoothing and hist_max tracks the smoothed grid
            self.assertLess(
                float(w_sm.max()), float(stat_none.two_det_weights["H1"].max())
            )
            self.assertAlmostEqual(stat_sm.hist_max, float(w_sm.max()), places=6)

            # Adjacent time-mismatch trigger (dt = 1 bin) scores higher with smoothing
            trigs = {
                "H1": {
                    "snr": np.array([6.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                    "sigmasq": np.array([1.0]),
                },
                "L1": {
                    "snr": np.array([6.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([999.999]),
                    "sigmasq": np.array([1.0]),
                },
            }
            r_none = stat_none.rank_stat_coinc(
                [(d, stat_none.single(trigs[d])) for d in ("H1", "L1")],
                np.array([0]),
                0.1,
                (0, 0),
            )
            r_sm = stat_sm.rank_stat_coinc(
                [(d, stat_sm.single(trigs[d])) for d in ("H1", "L1")],
                np.array([0]),
                0.1,
                (0, 0),
            )
            self.assertGreater(float(r_sm[0]), float(r_none[0]))


class ExpFitFeaturesTest(unittest.TestCase):
    """Test PR 4 modular ExpFitStatistic feature plugins."""

    def _write_fit_file(self, path, ifo):
        with h5py.File(path, "w") as hf:
            hf.attrs["stat"] = f"{ifo}-fit_coeffs"
            hf.attrs["stat_threshold"] = 5.5
            hf.attrs["analysis_time"] = 1000.0
            hf["template_id"] = np.array([1, 0], dtype=np.int32)
            hf["fit_coeff"] = np.array([6.0, 6.0], dtype=np.float32)
            hf["count_above_thresh"] = np.array([10.0, 10.0], dtype=np.float32)
            hf["count_in_template"] = np.array([100.0, 100.0], dtype=np.float32)
            hf["template_duration"] = np.array([0.2, 5.0], dtype=np.float32)
            hf["chi_eff"] = np.array([0.7, 0.0], dtype=np.float32)
            hf["eta"] = np.array([0.10, 0.25], dtype=np.float32)
            fbt = hf.create_group("fit_by_template")
            # Template 0: low count (noisy); Template 1: high count glitchy tail
            fbt["fit_coeff"] = np.array([2.0, 3.0], dtype=np.float32)
            fbt["count_above_thresh"] = np.array([90.0, 1.0], dtype=np.float32)
            fbt["count_in_template"] = np.array([200.0, 50.0], dtype=np.float32)

    def test_empirical_bayes_and_priors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            h1_fit = os.path.join(tmpdir, "H1-fit.hdf")
            l1_fit = os.path.join(tmpdir, "L1-fit.hdf")
            self._write_fit_file(h1_fit, "H1")
            self._write_fit_file(l1_fit, "L1")

            base = ExpFitStatistic(
                "snr",
                files=[h1_fit, l1_fit],
                ifos=["H1", "L1"],
                phasetd=False,
                sensitive_volume=False,
                normalize_fit_rate=False,
                kde=False,
                dq=False,
                chirp_mass=False,
            )
            np.testing.assert_allclose(
                base.fits_by_tid["H1"]["smoothed_fit_coeff"], [6.0, 6.0]
            )

            eb = ExpFitStatistic(
                "snr",
                files=[h1_fit, l1_fit],
                ifos=["H1", "L1"],
                phasetd=False,
                sensitive_volume=False,
                normalize_fit_rate=False,
                kde=False,
                dq=False,
                chirp_mass=False,
                empirical_bayes_fit=True,
                eb_alpha_prior_count=10.0,
                eb_rate_prior_count=10.0,
                template_prior=True,
                prior_short_dur_weight=1.5,
                prior_dur_knee=0.5,
                prior_dur_scale=0.25,
                prior_chi_eff_coeff=2.0,
                prior_chi_eff_knee=0.3,
                prior_eta_coeff=10.0,
                prior_eta_knee=0.16,
                network_l2_newsnr=True,
                l2_newsnr_coeff=1.0,
                l2_newsnr_min=0.0,
                l2_newsnr_max=20.0,
                l2_newsnr_ref=10.0,
                l2_l1_excess_coeff=0.5,
                l2_l1_excess_max=2.0,
            )
            # Template 0 (fbt_count=1) stays near prior; template 1 shrinks toward 2.0
            a_post = eb.fits_by_tid["H1"]["smoothed_fit_coeff"]
            self.assertGreater(a_post[0], 5.0)
            self.assertLess(a_post[1], 2.5)
            # Template 1 (short dur, high spin, low eta) gets negative log-prior
            self.assertAlmostEqual(float(eb.template_prior[0]), 0.0, places=6)
            self.assertLess(float(eb.template_prior[1]), -2.0)

            trigs_eq = {
                "H1": {
                    "template_id": np.array([0]),
                    "snr": np.array([7.5]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
                "L1": {
                    "template_id": np.array([0]),
                    "snr": np.array([7.5]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
            }
            trigs_asym = {
                "H1": {
                    "template_id": np.array([0]),
                    "snr": np.array([10.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
                "L1": {
                    "template_id": np.array([0]),
                    "snr": np.array([5.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
            }
            s_eq = [(d, eb.single(trigs_eq[d])) for d in ("H1", "L1")]
            s_asym = [(d, eb.single(trigs_asym[d])) for d in ("H1", "L1")]
            stat_eq = float(
                eb.rank_stat_coinc(
                    s_eq, np.array([0]), 0.1, (0, 0), time_addition=0.01
                )[0]
            )
            stat_asym = float(
                eb.rank_stat_coinc(
                    s_asym, np.array([0]), 0.1, (0, 0), time_addition=0.01
                )[0]
            )
            self.assertGreater(stat_asym, stat_eq)

            # coinc_lim_for_thresh handles both full s and 1-element fixed_single_info
            trigs_loud = {
                "H1": {
                    "template_id": np.array([0]),
                    "snr": np.array([12.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
                "L1": {
                    "template_id": np.array([0]),
                    "snr": np.array([12.0]),
                    "coa_phase": np.array([0.0]),
                    "end_time": np.array([1000.0]),
                },
            }
            s_loud = [(d, eb.single(trigs_loud[d])) for d in ("H1", "L1")]
            stat_loud = float(
                eb.rank_stat_coinc(
                    s_loud, np.array([0]), 0.1, (0, 0), time_addition=0.01
                )[0]
            )
            lim_fixed = float(
                eb.coinc_lim_for_thresh(
                    [s_loud[1]], stat_loud - 0.1, "H1", time_addition=0.01
                )[0]
            )
            self.assertGreaterEqual(lim_fixed, float(s_loud[0][1]["snglstat"][0]))

    def test_kde_sorted_and_safe(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            h1_fit = os.path.join(tmpdir, "H1-fit.hdf")
            l1_fit = os.path.join(tmpdir, "L1-fit.hdf")
            sig_kde = os.path.join(tmpdir, "sig-kde.hdf")
            tmp_kde = os.path.join(tmpdir, "tmp-kde.hdf")
            self._write_fit_file(h1_fit, "H1")
            self._write_fit_file(l1_fit, "L1")
            with h5py.File(sig_kde, "w") as hf:
                hf.attrs["stat"] = "signal-kde_file"
                hf["template_id"] = np.array([1, 0], dtype=np.int32)
                hf["data_kde"] = np.array([0.1, 2.0], dtype=np.float64)
            with h5py.File(tmp_kde, "w") as hf:
                hf.attrs["stat"] = "template-kde_file"
                hf["template_id"] = np.array([1, 0], dtype=np.int32)
                hf["data_kde"] = np.array([0.0, 1.0], dtype=np.float64)

            kde_stat = ExpFitStatistic(
                "snr",
                files=[h1_fit, l1_fit, sig_kde, tmp_kde],
                ifos=["H1", "L1"],
                phasetd=False,
                sensitive_volume=False,
                normalize_fit_rate=False,
                kde=True,
                dq=False,
                chirp_mass=False,
            )
            kde_stat.curr_tnum = np.array([0, 1])
            ratios = kde_stat.kde_ratio()
            self.assertTrue(np.all(np.isfinite(ratios)))
            self.assertGreater(float(ratios[1]), float(ratios[0]))


# create and populate unittest's test suite
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(CoincStatTest))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(PhaseTDSmoothingTest))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(ExpFitFeaturesTest))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
