"""Unit test for coincident ranking statistic implementations."""

import os
import tempfile
import unittest

import h5py
import numpy as np
from utils import parse_args_cpu_only, simple_exit

from pycbc.events.stat import PhaseTDStatistic, statistic_dict

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


# create and populate unittest's test suite
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(CoincStatTest))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(PhaseTDSmoothingTest))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
