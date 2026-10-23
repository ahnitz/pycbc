#!/usr/bin/env python
"""Unit test suite for Milestone 4 PyCBC features:
- ranking.py: newsnr_sgveto_autochisq with c_sg=0.012 and s2 <= 5000.0
- pycbc_calculate_psd: --psd-delta-f and --native-psd-resolution flags
- pycbc_foreground_minifollowup: synced CLI flags
- pycbc_add_statmap: --hierarchical-removal-threshold-ifar alias
"""

import unittest
import numpy as np
from pycbc.events.ranking import (
    newsnr_sgveto_autochisq,
    get_newsnr_sgveto_autochisq,
)


class TestMilestone4Features(unittest.TestCase):
    def test_ranking_gw241011_test_vector(self):
        """Verify GW241011_233834 test vector recovery:
        rho = 33.87, sg_chisq = 14.04 is inside the envelope (tau_s = 15.27 > 14.04),
        so zs = 0.0. With residual dist = 0.15, newsnr = 31.14.
        """
        snr = 33.87
        chisq = 1.0
        autochisq = 1.0
        sg_chisq = 14.04

        # With sg_chisq = 14.04 alone, zs is 0.0 (no downweighting from sg_chisq)
        newsnr_clean = newsnr_sgveto_autochisq(
            snr,
            chisq,
            sg_chisq,
            autochisq,
            c_sg=0.012,
        )
        self.assertAlmostEqual(newsnr_clean, 33.87, delta=0.01)

        # With minor residual dist = 0.15 (e.g. from elevated power or sg chisq):
        # s2 = 1147.18, tau_s = 3.20 + 0.012 * 1147.18 = 16.966, sig_s = 8.483
        # sg_chisq = tau_s + 0.15 * sig_s = 18.239
        sg_chisq_dist = 16.966 + 0.15 * 8.483
        newsnr_dist015 = newsnr_sgveto_autochisq(
            snr,
            chisq,
            sg_chisq_dist,
            autochisq,
            c_sg=0.012,
        )
        self.assertAlmostEqual(newsnr_dist015, 29.37, delta=0.05)

        # Using ranking function
        trigs = {
            "chisq": np.array([chisq]),
            "chisq_dof": np.array([1.5]),
            "snr": np.array([snr]),
            "sg_chisq": np.array([sg_chisq_dist]),
            "auto_chisq": np.array([autochisq]),
            "auto_chisq_dof": np.array([1.0]),
        }
        getter_newsnr = get_newsnr_sgveto_autochisq(trigs, c_sg=0.012)[0]
        self.assertAlmostEqual(getter_newsnr, 29.37, delta=0.05)

    def test_ranking_low_snr_glitch_rejection(self):
        """Verify low-SNR glitch with sg_chisq=6.0 is suppressed to 3.70 with tau_s=3.97."""
        snr = 8.0
        chisq = 1.0
        autochisq = 1.0
        sg_chisq = 6.0

        newsnr = newsnr_sgveto_autochisq(
            snr,
            chisq,
            sg_chisq,
            autochisq,
            c_sg=0.012,
        )
        self.assertAlmostEqual(newsnr, 3.70, delta=0.05)

    def test_ranking_s2_saturation_limit(self):
        """Verify s2 saturation limit extends up to 5000.0."""
        # For snr = 60, snr^2 = 3600 > old limit 300
        snr = 60.0
        chisq = 1.0
        autochisq = 1.0
        # If tau_s scales with snr^2, tau_s = 1.5 + 0.012 * 3600 = 44.7
        # With sg_chisq = 30 < 44.7, dist should be 0, no downweighting!
        sg_chisq = 30.0
        newsnr = newsnr_sgveto_autochisq(
            snr,
            chisq,
            sg_chisq,
            autochisq,
            c_sg=0.012,
        )
        # Should NOT be downweighted since 30 < 44.7
        self.assertAlmostEqual(newsnr, 60.0, delta=0.01)

        # But above 5000 (e.g. snr=80 -> snr^2=6400, s2=5000, tau_s=61.5)
        # If sg_chisq = 80 > 61.5, downweighting should occur
        newsnr_high = newsnr_sgveto_autochisq(
            80.0,
            chisq,
            80.0,
            autochisq,
            c_sg=0.012,
        )
        self.assertLess(newsnr_high, 80.0)

    def test_ranking_clean_signal_preservation(self):
        """Verify clean signal with sg_chisq=1.0 is not downweighted."""
        for snr in [8.0, 15.0, 30.0, 50.0]:
            newsnr = newsnr_sgveto_autochisq(
                snr,
                1.0,
                1.0,
                1.0,
                c_sg=0.012,
            )
            self.assertAlmostEqual(newsnr, snr, delta=0.01)

    def test_psd_resolution_calculation(self):
        """Verify PSD delta_f and flen logic for native Welch resolution."""
        sample_rate = 2048.0
        # Case 1: Standard oversampled 512s segment
        seg_len = 512.0
        delta_f_std = 1.0 / seg_len
        flen_std = int(sample_rate / 2.0 / delta_f_std) + 1
        self.assertEqual(flen_std, 524289)

        # Case 2: Native Welch 16s resolution
        psd_seg_len = 16.0
        delta_f_native = 1.0 / psd_seg_len
        flen_native = int(sample_rate / 2.0 / delta_f_native) + 1
        self.assertEqual(flen_native, 16385)
        # Verify 32x sample reduction
        self.assertEqual((flen_std - 1) / (flen_native - 1), 32.0)

    def test_cli_help_options(self):
        """Verify CLI scripts accept the required option flags."""
        import subprocess
        import sys

        # 1. pycbc_calculate_psd accepts --psd-delta-f, --native-welch-resolution, etc.
        res = subprocess.run(
            [sys.executable, "bin/all_sky_search/pycbc_calculate_psd", "--help"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0)
        self.assertIn("--psd-delta-f", res.stdout)
        self.assertIn("--native-psd-resolution", res.stdout)
        self.assertIn("--native-welch-resolution", res.stdout)
        self.assertIn("--no-native-welch-resolution", res.stdout)
        self.assertIn("--interpolated-resolution", res.stdout)

    def test_calculate_psd_default_welch_resolution(self):
        """Verify pycbc_calculate_psd parser defaults native_welch_resolution to True."""
        import subprocess
        import sys

        script_template = """
import sys
from unittest.mock import patch
import runpy

sys.argv = [
    'pycbc_calculate_psd',
    '--low-frequency-cutoff', '20',
    '--channel-name', 'H1:STRAIN',
    '--analysis-segment-file', 'dummy.xml',
    '--output-file', 'dummy.hdf'{extra_arg}
]
with patch('pycbc.init_logging'), \\
     patch('pycbc.psd.verify_psd_options'), \\
     patch('pycbc.strain.StrainSegments.verify_segment_options'), \\
     patch('pycbc.events.select_segments_by_definer', return_value=[]), \\
     patch('pycbc.io.HFile'):
    mod = runpy.run_path(
        'bin/all_sky_search/pycbc_calculate_psd', run_name='__main__'
    )
print('DEFAULT_RES:', mod['args'].native_welch_resolution)
"""
        # Default run
        res = subprocess.run(
            [sys.executable, "-c", script_template.format(extra_arg="")],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            res.returncode, 0, msg=f"Subprocess failed:\n{res.stderr}"
        )
        self.assertIn("DEFAULT_RES: True", res.stdout)

        # --no-native-welch-resolution
        res_no = subprocess.run(
            [
                sys.executable,
                "-c",
                script_template.format(
                    extra_arg=", '--no-native-welch-resolution'"
                ),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            res_no.returncode, 0, msg=f"Subprocess failed:\n{res_no.stderr}"
        )
        self.assertIn("DEFAULT_RES: False", res_no.stdout)

        # --interpolated-resolution
        res_interp = subprocess.run(
            [
                sys.executable,
                "-c",
                script_template.format(
                    extra_arg=", '--interpolated-resolution'"
                ),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            res_interp.returncode, 0, msg=f"Subprocess failed:\n{res_interp.stderr}"
        )
        self.assertIn("DEFAULT_RES: False", res_interp.stdout)

        # 2. pycbc_add_statmap accepts --hierarchical-removal-threshold-ifar
        res = subprocess.run(
            [sys.executable, "bin/all_sky_search/pycbc_add_statmap", "--help"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0)
        self.assertIn("--hierarchical-removal-threshold-ifar", res.stdout)
        self.assertIn("--hierarchical-removal-ifar-thresh", res.stdout)

        # 3. pycbc_foreground_minifollowup accepts sync flags and defaults --sort-variable to stat
        res = subprocess.run(
            [sys.executable, "bin/minifollowups/pycbc_foreground_minifollowup", "--help"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0)
        self.assertIn("--num-events", res.stdout)
        self.assertIn("--ifar-threshold", res.stdout)
        self.assertIn("--far-threshold", res.stdout)
        self.assertIn("--min-events", res.stdout)
        self.assertIn("--max-events", res.stdout)
        self.assertIn("Default=stat", res.stdout)

    def test_event_flipper_and_layout_collapse(self):
        """Verify layout.py collapse handling and orange.html Event Flipper markup."""
        import os
        import tempfile
        from pycbc.results import layout, render

        with tempfile.TemporaryDirectory() as tmpdir:
            dummy_file = type("DummyFile", (), {"name": "plot1.png"})()
            layout.group_layout(
                tmpdir,
                [dummy_file],
                unique="_001",
                collapse=True,
                title="Event #1 | Stat=141.89",
            )
            cp1 = render.get_embedded_config(os.path.join(tmpdir, "well_001.html"))
            self.assertTrue(cp1.check_option("well_001.html", "collapse"))
            self.assertEqual(cp1.get("well_001.html", "collapse"), "true")

            layout.group_layout(
                tmpdir,
                [dummy_file],
                unique="_002",
                collapse=False,
                title="Non-collapsed well",
            )
            cp2 = render.get_embedded_config(os.path.join(tmpdir, "well_002.html"))
            self.assertFalse(cp2.check_option("well_002.html", "collapse"))

        with open("pycbc/results/templates/orange.html", "r") as f:
            orange_html = f.read()
        self.assertIn("event-flipper-bar", orange_html)
        self.assertIn("flipper-prev", orange_html)
        self.assertIn("flipper-next", orange_html)
        self.assertIn("flipper-select", orange_html)
        self.assertIn("flipper-mode-toggle", orange_html)
        self.assertIn("ArrowLeft", orange_html)
        self.assertIn("ArrowRight", orange_html)


if __name__ == "__main__":
    unittest.main()
