"""Unit tests for PR 5 workflow, statmap, plotting, and HTML scalability fixes."""

import os
import subprocess
import sys
import tempfile
import unittest
from configparser import ConfigParser

import h5py
import numpy as np
from jinja2 import Environment, FileSystemLoader

import pycbc.results

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


class PR5ScalabilityTest(unittest.TestCase):
    """Focused tests for PR 5 bugfixes and HDF5/HTML bloat elimination."""

    def test_plot_gating_multi_ifo(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gate_hdf = os.path.join(tmpdir, "gates.hdf")
            out_html = os.path.join(tmpdir, "gates.html")
            with h5py.File(gate_hdf, "w") as hf:
                for ifo, t_off in (("H1", 1000.0), ("L1", 1100.0)):
                    grp = hf.create_group(f"{ifo}/gating/auto")
                    grp["time"] = np.array([t_off], dtype=np.float64)
                    grp["width"] = np.array([1.0], dtype=np.float64)
                    grp["pad"] = np.array([0.5], dtype=np.float64)

            script = os.path.join(REPO_ROOT, "bin", "plotting", "pycbc_plot_gating")
            env = os.environ.copy()
            env["PYTHONPATH"] = REPO_ROOT
            subprocess.run(
                [
                    sys.executable,
                    script,
                    "--input-file",
                    gate_hdf,
                    "--output-file",
                    out_html,
                ],
                check=True,
                env=env,
            )
            with open(out_html, encoding="utf-8") as fp:
                html_text = fp.read()
            self.assertIn("H1 auto gates:", html_text)
            self.assertIn("L1 auto gates:", html_text)

    def test_page_snrchi_zero_chisq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            trig_hdf = os.path.join(tmpdir, "trigs.hdf")
            out_png = os.path.join(tmpdir, "snrchi.png")
            with h5py.File(trig_hdf, "w") as hf:
                grp = hf.create_group("H1")
                grp["snr"] = np.array([6.5, 7.5, 9.0, 8.0], dtype=np.float32)
                grp["chisq"] = np.array([0.0, 14.0, 28.0, 0.0], dtype=np.float32)
                grp["chisq_dof"] = np.array([8, 8, 8, 0], dtype=np.int32)
                grp["end_time"] = np.array(
                    [1000.0, 1001.0, 1002.0, 1003.0], dtype=np.float64
                )
                hf["H1/segments/start"] = np.array([990.0], dtype=np.float64)
                hf["H1/segments/end"] = np.array([1010.0], dtype=np.float64)

            script = os.path.join(REPO_ROOT, "bin", "plotting", "pycbc_page_snrchi")
            env = os.environ.copy()
            env["PYTHONPATH"] = REPO_ROOT
            subprocess.run(
                [
                    sys.executable,
                    script,
                    "--min-snr",
                    "6.0",
                    "--chisq-choice",
                    "traditional",
                    "--newsnr-contours",
                    "6",
                    "--trigger-file",
                    trig_hdf,
                    "--ifo",
                    "H1",
                    "--output-file",
                    out_png,
                ],
                check=True,
                env=env,
            )
            self.assertTrue(os.path.exists(out_png))
            self.assertGreater(os.path.getsize(out_png), 0)

    def test_coinc_statmap_no_intermediate_h_bloat(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            coinc_hdf = os.path.join(tmpdir, "coinc.hdf")
            statmap_hdf = os.path.join(tmpdir, "statmap.hdf")
            n_bg = 20
            with h5py.File(coinc_hdf, "w") as hf:
                hf.attrs["ifos"] = "H1 L1"
                hf.attrs["num_of_ifos"] = 2
                hf.attrs["pivot"] = "H1"
                hf.attrs["fixed"] = "L1"
                hf.attrs["H1_foreground_time"] = 100.0
                hf.attrs["L1_foreground_time"] = 100.0
                hf.attrs["coinc_time"] = 100.0
                hf.attrs["timeslide_interval"] = 5.0

                stats = np.concatenate(
                    [np.array([25.0, 20.0]), np.linspace(1.0, 5.0, n_bg)]
                )
                tsid = np.concatenate([np.array([0, 0]), np.arange(1, n_bg + 1)])
                t_h1 = np.concatenate(
                    [np.array([1010.0, 1050.0]), np.linspace(1000.0, 1090.0, n_bg)]
                )
                t_l1 = np.concatenate(
                    [np.array([1010.0, 1050.0]), np.linspace(1000.0, 1090.0, n_bg)]
                )

                hf["stat"] = stats.astype(np.float32)
                hf["timeslide_id"] = tsid.astype(np.int32)
                hf["template_id"] = np.zeros(len(stats), dtype=np.int32)
                hf["decimation_factor"] = np.ones(len(stats), dtype=np.int32)
                hf["H1/time"] = t_h1
                hf["H1/trigger_id"] = np.arange(len(stats), dtype=np.int64)
                hf["L1/time"] = t_l1
                hf["L1/trigger_id"] = np.arange(len(stats), dtype=np.int64)
                for key in ("H1", "L1", "H1L1"):
                    hf[f"segments/{key}/start"] = np.array([1000.0], dtype=np.float64)
                    hf[f"segments/{key}/end"] = np.array([1100.0], dtype=np.float64)

            script = os.path.join(
                REPO_ROOT, "bin", "all_sky_search", "pycbc_coinc_statmap"
            )
            env = os.environ.copy()
            env["PYTHONPATH"] = REPO_ROOT
            subprocess.run(
                [
                    sys.executable,
                    script,
                    "--coinc-files",
                    coinc_hdf,
                    "--output-file",
                    statmap_hdf,
                    "--ifos",
                    "H1",
                    "L1",
                    "--cluster-window",
                    "1.0",
                    "--max-hierarchical-removal",
                    "-1",
                    "--hierarchical-removal-against",
                    "exclusive",
                ],
                check=True,
                env=env,
            )
            with h5py.File(statmap_hdf, "r") as hf:
                self.assertEqual(int(hf.attrs["hierarchical_removal_iterations"]), 2)
                self.assertIn("background_h0", hf)
                self.assertNotIn("background_h1", hf)
                self.assertIn("background_h2", hf)
                self.assertIn("background_h2/decimation_factor", hf)
                self.assertIn("foreground_h1", hf)
                self.assertIn("foreground_h2", hf)

    def test_orange_template_lazy_media(self):
        tmpl_dir = os.path.join(os.path.dirname(pycbc.results.__file__), "templates")
        env = Environment(loader=FileSystemLoader(tmpl_dir))
        tmpl = env.get_template("orange.html")

        class _FakeFile:
            def __init__(self, fname, html_body):
                self._fname = fname
                self.path = fname
                self._body = html_body

            def filename(self):
                return self._fname

            def render(self):
                return self._body

        class _FakeDir:
            path = "sec"
            files = [
                _FakeFile("well_1.html", ""),
                _FakeFile("well_2.html", ""),
                _FakeFile("extra.html", '<img src="extra.png" />'),
            ]

            def title(self):
                return "Section"

        class _FakeConfigParser(ConfigParser):
            def check_option(self, section, name):
                return self.has_section(section) and self.has_option(section, name)

        def _fake_cp(_path):
            cp = _FakeConfigParser()
            fn = os.path.basename(_path)
            cp.add_section(fn)
            cp.set(fn, "collapse", "True")
            cp.set(fn, "title", fn)
            return cp

        def _fake_render(path, _arg):
            return f'<img src="{os.path.basename(path)}.png" />'

        env.globals.update(
            path_exists=lambda _p: True,
            get_embedded_config=_fake_cp,
            setup_template_render=_fake_render,
        )
        rendered = tmpl.render(
            dir=_FakeDir(),
            dirs=[_FakeDir()],
            plots_dir="",
            dot_dot_str="",
            analysis_title="Test",
            analysis_subtitle="Sub",
            version_info=[],
            get_embedded_config=_fake_cp,
            path_exists=lambda _p: True,
            setup_template_render=_fake_render,
        )
        self.assertIn('<img src="well_1.html.png"', rendered)
        self.assertIn('<img data-src="well_2.html.png"', rendered)
        self.assertIn('<img data-src="extra.png"', rendered)


if __name__ == "__main__":
    unittest.main()
