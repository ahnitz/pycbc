# Copyright (C) 2026  Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""Unit tests for the unified DistMarg / SingleTemplate marginalization engine
and TimeDomainFilterBank-backed FIR SNR reconstruction.
"""

import unittest
import astropy.utils.iers
astropy.utils.iers.conf.auto_download = False

import numpy as np
from scipy.special import i0e

from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
    UniformRadius,
)
from pycbc.filter.matched_filter_ratio import reconstruct_fir_snr_span
from pycbc.inference.models.single_template import SingleTemplate
from pycbc.inference.models.tools import setup_distance_marg_interpolant
from pycbc.types import TimeSeries


class TestMarginalLogLR(unittest.TestCase):
    def test_distance_marg_interpolant_phase(self):
        """Verify setup_distance_marg_interpolant matches direct quadrature."""
        dist_locs = np.linspace(10.0, 1000.0, 1000)
        weights = dist_locs ** 2
        weights /= weights.sum()
        dist_rescale = 500.0 / dist_locs

        lut_phase = setup_distance_marg_interpolant(
            (dist_rescale, weights),
            phase=True,
            snr_range=(0.05, 100.0),
            density=(128, 128),
        )

        for sh, hh in [(12.0, 64.0), (25.0, 100.0), (4.5, 36.0)]:
            s = sh * dist_rescale
            h = hh * (dist_rescale ** 2)
            exact = float(
                np.max(s - 0.5 * h)
                + np.log(
                    np.sum(
                        np.exp((s - 0.5 * h) - np.max(s - 0.5 * h))
                        * i0e(s)
                        * weights
                    )
                )
            )
            approx = float(lut_phase(sh, hh, bounds_check=False))
            self.assertAlmostEqual(approx, exact, places=1)

    def test_single_template_snr_snippets(self):
        """Verify SingleTemplate + DistMarg evaluates marginal logLR directly
        on pre-reconstructed dual-detector SNR snippets around the peak.
        """
        dt = 1.0 / 4096.0
        half_P = 8
        half_K = 46
        P = 2 * half_P + 1
        K = 2 * half_K + 1

        z_p = np.zeros(P, dtype=np.complex64)
        z_p[half_P - 1 : half_P + 2] = [
            6.0 * np.exp(0.5j),
            8.5 * np.exp(0.7j),
            6.2 * np.exp(0.9j),
        ]
        z_s_sig = np.zeros(K, dtype=np.complex64)
        z_s_sig[half_K + 4 : half_K + 7] = [
            -4.5 * np.exp(0.5j),
            -6.5 * np.exp(0.7j),
            -4.6 * np.exp(0.9j),
        ]
        z_s_noise = np.zeros(K, dtype=np.complex64)
        z_s_noise[:] = 0.2

        sig_p = 1000.0
        sig_s = 1000.0
        sh_sig = {
            "H1": TimeSeries(z_p * sig_p, delta_t=dt, epoch=-half_P * dt),
            "L1": TimeSeries(z_s_sig * sig_s, delta_t=dt, epoch=-half_K * dt),
        }
        hh = {"H1": sig_p ** 2, "L1": sig_s ** 2}

        vparams = ["tc", "ra", "dec", "inclination", "polarization", "distance"]
        prior = JointDistribution(
            vparams,
            Uniform(tc=(-6 * dt, 6 * dt)),
            UniformAngle(ra=(0.0, 2.0)),
            CosAngle(dec=None),
            SinAngle(inclination=None),
            UniformAngle(polarization=(0.0, 2.0)),
            UniformRadius(distance=(1.0, 10000.0)),
        )

        model = SingleTemplate(
            vparams,
            data=None,
            low_frequency_cutoff=None,
            sample_rate=4096,
            sh=sh_sig,
            hh=hh,
            prior=prior,
            static_params={},
            marginalize_phase=True,
            marginalize_distance=True,
            marginalize_distance_interpolator=True,
            marginalize_distance_snr_range=(0.02, 500.0),
            marginalize_distance_density=(128, 128),
            marginalize_vector_params="tc,ra,dec,inclination,polarization",
            marginalize_vector_samples=128,
            marginalize_sky_initial_samples=10000,
        )
        model.update()
        model.keep_ifos = ["H1", "L1"]
        model.tstart = {"H1": -half_P * dt, "L1": -half_K * dt}
        model.tend = {"H1": half_P * dt, "L1": half_K * dt}
        lr_coinc = float(model._loglr())

        model.sh["L1"] = TimeSeries(z_s_noise * sig_s, delta_t=dt, epoch=-half_K * dt)
        model.snr["L1"] = TimeSeries(z_s_noise, delta_t=dt, epoch=-half_K * dt)
        model.reset_vector_params()
        lr_noise = float(model._loglr())

        self.assertGreater(lr_coinc, lr_noise + 8.0)

    def test_reconstruct_fir_snr_span_td_bank(self):
        """Verify TimeDomainFilterBank-backed reconstruct_fir_snr_span matches
        direct correlation to float32 precision.
        """
        rng = np.random.default_rng(1234)
        ref_series = rng.standard_normal(1048576, dtype=np.float32)
        comp_tap = (
            rng.standard_normal(433) + 1j * rng.standard_normal(433)
        ).astype(np.complex64)
        shift = -120

        z_span = reconstruct_fir_snr_span(ref_series, comp_tap, shift, 100000, 105000)
        z_short = reconstruct_fir_snr_span(ref_series, comp_tap, shift, 100500, 100510)
        z_pt = reconstruct_fir_snr_span(ref_series, comp_tap, shift, 100503, 100504)
        np.testing.assert_allclose(z_span[500:510], z_short, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(z_span[503:504], z_pt, rtol=1e-5, atol=1e-5)

    def test_multi_detector_coinc_hdf_output(self):
        """Verify write_coinc_hdf and evaluate_coinc_windows handle 1, 2, and 3
        active detectors in a combined HDF5 file with non-participating detectors
        marked by time = -1.0 and trigger_id = -1.
        """
        import importlib.machinery
        import importlib.util
        import os
        import tempfile
        from pycbc.io import HFile

        script_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "bin",
            "all_sky_search",
            "pycbc_followup_marg_loglr",
        )
        loader = importlib.machinery.SourceFileLoader("followup_mod", script_path)
        spec = importlib.util.spec_from_loader("followup_mod", loader)
        mod = importlib.util.module_from_spec(spec)
        loader.exec_module(mod)

        instruments = ["H1", "L1", "V1"]
        sr = 2048.0
        half_P = 8
        half_K = 66
        P = 2 * half_P + 1
        K = 2 * half_K + 1

        model = mod.build_single_template_evaluator(
            instruments,
            sample_rate=sr,
            half_P=half_P,
            half_K=half_K,
            vsamples=64,
            sky_initial_samples=5000,
        )

        z_p = np.zeros(P, dtype=np.complex64)
        z_p[half_P] = 8.0 * np.exp(0.3j)
        z_s1 = np.zeros((2, K), dtype=np.complex64)
        z_s1[:, half_K] = 5.0 * np.exp(0.3j)
        z_s2 = np.zeros((2, K), dtype=np.complex64)
        z_s2[:, half_K] = 4.0 * np.exp(0.3j)

        lr_1d = mod.evaluate_coinc_windows(model, "H1", z_p, 1e6, sr)
        lr_2d = mod.evaluate_coinc_windows(
            model,
            "H1",
            z_p,
            1e6,
            sr,
            sec_snips_by_ifo={"L1": z_s1},
            sigmasq_by_ifo={"L1": 1e6},
        )
        lr_3d = mod.evaluate_coinc_windows(
            model,
            "H1",
            z_p,
            1e6,
            sr,
            sec_snips_by_ifo={"L1": z_s1, "V1": z_s2},
            sigmasq_by_ifo={"L1": 1e6, "V1": 1e6},
        )
        self.assertTrue(np.isfinite(lr_1d))
        self.assertEqual(lr_2d.shape, (2,))
        self.assertEqual(lr_3d.shape, (2,))

        with tempfile.NamedTemporaryFile(suffix=".hdf", delete=False) as tf:
            out_fn = tf.name
        try:
            segs_by_ifo = {
                "H1": ([100.0], [200.0]),
                "L1": ([100.0], [180.0]),
                "V1": ([120.0], [200.0]),
            }
            mod.write_coinc_hdf(
                out_fn,
                [np.array([lr_1d, lr_2d[0], lr_3d[0]], dtype=np.float32)],
                [np.array([0, 0, 0], dtype=np.int32)],
                [np.array([10, 11, 12], dtype=np.int32)],
                [np.array([1, 1, 1], dtype=np.int32)],
                {
                    "H1": [np.array([190.0, 110.0, 150.0])],
                    "L1": [np.array([-1.0, 110.002, 150.002])],
                    "V1": [np.array([-1.0, -1.0, 150.005])],
                },
                {
                    "H1": [np.array([1, 2, 3], dtype=np.int32)],
                    "L1": [np.array([-1, -1, -1], dtype=np.int32)],
                    "V1": [np.array([-1, -1, -1], dtype=np.int32)],
                },
                segs_by_ifo,
                instruments,
            )
            with HFile(out_fn, "r") as f:
                combos = [c.decode("utf8") for c in f["ifo_combination"][:]]
                self.assertEqual(combos, ["H1", "H1L1", "H1L1V1"])
                for ck in ("H1", "L1", "V1", "H1L1", "H1V1", "L1V1", "H1L1V1"):
                    self.assertIn(f"segments/{ck}/start", f)
                    self.assertIn(f"segments/{ck}/end", f)
        finally:
            if os.path.exists(out_fn):
                os.remove(out_fn)

    def test_tiered_decimation_and_snr_capped_coherent_stat(self):
        """Verify parse_loudest_keep_values, compute_sngl_lognoiserate, and
        SNR-capped normalized SingleTemplate coherent statistic behavior.
        """
        import importlib.machinery
        import importlib.util
        import os

        script_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "bin",
            "all_sky_search",
            "pycbc_followup_marg_loglr",
        )
        loader = importlib.machinery.SourceFileLoader("followup_mod2", script_path)
        spec = importlib.util.spec_from_loader("followup_mod2", loader)
        mod = importlib.util.module_from_spec(spec)
        loader.exec_module(mod)

        tiers = mod.parse_loudest_keep_values("[7.5:1,6.0:10,4.0:50,2.0:250,0.0:1000]")
        self.assertIsNotNone(tiers)
        stat_vals, total_factors = tiers
        np.testing.assert_allclose(stat_vals, [np.inf, 7.5, 6.0, 4.0, 2.0, 0.0])
        np.testing.assert_array_equal(total_factors, [1, 10, 50, 250, 1000])

        class DummyRankMethod:
            alphabelow = 5.0
            cap_alpha_below = True
            fits_by_tid = {
                "H1": {
                    "smoothed_fit_coeff": np.array([5.5], dtype=np.float32),
                    "smoothed_rate_above_thresh": np.array([1e-3], dtype=np.float32),
                    "thresh": 6.0,
                }
            }

        rm = DummyRankMethod()
        ln_r_hi = float(mod.compute_sngl_lognoiserate(rm, "H1", 0, 8.0))
        ln_r_lo = float(mod.compute_sngl_lognoiserate(rm, "H1", 0, 4.0))
        self.assertLess(ln_r_hi, ln_r_lo)

        model = mod.build_single_template_evaluator(
            ["H1", "L1"],
            sample_rate=2048.0,
            half_P=8,
            half_K=46,
            vsamples=64,
            sky_initial_samples=5000,
            d_min=0.001,
            d_max=10.0,
        )
        z_p = np.zeros(17, dtype=np.complex64)
        z_p[8] = 20.0 * np.exp(0.4j)
        z_s_coh = np.zeros((1, 93), dtype=np.complex64)
        z_s_coh[0, 46] = -15.0 * np.exp(0.4j)
        z_s_inc = np.zeros((1, 93), dtype=np.complex64)
        z_s_inc[0, 46] = 15.0 * np.exp(0.4j + 1.57j)

        gamma = 9.0 / np.hypot(20.0, 15.0)
        lr_coh = float(
            mod.evaluate_coinc_windows(
                model, "H1", z_p * gamma, 1.0, 2048.0, "L1", z_s_coh * gamma, 1.0
            )[0]
        )
        lr_inc = float(
            mod.evaluate_coinc_windows(
                model, "H1", z_p * gamma, 1.0, 2048.0, "L1", z_s_inc * gamma, 1.0
            )[0]
        )
        self.assertGreater(lr_coh, lr_inc + 5.0)


if __name__ == "__main__":
    unittest.main()
