# Copyright (C) 2026 Alex Nitz
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
"""Unit tests for single-detector ranking functions and tolerances."""

import numpy as np

from pycbc.events import ranking, stat


def test_newsnr_and_sgveto_defaults_and_scalar_types():
    """Verify default tolerances match historical formulas and scalar types."""
    # Scalar inputs return scalar values without downweighting below thresholds
    assert np.isclose(ranking.newsnr(10.0, 1.0), 10.0)
    assert np.isscalar(ranking.newsnr(10.0, 1.0))
    assert np.isclose(ranking.newsnr_sgveto(10.0, 1.0, 4.0), 10.0)
    assert np.isscalar(ranking.newsnr_sgveto(10.0, 1.0, 4.0))

    # Above chisq=1 and sgchisq=4, match exact analytic penalty
    expected_newsnr = 10.0 * (0.5 * (1.0 + 2.0**3.0)) ** (-1.0 / 6.0)
    assert np.isclose(ranking.newsnr(10.0, 2.0), expected_newsnr)
    assert np.isclose(ranking.newsnr_sgveto(10.0, 2.0, 16.0), expected_newsnr / 2.0)


def test_snr2_scaled_chisq_and_sgchisq_tolerance():
    """Verify SNR^2-scaled tolerances preserve loud signals and reject glitches."""
    snr = np.array([12.0, 40.0])
    brchisq = np.array([1.2, 1.2])
    sgchisq = np.array([30.0, 25.0])

    # Without scaling, loud SNR=40 signal with sg_chisq=25 (HM residual) is penalized
    unscaled = ranking.newsnr_sgveto(snr, brchisq, sgchisq)
    assert unscaled[1] < 20.0

    # With sgchisq_scale=0.02 and chisq_thresh=1.25, loud signal stays unpenalized
    # while moderate SNR=12 glitch (sg_chisq=30 > 4 + 0.02*144 = 6.88) stays suppressed
    scaled = ranking.newsnr_sgveto(
        snr,
        brchisq,
        sgchisq,
        chisq_thresh=1.25,
        sgchisq_scale=0.02,
    )
    assert np.isclose(scaled[1], 40.0)
    assert scaled[0] < 6.0


def test_newsnr_sgveto_autochisq_and_stat_keyword_forwarding():
    """Verify auto-chisq reweighting, caps, fallback, and Stat kwarg forwarding."""
    # Scalar input preserves scalar output and chisq_scale/chisq_max_snr2 capping
    scalar_val = ranking.newsnr_sgveto_autochisq(10.0, 1.0, 2.0, 1.0)
    assert np.isscalar(scalar_val)
    assert np.isclose(scalar_val, 10.0)
    capped = ranking.newsnr(
        np.array([10.0, 30.0]),
        np.array([2.0, 2.0]),
        chisq_thresh=1.0,
        chisq_scale=0.005,
        chisq_max_snr2=200.0,
    )
    assert np.isclose(capped[0], ranking.newsnr(10.0, 2.0, chisq_thresh=1.5))
    assert np.isclose(capped[1], ranking.newsnr(30.0, 2.0, chisq_thresh=2.0))

    trigs = {
        "snr": np.array([12.0, 30.0, 15.0], dtype=np.float32),
        "chisq": np.array([16.0, 16.0, 16.0], dtype=np.float32),
        "chisq_dof": np.array([9, 9, 9], dtype=np.int32),
        "sg_chisq": np.array([2.0, 18.0, 2.0], dtype=np.float32),
        "cont_chisq": np.array([400.0, 40.0, 99.0], dtype=np.float32),
        "cont_chisq_dof": np.array([40, 40, 0], dtype=np.int32),
    }
    stat_obj = stat.Stat(
        "newsnr_sgveto_autochisq",
        sngl_ranking_sgchisq_scale="0.02",
        sngl_ranking_autochisq_thresh="1.2",
        sngl_ranking_autochisq_scale="0.001",
        sngl_ranking_autochisq_max_snr2="400.0",
    )
    vals = stat_obj.get_sngl_ranking(trigs)
    # Trigger 0 penalized by auto-chisq; 1 is clean; 2 uses dof==0 fallback
    assert vals[0] < 7.0
    assert np.isclose(vals[1], 30.0)
    assert np.isclose(vals[2], 15.0)
