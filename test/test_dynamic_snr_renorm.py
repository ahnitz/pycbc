"""The dynamic SNR renormalisation envelope, evaluated on demand.

DynamicSNRRenormFactor must give, at any requested indices, exactly what the
full envelope gives there; the envelope itself is checked against a direct
per-sample evaluation of its definition (hollow-window mean of 0.5|z|^2).
"""
import numpy as np
import pytest

from pycbc.filter.matched_filter_ratio import DynamicSNRRenormFactor, get_dynamic_snr_renorm_factor


def _direct(z, dt, window, hollow, floor=1.0, boost=None, scale=None):
    """The envelope's definition, one sample at a time (float64)."""
    p = 0.5 * np.abs(z.astype(np.complex128)) ** 2
    if scale is not None:
        p = p / scale ** 2
    n = len(z)
    w_out = int(round(0.5 * window / dt))
    w_in = int(round(hollow / dt))
    out = np.empty(n)
    for t in range(n):
        lo, hi = max(t - w_out, 0), min(t + w_out + 1, n)
        ilo, ihi = max(t - w_in, 0), min(t + w_in + 1, n)
        s = p[lo:hi].sum() - p[ilo:ihi].sum()
        c = (hi - lo) - (ihi - ilo)
        var = max(s / max(c, 1), max(floor, 1e-12))
        out[t] = 1.0 / np.sqrt(var)
    if boost is not None:
        out = np.minimum(out, boost)
    return out


@pytest.mark.parametrize("kw", [dict(), dict(variance_floor=0.25, max_boost_factor=1.3), dict(scale=3.0)])
def test_lazy_envelope_matches_definition_and_full_evaluation(kw):
    rng = np.random.default_rng(1)
    n, dt = 4000, 1.0 / 256
    z = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) * kw.get("scale", 1.0)).astype(np.complex64)
    z[1500:1520] *= 8.0                       # a loud excursion the hollow window must exclude
    lazy = DynamicSNRRenormFactor(z, dt=dt, window_duration=2.0, hollow_duration=0.25, **kw)
    full = get_dynamic_snr_renorm_factor(z, dt=dt, window_duration=2.0, hollow_duration=0.25, **kw)
    assert len(lazy) == n and full.dtype == np.float32
    idx = rng.integers(0, n, 300)
    np.testing.assert_array_equal(lazy[idx], full[idx])
    ref = _direct(z, dt, 2.0, 0.25, floor=kw.get("variance_floor", 1.0),
                  boost=kw.get("max_boost_factor"), scale=kw.get("scale"))
    np.testing.assert_allclose(full, ref, rtol=2e-5)


def test_disabled_and_degenerate_inputs_give_unit_factors(monkeypatch):
    z = np.ones(10, np.complex64)
    assert np.array_equal(DynamicSNRRenormFactor(z, dt=0.0)[np.arange(10)], np.ones(10, np.float32))
    monkeypatch.setenv("PYCBC_DISABLE_DYNAMIC_SNR_RENORM", "1")
    assert np.array_equal(DynamicSNRRenormFactor(z, dt=1.0)[[0, 5]], np.ones(2, np.float32))
