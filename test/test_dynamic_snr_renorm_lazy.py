"""The dynamic SNR renormalisation envelope, evaluated on demand.

DynamicSNRRenormFactor must give, at any requested indices, exactly what the
full envelope gives there; the envelope itself is checked against a direct
per-sample evaluation of its definition (hollow-window mean of 0.5|z|^2).
"""
import numpy as np
import pytest

from pycbc.filter.matched_filter_ratio import DynamicSNRRenormFactor, get_dynamic_snr_renorm_factor


def _direct(z, dt, window, hollow, floor=1.0, boost=None, scale=None):
    """The envelope's definition, one sample at a time (float64): the hollow-window
    mean of 0.5|z|^2, or of x^2 for a real series x, with non-finite samples' power
    taken as zero."""
    if np.iscomplexobj(z):
        p = 0.5 * np.abs(z.astype(np.complex128)) ** 2
    else:
        p = z.astype(np.float64) ** 2
    p = np.where(np.isfinite(p), p, 0.0)
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


def test_degenerate_inputs_give_unit_factors():
    z = np.ones(10, np.complex64)
    assert np.array_equal(DynamicSNRRenormFactor(z, dt=0.0)[np.arange(10)], np.ones(10, np.float32))
    assert np.array_equal(DynamicSNRRenormFactor(z[:0], dt=1.0)[[]], np.ones(0, np.float32))


@pytest.mark.parametrize("n", [4000, 4003])
def test_real_series_matches_definition(n):
    """A real series (the top filter's output in a real-real-complex search), with a
    length that is and is not a whole number of the envelope's blocks."""
    rng = np.random.default_rng(2)
    dt = 1.0 / 256
    x = (rng.standard_normal(n) * 3.0).astype(np.float32)
    x[1500:1520] *= 8.0
    full = get_dynamic_snr_renorm_factor(x, dt=dt, window_duration=2.0, hollow_duration=0.25, scale=3.0)
    assert full.dtype == np.float32
    ref = _direct(x, dt, 2.0, 0.25, scale=3.0)
    np.testing.assert_allclose(full, ref, rtol=2e-5)


def test_non_finite_samples_count_as_zero_power():
    rng = np.random.default_rng(3)
    n, dt = 4000, 1.0 / 256
    x = rng.standard_normal(n).astype(np.float32)
    x[[100, 2000, 3999]] = [np.nan, np.inf, -np.inf]
    full = get_dynamic_snr_renorm_factor(x, dt=dt, window_duration=2.0, hollow_duration=0.25)
    np.testing.assert_allclose(full, _direct(x, dt, 2.0, 0.25), rtol=2e-5)
