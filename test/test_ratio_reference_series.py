import unittest
import numpy as np

import pycbc.filter
from pycbc.fft import backend_cpu
from pycbc.types import FrequencySeries
from pycbc.filter import matched_filter_ratio as mfr


def _spectra(seed, n_freq, delta_f):
    rng = np.random.default_rng(seed)
    data = rng.standard_normal(n_freq) + 1j * rng.standard_normal(n_freq)
    return FrequencySeries(data.astype(np.complex64), delta_f=delta_f)


@unittest.skipIf(mfr._mf is None, "matchedfilter is not installed")
class TestRealReferenceSeries(unittest.TestCase):
    """The real reference series must not depend on earlier calls.

    The complex-to-real inverse transform may overwrite its input (FFTW does
    by default), so a cached input buffer must be fully rewritten on every
    call. Each backend is checked: a second call on new data must equal the
    real part of the complex series for that data, computed by a fresh engine.
    """

    def setUp(self):
        self.delta_f = 1.0 / 16
        self.n = 2048 * 16
        self.n_freq = self.n // 2 + 1
        rng = np.random.default_rng(3)
        psd = np.abs(rng.standard_normal(self.n_freq)) + 1.0
        self.psd = FrequencySeries(psd.astype(np.float32), delta_f=self.delta_f)
        self.template = _spectra(1, self.n_freq, self.delta_f)
        self.template.f_lower = 20.0
        # bank templates carry sigmasq(psd); give the test template the same
        self.template.sigmasq = lambda psd: float(pycbc.filter.sigmasq(
            self.template, psd, low_frequency_cutoff=20.0, high_frequency_cutoff=900.0))

    def _engine(self):
        return mfr.MatchedFilterRatioControl(6.0, self.delta_f, high_frequency_cutoff=900.0)

    def test_second_call_matches_fresh_engine(self):
        default = backend_cpu.cpu_backend
        for backend in backend_cpu._alist:
            with self.subTest(backend=backend):
                backend_cpu.set_backend([backend])
                try:
                    engine = self._engine()
                    first = _spectra(10, self.n_freq, self.delta_f)
                    second = _spectra(11, self.n_freq, self.delta_f)
                    engine.make_reference_series(first, self.psd, self.template, real=True)
                    got, _ = engine.make_reference_series(second, self.psd, self.template, real=True)
                    expect, _ = self._engine().make_reference_series(
                        second, self.psd, self.template, real=False)
                finally:
                    backend_cpu.set_backend([default])
                self.assertEqual(got.dtype, np.float32)
                scale = np.abs(expect.real).max()
                np.testing.assert_allclose(got, expect.real, rtol=0, atol=1e-5 * scale)


if __name__ == '__main__':
    unittest.main()
