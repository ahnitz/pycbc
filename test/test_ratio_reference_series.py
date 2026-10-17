import unittest
import numpy as np

import pycbc.filter
from pycbc.fft import backend_cpu
from pycbc.filter.matchedfilter import matched_filter_core
from pycbc.types import FrequencySeries
from pycbc.filter import matched_filter_ratio as mfr


def _spectra(seed, n_freq, delta_f):
    rng = np.random.default_rng(seed)
    data = rng.standard_normal(n_freq) + 1j * rng.standard_normal(n_freq)
    return FrequencySeries(data.astype(np.complex64), delta_f=delta_f)


@unittest.skipIf(mfr._mf is None, "matchedfilter is not installed")
class TestRealReferenceSeries(unittest.TestCase):
    """The top filter's real output series.

    It must equal the real part of pycbc's own complex SNR series, scaled the
    same way, and must not depend on earlier calls: the complex-to-real inverse
    transform may overwrite its input (FFTW does by default), so the cached
    input buffer has to be fully rewritten on every call. Every available CPU
    FFT backend is checked, each on a second call with new data.
    """

    f_lower = 20.0
    f_high = 900.0

    def setUp(self):
        self.delta_f = 1.0 / 16
        self.n = 2048 * 16
        self.n_freq = self.n // 2 + 1
        rng = np.random.default_rng(3)
        psd = np.abs(rng.standard_normal(self.n_freq)) + 1.0
        self.psd = FrequencySeries(psd.astype(np.float32), delta_f=self.delta_f)
        self.template = _spectra(1, self.n_freq, self.delta_f)
        self.template.f_lower = self.f_lower
        # bank templates carry sigmasq(psd); give the test template the same
        self.template.sigmasq = lambda psd: float(pycbc.filter.sigmasq(
            self.template, psd, low_frequency_cutoff=self.f_lower,
            high_frequency_cutoff=self.f_high))

    def _expected(self, stilde):
        """Re of pycbc's complex SNR series for overwhitened data (no PSD division)."""
        h_norm = self.template.sigmasq(self.psd)
        snr, _, norm = matched_filter_core(
            self.template, stilde, psd=None, low_frequency_cutoff=self.f_lower,
            high_frequency_cutoff=self.f_high, h_norm=h_norm)
        return (snr.numpy() * (norm * stilde.delta_t)).real

    def test_second_call_matches_pycbc_snr(self):
        default = backend_cpu.cpu_backend
        first = _spectra(10, self.n_freq, self.delta_f)
        second = _spectra(11, self.n_freq, self.delta_f)
        for backend in backend_cpu._alist:
            with self.subTest(backend=backend):
                backend_cpu.set_backend([backend])
                try:
                    engine = mfr.MatchedFilterRatioControl(
                        6.0, self.delta_f, high_frequency_cutoff=self.f_high)
                    engine.make_reference_series(first, self.psd, self.template)
                    got, h_norm = engine.make_reference_series(second, self.psd, self.template)
                    expect = self._expected(second)
                finally:
                    backend_cpu.set_backend([default])
                self.assertEqual(got.dtype, np.float32)
                self.assertEqual(h_norm, self.template.sigmasq(self.psd))
                scale = np.abs(expect).max()
                np.testing.assert_allclose(got, expect, rtol=0, atol=1e-5 * scale)


if __name__ == '__main__':
    unittest.main()
