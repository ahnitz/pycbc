import unittest
import numpy as np
import pycbc.types
import pycbc.waveform
import pycbc.filter

class TestRangeAndPSDCutoff(unittest.TestCase):
    def test_sigma_high_frequency_cutoff(self):
        # Create a frequency series with low values near Nyquist (stopband simulation)
        df = 0.5
        n_pts = 1025 # 0 to 512 Hz
        # Realistic scaled PSD (with dyn_range_fac applied)
        psd_data = np.ones(n_pts, dtype=np.float32) * 1e8
        # Stopband rolloff near Nyquist (strain filtered down by 80 dB -> power down by 160 dB)
        psd_data[int(460/df):] = 1e-4 # steep drop in PSD
        psd = pycbc.types.FrequencySeries(psd_data, delta_f=df)

        flow = 30.0
        delta_t = 1.0 / (2 * 512)
        out = pycbc.types.zeros(n_pts, dtype=np.complex64)
        htilde = pycbc.waveform.get_waveform_filter(
            out, mass1=1.4, mass2=1.4, approximant="TaylorF2",
            f_lower=flow, delta_f=df, delta_t=delta_t,
            distance=1.0 / pycbc.DYN_RANGE_FAC
        ).astype(np.complex64)

        # Without high_frequency_cutoff, integration enters stopband where 1/PSD blows up
        sig_unbounded = pycbc.filter.sigma(htilde, psd=psd, low_frequency_cutoff=flow)
        
        # With high_frequency_cutoff bounded away from stopband
        f_high = 450.0
        sig_bounded = pycbc.filter.sigma(htilde, psd=psd, low_frequency_cutoff=flow, high_frequency_cutoff=f_high)

        self.assertGreater(sig_unbounded, 10.0 * sig_bounded)
        self.assertTrue(np.isfinite(sig_bounded))
        self.assertGreater(sig_bounded, 0)

if __name__ == '__main__':
    unittest.main()
