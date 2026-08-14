# Copyright (C) 2026  Alex Nitz
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
"""Tests that the model's detectors are referenced at the time analyzed.

Detector estimates sidereal time by advancing it at a constant rate from a
reference time, so it is exact at that reference and drifts away from it.
The default reference is the time of GW150914, which is right for data from
2015 and progressively wrong for everything since. At the time of GW170817,
about two years later, the antenna pattern is off by up to 1.1e-4 on
factors of order one and the arrival time by up to 1.2 microseconds, both
taken over the whole sky; at the position used below it is 0.5 to 0.7
microseconds. The sky draws already reference the time being analyzed, so a
detector left at the default was both wrong and inconsistent with them.
"""

import unittest

import numpy

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.time import gmst_accurate
from pycbc.waveform import get_td_waveform
from pycbc.waveform.generator import (FDomainCBCGenerator,
                                      FDomainDetFrameGenerator,
                                      FDomainDetFrameTwoPolGenerator)

# well after the default reference, which is what makes the drift visible
TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 16, 2048
DEFAULT_REFERENCE = 1126259462.0


class TestRelbinDetectorReference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=60., inclination=0.5, coa_phase=1.1)
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        seed = 17
        for ifo in ['H1', 'L1']:
            noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                                   seed=seed)
            noise._epoch = TC - SEGLEN / 2
            seed += 53
            signal = Detector(ifo).project_wave(hp, hc, 1.7, -0.4, 0.3)
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.flow = {ifo: FLOW for ifo in cls.data}

    def model(self):
        variable = ['distance', 'inclination']
        return models.Relative(
            list(variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                               approximant='TaylorF2', tc=TC,
                               ra=1.7, dec=-0.4, polarization=0.3),
            prior=JointDistribution(list(variable),
                                    SinAngle(inclination=None),
                                    Uniform(distance=(10, 200))),
            fiducial_params={'mass1': 1.4}, epsilon=0.1)

    def test_detectors_are_referenced_at_the_analyzed_time(self):
        """Not left at the default, which is the time of GW150914."""
        model = self.model()
        for ifo, det in model.det.items():
            self.assertAlmostEqual(det.reference_time, TC, places=3,
                                   msg="%s referenced at %r" % (
                                       ifo, det.reference_time))
            self.assertNotEqual(det.reference_time, DEFAULT_REFERENCE)

    def test_sidereal_time_is_right_where_it_is_used(self):
        """The estimate must agree with the accurate one at the event.

        This is the property the reference time buys, and the reason the
        default is not good enough: the estimate is exact at its reference
        and drifts from it, so referencing it two years away leaves a real
        error at the time actually being analyzed.
        """
        model = self.model()
        exact = gmst_accurate(TC)
        for ifo, det in model.det.items():
            self.assertLess(abs(det.gmst_estimate(TC) - exact), 1e-9,
                            "%s sidereal time is off at the event" % ifo)

        # and the default really would have been wrong, so the test above
        # is not passing for free
        stale = Detector('H1').gmst_estimate(TC)
        self.assertGreater(abs(stale - exact), 1e-6)

    def test_antenna_response_matches_the_accurate_calculation(self):
        """What the drift costs, in the quantity the likelihood uses."""
        model = self.model()
        ra, dec, pol = 1.7, -0.4, 0.3
        for ifo, det in model.det.items():
            fp, fc = det.antenna_pattern(ra, dec, pol, TC)
            ref = Detector(ifo, reference_time=None)
            rfp, rfc = ref.antenna_pattern(ra, dec, pol, TC)
            scale = max(abs(rfp), abs(rfc))
            self.assertLess(max(abs(fp - rfp), abs(fc - rfc)) / scale, 1e-9,
                            "%s antenna response is off" % ifo)

    def test_marginalized_model_detectors_are_referenced_too(self):
        """The same drift was in the marginalized models.

        They build their detectors lazily, on the first evaluation, and
        cache them. The reference is therefore the coalescence time of
        whichever sample came first, which is right: the time varies only
        over the prior, a fraction of a second, while the drift being fixed
        here is measured over years.
        """
        variable = ['distance', 'inclination', 'tc']
        model = models.MarginalizedTime(
            list(variable), {k: v.copy() for k, v in self.data.items()},
            low_frequency_cutoff=self.flow, psds=self.psds,
            static_params=dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                               approximant='TaylorF2', ra=1.7, dec=-0.4,
                               polarization=0.3),
            prior=JointDistribution(list(variable),
                                    SinAngle(inclination=None),
                                    Uniform(distance=(10, 200)),
                                    Uniform(tc=(TC - 0.1, TC + 0.1))),
            marginalize_phase=True, marginalize_vector_params='tc',
            sample_rate=4096, marginalize_vector_samples=64)
        numpy.random.seed(3)
        model.update(distance=60., inclination=0.5)
        self.assertTrue(numpy.isfinite(model.loglr))

        self.assertTrue(model.dets, "no detectors were built")
        exact = gmst_accurate(TC)
        for det_name, det in model.dets.items():
            self.assertNotEqual(det.reference_time, DEFAULT_REFERENCE)
            self.assertLess(abs(det.gmst_estimate(TC) - exact), 1e-9,
                            "%s sidereal time is off at the event"
                            % det_name)


class TestGeneratorDetectorReference(unittest.TestCase):
    """The waveform generator builds detectors of its own, and used them
    at the default reference.

    Fixing the models was not enough. The detector-frame generators keep
    their own detectors and ask them for the arrival time that sets the
    time shift applied to the waveform, so a generator left at the default
    put the signal in the wrong place no matter what the model did with its
    own detectors. That shift is what these tests pin down: it has to agree
    with the accurate sidereal time at the time actually being analyzed.

    These are deliberately cheap. The generator is the whole subject here,
    so there is no data, no PSD and no model to build; a single waveform is
    generated per case at a coarse delta_f.
    """

    generator_classes = [FDomainDetFrameGenerator,
                         FDomainDetFrameTwoPolGenerator]

    def generator(self, cls):
        return cls(FDomainCBCGenerator, epoch=TC - SEGLEN / 2.,
                   detectors=['H1', 'L1'],
                   variable_args=['tc', 'ra', 'dec', 'polarization'],
                   delta_f=1. / SEGLEN, f_lower=FLOW,
                   approximant='TaylorF2', mass1=1.4, mass2=1.35)

    def test_detectors_are_referenced_at_the_analyzed_time(self):
        """Not left at the default, which is the time of GW150914.

        The reference cannot be set when the detectors are built, because
        tc is generally a variable parameter, so it is set on the first
        call to generate and kept from then on.
        """
        for cls in self.generator_classes:
            gen = self.generator(cls)
            gen.generate(tc=TC, ra=1.7, dec=-0.4, polarization=0.3)
            for ifo, det in gen.detectors.items():
                self.assertAlmostEqual(
                    det.reference_time, TC, places=3,
                    msg="%s %s referenced at %r" % (
                        cls.__name__, ifo, det.reference_time))
                self.assertNotEqual(det.reference_time, DEFAULT_REFERENCE)

    def test_arrival_time_matches_the_accurate_calculation(self):
        """The shift the generator applies must be the accurate one.

        ``generate`` shifts each detector's waveform by the arrival time its
        own detector reports, so that number is where the drift enters. Here
        it is checked against a detector that takes no shortcut and computes
        the sidereal time exactly at every call.
        """
        ra, dec = 1.7, -0.4
        for cls in self.generator_classes:
            gen = self.generator(cls)
            gen.generate(tc=TC, ra=ra, dec=dec, polarization=0.3)
            for ifo, det in gen.detectors.items():
                exact = Detector(ifo, reference_time=None).arrival_time(
                    TC, ra, dec, 'geocentric')
                self.assertLess(
                    abs(det.arrival_time(TC, ra, dec, 'geocentric') - exact),
                    1e-9, "%s %s arrival time is off" % (cls.__name__, ifo))

                # and the default really would have been wrong, so the check
                # above is not passing for free
                stale = Detector(ifo).arrival_time(TC, ra, dec, 'geocentric')
                self.assertGreater(abs(stale - exact), 1e-7)

    def test_a_marginalized_time_is_averaged(self):
        """tc arrives as a vector when it is being marginalized over.

        The reference has to be a single number, so it is the average, which
        is what the models and the sky draws already do. The spread is the
        width of the prior, far too small to matter to an estimate whose
        error is measured over years. This goes through the referencing
        directly rather than through ``generate``, because not every
        detector-frame generator accepts a vector time in the first place.
        """
        tc = TC + numpy.linspace(-0.1, 0.1, 16)
        for cls in self.generator_classes:
            gen = self.generator(cls)
            gen.reference_detectors(tc)
            for ifo, det in gen.detectors.items():
                self.assertAlmostEqual(
                    det.reference_time, tc.mean(), places=6,
                    msg="%s %s referenced at %r" % (
                        cls.__name__, ifo, det.reference_time))

    def test_the_reference_is_kept_once_it_is_set(self):
        """The detectors are built once and reused.

        The reference is therefore the time of whichever waveform was asked
        for first, and a later call does not move it. That is what makes the
        cost a single construction rather than one per likelihood
        evaluation, and it is accurate enough because the time only ever
        moves over its prior.
        """
        for cls in self.generator_classes:
            gen = self.generator(cls)
            gen.generate(tc=TC, ra=1.7, dec=-0.4, polarization=0.3)
            first = {ifo: det for ifo, det in gen.detectors.items()}
            gen.generate(tc=TC + 0.05, ra=1.2, dec=0.3, polarization=1.1)
            for ifo, det in gen.detectors.items():
                self.assertIs(det, first[ifo])
                self.assertAlmostEqual(det.reference_time, TC, places=3)


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestRelbinDetectorReference))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestGeneratorDetectorReference))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    from utils import simple_exit
    simple_exit(results)
