# Copyright (C) 2026 Alex Nitz, PyCBC Team
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

"""Unit and integration tests for GWOSC HDF5 strain reading support in PyCBC.
"""

import argparse
import glob
import os
import shutil
import tempfile
import unittest

import h5py
import numpy as np

from pycbc.frame import read_frame, locations_to_cache
from pycbc.frame.gwosc_hdf import (
    is_gwosc_hdf_file,
    parse_gwosc_hdf_filename,
    get_gwosc_hdf_metadata,
    resolve_gwosc_hdf_channel,
    read_frame_gwosc_hdf,
    extract_files_from_cache,
)
from pycbc.strain import (
    insert_strain_option_group,
    insert_strain_option_group_multi_ifo,
    from_cli,
    from_cli_multi_ifos,
)
from pycbc.types import TimeSeries


class TestGWOSCHDFSynthetic(unittest.TestCase):
    """Unit tests using synthetically generated GWOSC-format HDF5 files."""

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp(prefix='pycbc_test_gwosc_hdf_')
        cls.rate = 2048.0
        cls.dt = 1.0 / cls.rate
        cls.file_dur = 100.0  # 100s per file

        # File 1: [1000.0, 1100.0)
        cls.f1_start = 1000.0
        cls.f1_path = os.path.join(
            cls.temp_dir, f"H-H1_GWOSC_TEST_4KHZ-{int(cls.f1_start)}-{int(cls.file_dur)}.hdf5"
        )
        cls._create_mock_gwosc_hdf(cls.f1_path, 'H1', cls.f1_start, cls.file_dur, cls.rate)

        # File 2: [1100.0, 1200.0) - contiguous with File 1
        cls.f2_start = 1100.0
        cls.f2_path = os.path.join(
            cls.temp_dir, f"H-H1_GWOSC_TEST_4KHZ-{int(cls.f2_start)}-{int(cls.file_dur)}.hdf5"
        )
        cls._create_mock_gwosc_hdf(cls.f2_path, 'H1', cls.f2_start, cls.file_dur, cls.rate)

        # File 3: [1250.0, 1350.0) - creates a gap of 50s after File 2
        cls.f3_start = 1250.0
        cls.f3_path = os.path.join(
            cls.temp_dir, f"H-H1_GWOSC_TEST_4KHZ-{int(cls.f3_start)}-{int(cls.file_dur)}.hdf5"
        )
        cls._create_mock_gwosc_hdf(cls.f3_path, 'H1', cls.f3_start, cls.file_dur, cls.rate)

        # Non-standard filename file: [2000.0, 2050.0)
        cls.f_nonstd_start = 2000.0
        cls.f_nonstd_dur = 50.0
        cls.f_nonstd_path = os.path.join(cls.temp_dir, "custom_named_file.hdf5")
        cls._create_mock_gwosc_hdf(
            cls.f_nonstd_path, 'L1', cls.f_nonstd_start, cls.f_nonstd_dur, cls.rate
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)

    @classmethod
    def _create_mock_gwosc_hdf(cls, path, ifo, gps_start, duration, rate):
        n_pts = int(round(duration * rate))
        t = gps_start + np.arange(n_pts) / rate
        # Distinct continuous signal: sin(2 * pi * 5 * t)
        strain_data = np.sin(2.0 * np.pi * 5.0 * t).astype(np.float32)

        with h5py.File(path, 'w') as h5f:
            meta = h5f.create_group('meta')
            meta.create_dataset('GPSstart', data=int(gps_start))
            meta.create_dataset('Duration', data=int(duration))
            meta.create_dataset('Detector', data=ifo.encode('utf-8'))
            meta.create_dataset('Observatory', data=ifo[0].encode('utf-8'))
            meta.create_dataset('StrainChannel', data=f'{ifo}:GWOSC-TEST_STRAIN'.encode('utf-8'))
            meta.create_dataset('FrameType', data=f'{ifo}_GWOSC_TEST'.encode('utf-8'))

            strain_grp = h5f.create_group('strain')
            dset = strain_grp.create_dataset('Strain', data=strain_data)
            dset.attrs['Xstart'] = gps_start
            dset.attrs['Xspacing'] = 1.0 / rate
            dset.attrs['Npoints'] = n_pts
            dset.attrs['Xunits'] = 'second'
            dset.attrs['Yunits'] = ''

            # Mock quality mask
            quality = h5f.create_group('quality/simple')
            quality.create_dataset('DQmask', data=np.ones(int(duration), dtype=np.uint32))

    def test_format_detection(self):
        self.assertTrue(is_gwosc_hdf_file(self.f1_path))
        self.assertTrue(is_gwosc_hdf_file(self.f_nonstd_path))
        self.assertFalse(is_gwosc_hdf_file("some_frame.gwf"))
        self.assertFalse(is_gwosc_hdf_file("some_file.xml"))

    def test_filename_parsing(self):
        meta = parse_gwosc_hdf_filename(self.f1_path)
        self.assertIsNotNone(meta)
        obs, ifo, desc, start, dur = meta
        self.assertEqual(obs, 'H')
        self.assertEqual(ifo, 'H1')
        self.assertEqual(start, 1000.0)
        self.assertEqual(dur, 100.0)

        # Non-standard name should return None from parse_gwosc_hdf_filename
        self.assertIsNone(parse_gwosc_hdf_filename(self.f_nonstd_path))

    def test_metadata_extraction_fallback(self):
        # Header inspection for non-standard filename
        meta = get_gwosc_hdf_metadata(self.f_nonstd_path)
        self.assertEqual(meta['ifo'], 'L1')
        self.assertEqual(meta['start_time'], 2000.0)
        self.assertEqual(meta['duration'], 50.0)
        self.assertEqual(meta['end_time'], 2050.0)

    def test_channel_resolution(self):
        with h5py.File(self.f1_path, 'r') as h5f:
            self.assertEqual(resolve_gwosc_hdf_channel(h5f, 'H1:STRAIN'), 'strain/Strain')
            self.assertEqual(resolve_gwosc_hdf_channel(h5f, 'strain/Strain'), 'strain/Strain')
            self.assertEqual(resolve_gwosc_hdf_channel(h5f, 'H1:GWOSC-TEST_STRAIN'), 'strain/Strain')
            self.assertEqual(resolve_gwosc_hdf_channel(h5f, 'quality/simple/DQmask'), 'quality/simple/DQmask')

            # Mismatched detector should raise
            with self.assertRaises(ValueError):
                resolve_gwosc_hdf_channel(h5f, 'L1:STRAIN')

    def test_single_file_subinterval_slice(self):
        start = 1010.0
        end = 1030.0
        ts = read_frame_gwosc_hdf(self.f1_path, 'H1:STRAIN', start_time=start, end_time=end)
        self.assertIsInstance(ts, TimeSeries)
        self.assertEqual(len(ts), int((end - start) * self.rate))
        self.assertAlmostEqual(float(ts.start_time), start)
        self.assertAlmostEqual(float(ts.end_time), end)

        # Verify waveform values match sin(2 * pi * 5 * t)
        expected_t = start + np.arange(len(ts)) * self.dt
        expected_val = np.sin(2.0 * np.pi * 5.0 * expected_t).astype(np.float32)
        np.testing.assert_allclose(ts.numpy(), expected_val, atol=1e-5)

    def test_multi_file_boundary_crossing(self):
        # Spans across f1 (1000-1100) and f2 (1100-1200)
        start = 1080.0
        end = 1120.0
        ts = read_frame_gwosc_hdf(
            [self.f1_path, self.f2_path], 'H1:STRAIN', start_time=start, end_time=end
        )
        self.assertEqual(len(ts), int((end - start) * self.rate))
        self.assertAlmostEqual(float(ts.start_time), start)
        self.assertAlmostEqual(float(ts.end_time), end)

        # Verify signal continuity across boundary at t = 1100
        expected_t = start + np.arange(len(ts)) * self.dt
        expected_val = np.sin(2.0 * np.pi * 5.0 * expected_t).astype(np.float32)
        np.testing.assert_allclose(ts.numpy(), expected_val, atol=1e-5)

    def test_data_gap_detection(self):
        # Requested interval [1180, 1270) spans across gap [1200, 1250) between f2 and f3
        with self.assertRaises(ValueError) as ctx:
            read_frame_gwosc_hdf(
                [self.f2_path, self.f3_path], 'H1:STRAIN', start_time=1180.0, end_time=1270.0
            )
        self.assertIn("gap", str(ctx.exception).lower())

    def test_out_of_bounds_detection(self):
        # Requesting data before available span
        with self.assertRaises(ValueError):
            read_frame_gwosc_hdf(self.f1_path, 'H1:STRAIN', start_time=900.0, end_time=1050.0)

        # Requesting data after available span
        with self.assertRaises(ValueError):
            read_frame_gwosc_hdf(self.f1_path, 'H1:STRAIN', start_time=1050.0, end_time=1150.0)

    def test_read_via_cache_file(self):
        cache_path = os.path.join(self.temp_dir, "test_cache.cache")
        with open(cache_path, 'w') as cf:
            cf.write(f"H H1_GWOSC 1000 100 file://localhost{self.f1_path}\n")
            cf.write(f"H H1_GWOSC 1100 100 file://localhost{self.f2_path}\n")

        ts = read_frame(cache_path, 'H1:STRAIN', start_time=1050.0, end_time=1150.0)
        self.assertEqual(len(ts), int(100.0 * self.rate))
        self.assertAlmostEqual(float(ts.start_time), 1050.0)

    def test_read_via_glob(self):
        glob_pattern = os.path.join(self.temp_dir, "H-H1_GWOSC_TEST_4KHZ-1*.hdf5")
        ts = read_frame(glob_pattern, 'H1:STRAIN', start_time=1050.0, end_time=1150.0)
        self.assertEqual(len(ts), int(100.0 * self.rate))

    def test_locations_to_cache(self):
        cache = locations_to_cache([self.f1_path, self.f2_path])
        self.assertIsNotNone(cache)
        self.assertEqual(cache.length, 2)


class TestGWOSCHDFRealData(unittest.TestCase):
    """Integration tests running against real GWOSC data in /home/ahnitz/gwdata."""

    def test_o3b_data_reading(self):
        o3b_dir = '/home/ahnitz/gwdata/O3b/H1/1256194048'
        if not os.path.exists(o3b_dir):
            self.skipTest(f"O3b data not found at {o3b_dir}")

        files = sorted(glob.glob(f"{o3b_dir}/*.hdf5"))[:2]
        self.assertGreaterEqual(len(files), 2)

        # Single file read
        ts = read_frame(files[0], 'H1:STRAIN', start_time=1256661000, end_time=1256661010)
        self.assertEqual(len(ts), 10 * 2048)
        self.assertEqual(ts.delta_t, 1.0 / 2048.0)

        # Boundary crossing between consecutive files (boundary at 1256665088)
        ts_cross = read_frame(files, 'H1:STRAIN', start_time=1256665000, end_time=1256665200)
        self.assertEqual(len(ts_cross), 200 * 2048)
        self.assertAlmostEqual(float(ts_cross.start_time), 1256665000.0)

    def test_o4a_multi_detector_reading(self):
        h1_file = '/home/ahnitz/gwdata/O4a/H1/1367343104/H-H1_DATASAMPLER_O4a_2048HZ-1368195072-4096.hdf5'
        l1_file = '/home/ahnitz/gwdata/O4a/L1/1367343104/L-L1_DATASAMPLER_O4a_2048HZ-1368195072-4096.hdf5'
        if not (os.path.exists(h1_file) and os.path.exists(l1_file)):
            self.skipTest("O4a real data files not found")

        parser = argparse.ArgumentParser()
        insert_strain_option_group_multi_ifo(parser)
        args = parser.parse_args([
            '--frame-files', f'H1:{h1_file}', f'L1:{l1_file}',
            '--channel-name', 'H1:H1:STRAIN', 'L1:L1:STRAIN',
            '--gps-start-time', '1368195100',
            '--gps-end-time', '1368195120',
            '--sample-rate', 'H1:2048', 'L1:2048'
        ])

        strains = from_cli_multi_ifos(args, ['H1', 'L1'])
        self.assertIn('H1', strains)
        self.assertIn('L1', strains)
        self.assertEqual(len(strains['H1']), 20 * 2048)
        self.assertEqual(len(strains['L1']), 20 * 2048)

    def test_o2_data_reading(self):
        o2_file = '/home/ahnitz/gwdata/O2/H1/1163919360/H-H1_DATASAMPLER_O2_2048HZ-1164558336-4096.hdf5'
        if not os.path.exists(o2_file):
            self.skipTest("O2 real data file not found")

        parser = argparse.ArgumentParser()
        insert_strain_option_group(parser)
        args = parser.parse_args([
            '--frame-files', o2_file,
            '--channel-name', 'H1:STRAIN',
            '--gps-start-time', '1164559000',
            '--gps-end-time', '1164559030',
            '--sample-rate', '2048'
        ])

        strain = from_cli(args)
        self.assertEqual(len(strain), 30 * 2048)
        self.assertAlmostEqual(float(strain.start_time), 1164559000.0)


if __name__ == '__main__':
    unittest.main()
