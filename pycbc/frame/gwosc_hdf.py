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

"""Module for reading gravitational wave strain from GWOSC HDF5 files.
"""

import glob
import logging
import math
import os
import re
import urllib.parse
import h5py
import numpy as np

from pycbc.types import TimeSeries

logger = logging.getLogger('pycbc.frame.gwosc_hdf')

HDF5_EXTENSIONS = ('.hdf5', '.h5', '.hdf')


def is_gwosc_hdf_file(file_path):
    """Check if the given path or URL points to an HDF5 file.

    Parameters
    ----------
    file_path : str
        File path or URL to check.

    Returns
    -------
    bool
        True if the file has an HDF5 extension or is recognized as HDF5 by h5py.
    """
    if not isinstance(file_path, str):
        return False

    path = urllib.parse.urlsplit(file_path).path
    ext = os.path.splitext(path)[1].lower()
    if ext in HDF5_EXTENSIONS:
        return True

    try:
        return bool(h5py.is_hdf5(path))
    except Exception:
        return False


def parse_gwosc_hdf_filename(file_path):
    """Parse standard LIGO/GWOSC filename convention to extract metadata.

    Standard naming conventions:
        [OBS]-[IFO]_[DESC]-[GPS]-[DUR].[EXT]
        [IFO]-[DESC]-[GPS]-[DUR].[EXT]

    Parameters
    ----------
    file_path : str
        Path or URL of the file.

    Returns
    -------
    tuple or None
        (obs, ifo, desc, start_time, duration) if successfully parsed, else None.
    """
    path = urllib.parse.urlsplit(file_path).path
    filename = os.path.basename(path)
    root, _ = os.path.splitext(filename)
    parts = root.split('-')

    if len(parts) >= 4:
        try:
            start_time = float(parts[-2])
            duration = float(parts[-1])
        except ValueError:
            return None

        part0 = parts[0]
        desc_part = parts[1]

        # Determine ifo and obs
        known_ifos = {'H1', 'L1', 'V1', 'K1', 'G1'}
        if part0 in known_ifos:
            ifo = part0
            obs = ifo[0]
            desc = desc_part
        elif '_' in desc_part and desc_part.split('_')[0] in known_ifos:
            ifo = desc_part.split('_')[0]
            obs = part0
            desc = desc_part
        else:
            ifo = part0
            obs = part0[0] if part0 else 'U'
            desc = desc_part

        return obs, ifo, desc, start_time, duration

    return None


def get_gwosc_hdf_metadata(file_path):
    """Get metadata for a GWOSC HDF5 file.

    Tries fast filename parsing first; falls back to inspecting HDF5 attributes.

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file.

    Returns
    -------
    dict
        Dictionary with keys 'obs', 'ifo', 'desc', 'start_time', 'duration',
        'sample_rate', and 'path'.
    """
    path = urllib.parse.urlsplit(file_path).path

    meta_from_name = parse_gwosc_hdf_filename(path)
    if meta_from_name is not None:
        obs, ifo, desc, start_time, duration = meta_from_name
        # Filename parsed, but we may still verify or load sample_rate on demand
        return {
            'obs': obs,
            'ifo': ifo,
            'desc': desc,
            'start_time': start_time,
            'duration': duration,
            'end_time': start_time + duration,
            'path': path,
        }

    # Fallback: inspect HDF5 attributes
    with h5py.File(path, 'r') as h5file:
        start_time = None
        duration = None
        ifo = None
        obs = None
        desc = 'GWOSC_STRAIN'

        if 'meta' in h5file:
            meta = h5file['meta']
            if 'GPSstart' in meta:
                start_time = float(meta['GPSstart'][()])
            if 'Duration' in meta:
                duration = float(meta['Duration'][()])
            if 'Detector' in meta:
                det = meta['Detector'][()]
                ifo = det.decode('utf-8') if isinstance(det, bytes) else str(det)
            if 'Observatory' in meta:
                ob = meta['Observatory'][()]
                obs = ob.decode('utf-8') if isinstance(ob, bytes) else str(ob)
            if 'FrameType' in meta:
                ft = meta['FrameType'][()]
                desc = ft.decode('utf-8') if isinstance(ft, bytes) else str(ft)

        if 'strain/Strain' in h5file:
            dset = h5file['strain/Strain']
            if start_time is None and 'Xstart' in dset.attrs:
                start_time = float(dset.attrs['Xstart'])
            if duration is None and 'Xspacing' in dset.attrs and 'Npoints' in dset.attrs:
                duration = float(dset.attrs['Npoints']) * float(dset.attrs['Xspacing'])

        if ifo is None:
            ifo = 'H1'
        if obs is None:
            obs = ifo[0]

        if start_time is None or duration is None:
            raise ValueError(f"Could not determine time span for GWOSC HDF file: {path}")

        return {
            'obs': obs,
            'ifo': ifo,
            'desc': desc,
            'start_time': start_time,
            'duration': duration,
            'end_time': start_time + duration,
            'path': path,
        }


def resolve_gwosc_hdf_channel(h5file, requested_channel):
    """Resolve a requested channel string to an internal HDF5 dataset path.

    Parameters
    ----------
    h5file : h5py.File
        Open HDF5 file handle.
    requested_channel : str
        The requested channel name (e.g. 'H1:STRAIN', 'H1:GWOSC-16KHZ_R1_STRAIN',
        'strain/Strain', or 'quality/simple/DQmask').

    Returns
    -------
    str
        Path of the dataset inside the HDF5 file to read.

    Raises
    ------
    ValueError
        If the channel cannot be resolved or is incompatible with the file.
    """
    # 1. Direct dataset path match
    norm_path = requested_channel.lstrip('/')
    if norm_path in h5file and isinstance(h5file[norm_path], h5py.Dataset):
        return norm_path

    # Extract detector prefix if present (e.g. 'H1' from 'H1:STRAIN')
    req_ifo = None
    req_name = requested_channel
    if ':' in requested_channel:
        parts = requested_channel.split(':')
        req_ifo = parts[0]
        req_name = parts[-1]

    # Validate detector compatibility
    file_ifo = None
    if 'meta/Detector' in h5file:
        det = h5file['meta/Detector'][()]
        file_ifo = det.decode('utf-8') if isinstance(det, bytes) else str(det)

    if req_ifo and file_ifo and req_ifo != file_ifo:
        raise ValueError(
            f"Requested channel '{requested_channel}' detector '{req_ifo}' "
            f"does not match file detector '{file_ifo}'"
        )

    # 2. Check meta/StrainChannel
    if 'meta/StrainChannel' in h5file:
        sc = h5file['meta/StrainChannel'][()]
        sc_str = sc.decode('utf-8') if isinstance(sc, bytes) else str(sc)
        if requested_channel == sc_str or req_name == sc_str or requested_channel == sc_str.split(':')[-1]:
            if 'strain/Strain' in h5file:
                return 'strain/Strain'

    # 3. Standard strain aliases (H1:STRAIN, GWOSC-*, HOFT, etc.)
    standard_strain_keywords = ('STRAIN', 'HOFT', 'GWOSC')
    req_upper = requested_channel.upper()
    if any(kw in req_upper for kw in standard_strain_keywords):
        if 'strain/Strain' in h5file:
            return 'strain/Strain'

    # 4. Default: If strain/Strain exists, return it
    if 'strain/Strain' in h5file:
        return 'strain/Strain'

    raise ValueError(
        f"Channel '{requested_channel}' not found or recognized in GWOSC HDF5 file. "
        f"Available datasets: {list(h5file.keys())}"
    )


def extract_files_from_cache(cache_file_path):
    """Parse a LAL cache or text file to extract file paths/URLs.

    Supports standard LAL cache format (OBS DESC START DUR URL) or
    simple line-by-line lists of file paths/URLs.

    Parameters
    ----------
    cache_file_path : str
        Path to the cache file.

    Returns
    -------
    list of str
        List of file paths/URLs.
    """
    file_list = []
    with open(cache_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # In LAL cache format, the last column is the URL/path
            url = parts[-1]
            if url.startswith('file://localhost/'):
                path = url.replace('file://localhost', '')
            elif url.startswith('file://'):
                path = url.replace('file://', '')
            else:
                path = url
            file_list.append(path)
    return file_list


def read_frame_gwosc_hdf(locations, channels, start_time=None, end_time=None,
                         duration=None, check_integrity=False, sieve=None):
    """Read TimeSeries data from GWOSC HDF5 files or cache.

    Parameters
    ----------
    locations : str or list
        Filename, glob pattern, cache file path, or list of files/caches.
    channels : str or list of str
        Channel name or list of channel names to read.
    start_time : float or int, optional
        GPS start time. If None, defaults to the start of the available data.
    end_time : float or int, optional
        GPS end time. Incompatible with `duration`.
    duration : float or int, optional
        Duration in seconds. Incompatible with `end_time`.
    check_integrity : bool, optional
        Validate file headers and dataset shapes.
    sieve : str, optional
        Regular expression to filter file URLs/paths.

    Returns
    -------
    TimeSeries or list of TimeSeries
        The requested strain time series.
    """
    if end_time is not None and duration is not None:
        raise ValueError("end_time and duration are mutually exclusive")

    # Normalize locations into a list of file paths
    raw_locations = locations if isinstance(locations, list) else [locations]
    file_paths = []

    for loc in raw_locations:
        # Handle lal.Cache or glue.lal.Cache objects if passed
        if hasattr(loc, 'tofile') or hasattr(loc, 'calc_srl'):
            # glue.lal.Cache or SWIG lal.Cache
            try:
                for entry in loc:
                    url = getattr(entry, 'path', None) or getattr(entry, 'url', None) or str(entry)
                    file_paths.append(urllib.parse.urlsplit(str(url)).path)
            except TypeError:
                import tempfile
                import lal
                with tempfile.NamedTemporaryFile('w') as tf:
                    lal.CacheExport(loc, tf.name)
                    file_paths.extend(extract_files_from_cache(tf.name))
            continue

        if not isinstance(loc, str):
            continue

        # Check if loc is a cache file (.cache, .lcf)
        loc_path = urllib.parse.urlsplit(loc).path
        _, ext = os.path.splitext(loc_path)
        if ext.lower() in ('.cache', '.lcf') and os.path.isfile(loc_path):
            file_paths.extend(extract_files_from_cache(loc_path))
        else:
            # Could be a glob or single file path
            expanded = glob.glob(loc_path)
            if expanded:
                file_paths.extend(expanded)
            else:
                file_paths.append(loc_path)

    # Remove duplicates preserving order
    unique_paths = []
    seen = set()
    for fp in file_paths:
        if fp not in seen:
            seen.add(fp)
            unique_paths.append(fp)
    file_paths = unique_paths

    if not file_paths:
        raise ValueError("No files found from the specified location(s)")

    # Apply regex sieve if requested
    if sieve:
        regex = re.compile(sieve)
        file_paths = [fp for fp in file_paths if regex.search(fp)]
        if not file_paths:
            raise ValueError(f"No files matched sieve regex '{sieve}'")

    # Index files by time coverage
    file_entries = []
    for fp in file_paths:
        meta = get_gwosc_hdf_metadata(fp)
        file_entries.append(meta)

    # Sort files by GPS start time
    file_entries.sort(key=lambda x: x['start_time'])

    # Determine global time bounds if not explicitly given
    if start_time is None:
        start_time = file_entries[0]['start_time']
    else:
        start_time = float(start_time)

    if duration is not None:
        if duration <= 0:
            raise ValueError("Duration must be strictly positive")
        end_time = start_time + float(duration)
    elif end_time is None:
        end_time = file_entries[-1]['end_time']
    else:
        end_time = float(end_time)

    if end_time <= start_time:
        raise ValueError(f"Invalid time interval: start_time ({start_time}) >= end_time ({end_time})")

    # Filter to files that overlap [start_time, end_time]
    overlapping = [
        f for f in file_entries
        if f['end_time'] > start_time and f['start_time'] < end_time
    ]

    if not overlapping:
        raise ValueError(
            f"No GWOSC HDF files cover the requested interval [{start_time}, {end_time})"
        )

    # Verify contiguous coverage across [start_time, end_time]
    # 1. Earliest file must start at or before start_time
    if overlapping[0]['start_time'] > start_time:
        raise ValueError(
            f"Data gap: requested start_time {start_time} is earlier than "
            f"earliest available data {overlapping[0]['start_time']}"
        )
    # 2. Latest file must end at or after end_time
    if overlapping[-1]['end_time'] < end_time:
        raise ValueError(
            f"Data gap: requested end_time {end_time} is later than "
            f"latest available data {overlapping[-1]['end_time']}"
        )
    # 3. Check for gaps between adjacent files
    current_covered = overlapping[0]['end_time']
    for entry in overlapping[1:]:
        if entry['start_time'] > current_covered:
            raise ValueError(
                f"Data gap detected: missing data between {current_covered} "
                f"and {entry['start_time']} in GWOSC HDF files"
            )
        current_covered = max(current_covered, entry['end_time'])

    # Prepare channels list
    is_single_channel = not isinstance(channels, list)
    channel_list = [channels] if is_single_channel else channels

    results = []

    for chan in channel_list:
        chunks = []
        delta_t = None
        data_dtype = None

        for entry in overlapping:
            fp = entry['path']
            with h5py.File(fp, 'r') as h5f:
                if check_integrity:
                    if 'strain/Strain' in h5f:
                        dset = h5f['strain/Strain']
                        if 'Npoints' in dset.attrs:
                            if len(dset) != dset.attrs['Npoints']:
                                raise ValueError(f"Corrupt HDF5 dataset length in {fp}")

                dset_path = resolve_gwosc_hdf_channel(h5f, chan)
                dset = h5f[dset_path]

                # Extract time-series attributes
                if 'Xstart' in dset.attrs:
                    f_start = float(dset.attrs['Xstart'])
                else:
                    f_start = float(entry['start_time'])

                if 'Xspacing' in dset.attrs:
                    f_dt = float(dset.attrs['Xspacing'])
                else:
                    f_dt = float(entry['duration']) / float(len(dset))

                if delta_t is None:
                    delta_t = f_dt
                    data_dtype = dset.dtype
                elif not math.isclose(delta_t, f_dt, rel_tol=1e-7):
                    raise ValueError(
                        f"Inconsistent sample rate across files: {1.0/delta_t} vs {1.0/f_dt}"
                    )

                f_len = len(dset)
                f_end = f_start + f_len * f_dt

                # Compute slice overlapping [start_time, end_time]
                sub_start = max(start_time, f_start)
                sub_end = min(end_time, f_end)

                idx_start = max(0, int(round((sub_start - f_start) / f_dt)))
                idx_end = min(f_len, int(round((sub_end - f_start) / f_dt)))

                if idx_end > idx_start:
                    # Slice on disk
                    chunk_arr = dset[idx_start:idx_end]
                    chunks.append(chunk_arr)

        if not chunks:
            raise ValueError(f"No data retrieved for channel '{chan}' in [{start_time}, {end_time})")

        concatenated = np.concatenate(chunks) if len(chunks) > 1 else chunks[0]
        ts = TimeSeries(concatenated, delta_t=delta_t, epoch=start_time)
        results.append(ts)

    if is_single_channel:
        return results[0]
    return results
