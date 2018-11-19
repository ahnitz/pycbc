# Copyright (C) 2018 Alex Nitz

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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" Utitlities to query the status of gravitational-wave instruments
from public sources as well as dqsegdb.
"""

from glue.segments import segmentlist
import json, urllib

def query_flag(ifo, segment_name, start_time, end_time):
    """Return the times where the flag is active 

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    segment_name: string
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    ### Special cases as the LOSC convention is backwards from normal
    ### LIGO / Virgo operation!!!!
    if (('_HW_INJ' in segment_name and not 'NO' in segment_name) or 
       'VETO' in segment_name):
        data = query_flag(ifo, 'DATA', start_time, end_time)
         
        if '_HW_INJ' in segment_name:
            name = 'NO_' + segment_name    
        else:
            name = segment_name.replace('_VETO', '')

        negate = query_flag(ifo, name, start_time, end_time)
        return (data - negate).coalesce()

    duration = end_time - start_time
    url = 'https://www.gw-openscience.org/timeline/segments/json/O1/{}_{}/{}/{}/'
    url = url.format(ifo, segment_name, start_time, duration)
    
    try:
        segments = json.loads(urllib.urlopen(url).read())['segments']
    except:
        raise ValueError('Unable to find segments, check flag name or times')

    return segmentlist(segments)

def query_cumulative_flags(ifo, segment_names, start_time, end_time):
    """Return the times where any flag is active 

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    segment_name: list of strings
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    total_segs = None
    for flag_name in segment_names:
        segs = query_flag(ifo, flag_name, start_time, end_time)
        if total_segs is not None:
            total_segs = (total_segs + segs).coalesce()
        else:
            total_segs = segs
    return total_segs

def query_combined_flags(ifo, up_flags, start_time, end_time, down_flags=None):
    """Return the times where any "up" flag is active minus "down" flag times

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    up flags: list of strings
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query
    down flags: list of strings
        Flags which indicate times to subtract from the combined segments

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    segments = query_cumulative_flags(ifo, up_flags, start_time, end_time)
    down_flags = [] if down_flags is None else down_flags
    for flag_name in down_flags:
        mseg = query_flag(ifo, flag_name, start_time, end_time)
        segments = (segments - mseg).coalesce()
    return segments
