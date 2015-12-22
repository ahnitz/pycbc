#!/usr/bin/env python

# Copyright (C) 2011 Duncan M. Macleod
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

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import numpy as np,re,sys

import lal
from lal import CachedDetectors as _lalCachedDetectors

from glue.ligolw import table

cached_detector_by_name = dict((cd.frDetector.name, cd) for cd in _lalCachedDetectors)
cached_detector = cached_detector_by_name

name_to_prefix = dict((name, detector.frDetector.prefix) for name, detector in cached_detector_by_name.items())
prefix_to_name = dict((prefix, name) for name, prefix in name_to_prefix.items())

# =============================================================================
# Sky Positions structure
# =============================================================================

class SkyPositionTable(table.Table):
    """
    glue.ligolw.table.Table holding SkyPosition objects.
    """

    tableName     = "sky_positions:table"
    valid_columns = {"process_id":  "ilwd:char",\
                     "latitude":    "real_4",\
                     "longitude":   "real_4",\
                     "probability": "real_4",\
                     "solid_angle": "real_4",\
                     "system":      "lstring"}

    def rotate(self, R):
        rotated = table.new_from_template(self)
        for row in self:
            rotated.append(row.rotate(R))
        return rotated

    def normalize(self):
        for row in self:
            row.normalize()



class SkyPosition(object):

    __slots__ = SkyPositionTable.valid_columns.keys()

    def __iter__(self):
        return iter(tuple((self.longitude, self.latitude)))

    def __repr__(self):
        return "SkyPosition(%s, %s)" % (self.longitude, self.latitude)

    def __str__(self):
        return "(%s, %s)" % (self.longitude, self.latitude)

    def rotate(self, R):
        """
        Apply the 3x3 rotation matrix R to this SkyPosition.
        """

        cart = SphericalToCartesian(self)
        rot  = np.dot(np.asarray(R), cart)
        new  = CartesianToSpherical(rot, system=self.system)
        new.normalize()
        try:
            new.probability = self.probability
        except AttributeError:
            new.probability = None
        try:
            new.solid_angle = self.solid_angle
        except AttributeError:
            new.solid_angle = None
            
        return new

    def normalize(self):
        """
        Normalize this SkyPosition to be within standard limits
        [0 <= longitude < 2pi] and [-pi < latitude <= pi]
        """

        # first unwind positions, i.e. make sure that:
        # [0 <= alpha < 2pi] and [-pi < delta <= pi]

        # normalise longitude
        while self.longitude < 0:
            self.longitude += lal.TWOPI
        while self.longitude >= lal.TWOPI:
            self.longitude -= lal.TWOPI

        # normalise latitude
        while self.latitude <= -lal.PI:
            self.latitude += lal.TWOPI
        while self.latitude > lal.TWOPI:
            self.latitude -= lal.TWOPI

        #
        # now get latitude into canonical interval [-pi/2, pi/2)
        #

        if self.latitude > lal.PI_2:
            self.latitude = lal.PI - self.latitude
            if self.latitude < lal.PI:
                self.longitude += lal.PI
            else:
                self.longitude -= lal.PI

        if self.latitude < -lal.PI_2:
            self.latitude = -lal.PI - self.latitude
            if self.longitude < lal.PI:
                self.longitude += lal.PI
            else:
                self.longitude -= lal.PI

    def get_time_delay(self, ifo1, ifo2, gpstime):
        """
        Return the time delay for this SkyPosition between arrival at ifo1
        relative to ifo2, for the given gpstime.
        """
        det1 = cached_detector.get(prefix_to_name[ifo1])
        det2 = cached_detector.get(prefix_to_name[ifo2])

        return lal.ArrivalTimeDiff(det1.location, det2.location,\
                                        self.longitude, self.latitude,\
                                        lal.LIGOTimeGPS(gpstime))


# =============================================================================
# Convert between coordinate systems
# =============================================================================

# a few last definitions (taken from SkyCoordinates.c)
lal.LAL_ALPHAGAL = 3.366032942
lal.LAL_DELTAGAL = 0.473477302
lal.LAL_LGAL     = 0.576


def GeographicToEquatorial(input, gpstime):
    """
    Convert the SkyPosition object input from the inherited 'geographical'
    system to  to 'equatorial'.
    """

    if input.system.lower() != 'geographic':
        raise AttributeError("Input system!='geographic'")

    # get GMST
    gmst = np.fmod(lal.GreenwichSiderealTime(gpstime, 0),\
                   lal.TWOPI)

    # output
    output = SkyPosition()
    output.system = 'equatorial'
    output.longitude = np.fmod(input.longitude + gmst, lal.TWOPI)
    output.latitude = input.latitude
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def EquatorialToGeographic(input, gpstime):
    """
    Convert the SkyPosition object input from the inherited 'equatorial'
    system to  to 'geographic'.
    """

    if input.system.lower() != 'equatorial':
        raise AttributeError("Input system is not 'equatorial'")

    # get GMST
    gmst = np.fmod(lal.GreenwichSiderealTime(gpstime, 0),\
                   lal.TWOPI)

    # output
    output = SkyPosition()
    output.system = 'geographic'
    output.longitude = np.fmod(input.longitude - gmst, lal.TWOPI)
    output.latitude = input.latitude
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def ISOTimeDelayLine(ifos, ra, dec, gpstime=None, n=3):
    """
    Returns the n-point SkyPositionTable describing a line of constant time
    delay through the given ra and dec. for the given 2-tuple ifos.
    """

    if gpstime:
        gpstime = lal.LIGOTimeGPS(gpstime)

    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec
    p.system    = 'equatorial'
    p.probability = None
    p.solid_angle = None
    if gpstime:
        p = EquatorialToGeographic(p, gpstime)
    cart = SphericalToCartesian(p)

    # get baseline
    detectors = [cached_detector.get(prefix_to_name[ifo])\
                 for ifo in ifos]
    baseline = detectors[0].location - detectors[1].location
    baseline = baseline / np.linalg.norm(baseline)

    # get angle to baseline
    angle = np.arccos(np.dot(cart, baseline))

    # get evenly spaced ring over north pole
    longitudes = np.linspace(0, lal.TWOPI, n, endpoint=False)
    latitudes  = [lal.PI_2-angle]*len(longitudes)
    # get axis and angle of rotation
    north = np.array([0,0,1])
    axis  = np.cross(north, baseline)
    axis  = axis / np.linalg.norm(axis)
    angle = np.arccos(np.dot(north, baseline))
    R     = _rotation(axis, angle)

    # construct sky table
    iso = SkyPositionTable()
    for lon,lat in zip(longitudes, latitudes):
        e             = SkyPosition()
        e.longitude   = lon
        e.latitude    = lat
        e.probability = None
        e.solid_angle = None
        e.system      = 'geographic'
        e.normalize()
        e             = e.rotate(R)
        if gpstime:
            e         = GeographicToEquatorial(e, gpstime)
        iso.append(e)

    return iso


def MaxTimeDelayLine3(ifo1, ifo2, ra, dec, gpstime=None, n=3):
    """
        Alternative implementation of MaxTimeDelayLine.
    """

    if gpstime:
        gpstime = lal.LIGOTimeGPS(gpstime)

    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec
    p.system    = 'equatorial'
    p.probability = None
    p.solid_angle = None
    if gpstime:
        p = EquatorialToGeographic(p, gpstime)
    cart = SphericalToCartesian(p)

    # get baseline
    detectors = [cached_detector.get(prefix_to_name[ifo])\
                 for ifo in [ifo1,ifo2]]
    baseline = detectors[0].location - detectors[1].location
    baseline = baseline / np.linalg.norm(baseline)

    # get angular spacing
    dtheta = lal.TWOPI/n

    # get axis and angle of rotation
    north = np.array([0,0,1])
    axis  = np.cross(cart, baseline)
    axis  = axis / np.linalg.norm(axis)
    R     = _rotation(axis, dtheta)

    # set up list
    l = SkyPositionTable()
    # append central point
    l.append(p)

    for i in xrange(1, n):
        l.append(l[i-1].rotate(R))

    if gpstime:
        for i,p in enumerate(l):
            l[i] = GeographicToEquatorial(p, gpstime)
    return l
            
# =============================================================================
# Helpers
# =============================================================================

def SphericalToCartesian(skyPoint):
    """
    Convert SkyPosition object into Cartesian 3-vector
    """ 
 
    p     = np.zeros(3)
    phi   = skyPoint.longitude
    theta = lal.PI_2 - skyPoint.latitude
    a     = np.sin(phi)
    b     = np.cos(phi)
    c     = np.sin(theta)
    d     = np.cos(theta)

    p[0] = c*b
    p[1] = c*a
    p[2] = d

    return p 

def CartesianToSpherical(x, system='equatorial'):
    """
    Convert 3-tuple Cartesian sky position x to SkyPosition object in the
    given system
    """

    assert len(x)==3

    p = SkyPosition()
    p.system = system
    p.longitude   = np.arctan2(x[1], x[0])
    p.latitude    = lal.PI_2 - np.arccos(x[2])
    p.probability = None
    p.solid_angle = None
    p.normalize()

    return p

def _rotation(axis, angle):
    """
    Form 3x3 rotation matrix to rotate about a given 3-tuple axis by a given
    angle
    """

    R = np.zeros((3,3))

    ux,uy,uz = axis
    c        = np.cos(angle)
    s        = np.sin(angle)

    R[0][0] = c + ux**2*(1-c)
    R[0][1] = ux*uy*(1-c) - uz*s
    R[0][2] = ux*uz*(1-c) + uy*s
    
    R[1][0] = uy*ux*(1-c) + uz*s
    R[1][1] = c + uy**2*(1-c)
    R[1][2] = uy*uz*(1-c) - ux*s

    R[2][0] = uz*ux*(1-c) - uy*s
    R[2][1] = uz*uy*(1-c) + ux*s
    R[2][2] = c + uz**2*(1-c)

    return R



