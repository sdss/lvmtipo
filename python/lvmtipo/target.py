# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: target.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for siderostat field angles using homogeneous coordinates
"""

import sys
import math
import numpy
import astropy.coordinates
import astropy.time
import astropy.units

from .ambient import Ambient
from .site import Site

__all__ = ['Target']

class Target():
    """ sidereal astronomical target
    """
    def __init__(self, targ) :
        """ target coordinates
        :param targ Position in equatorial coordinates
        :type targ astropy.coordinates.Skycoord
        """

        if isinstance(targ, astropy.coordinates.SkyCoord) :
            self.targ = targ
        else :
            raise TypeError("invalid data types")
        # print(self.targ)

    def toHoriz(self, site, ambi = None, time = None) :
        """ convert from equatorial to horizontal coordinates
        :param site Observatory location
        :type site fieldrotation.Site
        :param ambi Ambient parameters used for refraction correction
        :type ambi fieldrotation.Ambient
        :param time time of the observation
        :type time astropy.time.Time
        :return alt-az coordinates
        :return type astropy.coordinates.AltAz
        """
        if isinstance(time, astropy.time.Time) :
            now = time 
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()

        if isinstance(site, Site) :
            if ambi is None:
               refr = Ambient(site = site)
            elif isinstance(site, Ambient) :
               refr = ambi
            else :
                raise TypeError("invalid ambi data type")
            earthloc = site.toEarthLocation()
        else :
            raise TypeError("invalid site  data type")

        # print(earthloc)
        # print(astropy.units.Quantity(100.*refr.press,unit=astropy.units.Pa))
        # print(astropy.units.Quantity(refr.wlen,unit= astropy.units.um))
        # todo: render also proper motions (all 3 coords)
        # This is a blank form of Alt/aZ because the two angles are yet unknown
        # altaz = astropy.coordinates.builtin_frames.AltAz
        altaz = astropy.coordinates.AltAz(
                 location = earthloc,
                 obstime=now,
                 pressure= astropy.units.Quantity(100.*refr.press,unit=astropy.units.Pa),
                 temperature = astropy.units.Quantity(refr.temp, unit = astropy.units.deg_C),
                 relative_humidity = refr.rhum,
                 obswl = astropy.units.Quantity(refr.wlen,unit= astropy.units.um))

        try:
            horiz = self.targ.transform_to(altaz)

        except ValueError as ex:
            # This is sometimes triggered by being offline or the
            # IERS data server being unreachable.
            # Try again with a sort of offline attempt of the IERS tables
            from astropy.utils.iers import conf
            conf.auto_download = False
            horiz = self.targ.transform_to(altaz)

        return horiz

      
