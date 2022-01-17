# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: site.py
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


__all__ = ['Site']


class Site():
    """ Geolocation of observatory
    """
    def __init__(self, long= -70.70056, lat= -29.01091, alt=2280, name=None) :
        """ Geolocation of observatory
        :param long geodetic longitude in degrees, E=positive
               Default is the LCO parameter.
        :type long float
        :param lat geodetic latitude in degrees, N=positive
               Default is the LCO parameter.
        :type lat float
        :param alt altitude above sea level
               Default is the LCO parameter.
        :type alt float
        :param name of one of the 4 LVM site acronyms, {LCO|APO|MPIA|KHU}
               If this parameter is present, it overrides the values
               of the other 3 numerical parameters.
        :type name string
        """

        if name is not None:
            if name == 'LCO' :
                self.long = -70.70056
                self.lat = -29.01091
                self.alt = 2280.
            elif name == 'APO' :
                self.long = -105.8202778 # -105d49m13s
                self.lat =  32.78028 # 32d46m49s
                self.alt = 2788.
            elif name == 'MPIA' :
                self.long = 8.724
                self.lat =  49.3965
                self.alt = 560.
            elif name == 'KHU' :
                self.long = 127.0533
                self.lat =  37.5970
                self.alt = 80.
        elif isinstance(long, (int,float) ) and isinstance(lat, (int,float)) and isinstance(alt,(int,float)) :
            self.long = long
            self.lat = lat
            self.alt = alt
        else :
            raise TypeError("invalid data types")

        # print("site" +str(self.long)+ "deg " + str(self.lat) + "deg")

    def toEarthLocation(self) :
        return astropy.coordinates.EarthLocation.from_geodetic(self.long, self.lat, height=self.alt)
