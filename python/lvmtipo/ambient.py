# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: ambient.py
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

from .site import Site

__all__ = ['Ambient']


class Ambient():
    """ ambient parameters relevant to atmospheric refraction
    :param press pressure in hPa. Can be None if the site parameter
                 is not-None so we can get a pressure estimate from sea level altitude.
    :type press float
    :param temp Temperature in deg Celsius.
    :type temp float
    :param rhum relative humidity in the range 0..1.
    :type rhum float
    :param wlen Observing wavelength in microns.
    :type wlen float
    :param site location of observatory
    :type site fieldrotation.Site
    """
    def __init__(self, press = None, temp = 15, rhum = 0.2, wlen=0.5,  site=None) :
        """ 
        """
        self.temp = temp
        self.rhum = rhum
        self.wlen = wlen
        if press is None :
            if site is None :
                # no information whatsoever: assume sea level
                self.press = 1013.0
            elif isinstance(site, Site) :
                # assume 8.1 km scale height
                self.press = 1013.0*math.exp(-site.alt/8135.0)
            else :
                self.press = 1013.0
        else:
            self.press = press

