# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: mirror.py
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

from .homtrans import HomTrans

__all__ = ['Mirror']


class Mirror():
    """ A flat mirror.
    This represents an infintely large flat plane. The internal
    representation is the surface normal and the standard 
    equation that the dot product of points on the surface by
    the surface normal equals the distance (of the plane to the
    origina of coordinates).
    """

    def __init__(self, normal, disttoorg):
        """ 
        :param normal The 3 Cartesian coordinates of the surface normal.
               It must have nonzero length, but
               does not need to be normalized to unit length.
        :type numpy.ndarray with 3 (xyz) numbers
        :param disttoorg The distance of the mirror to the origin of coordinates
               As in usual geometry, the distance is the shortest distance of the origin
               to the infinitely extended mirror plane.
        :type float
        """
        if isinstance(normal, numpy.ndarray) and isinstance(disttoorg, (int,float)) :
            self.d = float(disttoorg)
            if normal.ndim == 1 and normal.shape[0] == 3:
                len = numpy.linalg.norm(normal)
                normal /= len
                self.n = normal
            else :
                raise TypeError("invalid data types")
        else :
            raise TypeError("invalid data types")

    def toHomTrans(self) :
        """ 
        :return The homogeneous transformation that represents
              the reflection of rays off this mirror surface.
        """ 
        matr = numpy.zeros((4,4))
        for r in range(4):
            for c in range(4):
                if r == c :
                    matr[r,c] += 1.0
                if r < 3 :
                    if c < 3 :
                        matr[r,c] -= 2.0*self.n[r]*self.n[c]
                    else :
                        matr[r,c] = 2.0*self.d*self.n[r]
        return HomTrans(matr)
