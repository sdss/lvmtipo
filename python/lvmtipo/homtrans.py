# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: homtrans.py
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


__all__ = ['HomTrans']


class HomTrans():
    """ A single affine coordinate transformation.
    Represented internally by a 4x4 matrix as a projective
    "homogeneous" transformation.
    """

    def __init__(self, entries):

        if isinstance(entries, numpy.ndarray) :
            self.matr = entries
        else :
            self.matr = numpy.array(entries, numpy.double)


    def multiply(self, rhs):
        """ 
        :param rhs The transformation to the right of the multiplication
                   sign. So rhs is applied before this transformation.
        :type HomTrans
        :return The homogeneous transformation which is the (non-communtative)
                product of self with rhs, representing the consecutive
                application of rhs, then self.
        """
        if isinstance(rhs, HomTrans) :
            prod = numpy.matmul(self.matr,rhs.matr)
            return HomTrans(prod)
        elif isinstance(rhs, numpy.ndarray) :
            prod = numpy.matmul(self.matr,rhs)
            return HomTrans(prod)
        else :
            raise TypeError("invalid data types")

    def apply(self, rhs):
        """ 
        Apply self transformation to a vector of coordinates.
        :param rhs The vector. If it has only the standard 3 coordinates,
                   a virtual 1 is appended before applying the transformation. 
        :type numpy.ndarray of dimension 1
        :return a numpy.ndarray with a vector of 3 (standard, projected) Cartesian coordinates.
        """
        if isinstance(rhs, numpy.ndarray) :
            if rhs.ndim == 1 :
                if rhs.shape[0] == 4 :
                    prod = numpy.dot(self.matr,rhs)
                elif rhs.shape[0] == 3 :
                    w = numpy.append(rhs,[1])
                    prod = numpy.dot(self.matr,w)

                prod /= prod[3]
                return numpy.array([prod[0],prod[1],prod[2]],numpy.double)
            else:
                raise TypeError("rhs not  a vector")
