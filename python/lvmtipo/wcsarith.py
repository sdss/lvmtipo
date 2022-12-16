# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-11-28
# @Filename: wcsarith.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for extended arithmetics of WCS mappers.
"""

import math
import numpy
import astropy.coordinates
import astropy.units
import astropy.wcs
import astropy.io.fits

from lvmtipo.target import Target
from lvmtipo.site import Site
from lvmtipo.ambient import Ambient

__all__ = ['Wcsarith','Affine2d']


# class Wcsarith(astropy.wcs.WCS):
class Wcsarith():
    """ A astropy.wcs.Wcs class extended for operations like multiplications, chaining...
    Usage example:

    >>> from lvmtipo.siderostat import Siderostat
    >>> from lvmtipo.site import Site
    >>> from lvmtipo.wcsarith import Wcsarith
    >>> import astropy.io.fits

    >>> fits_in="/home/mathar/lvm.sci.agcam.west_00000243.fits"
    >>> sid=Siderostat(azang=180.0)
    >>> wcs = sid.update_fits(fits_in,None)
    >>> lvmwcs = Wcsarith(wcs)
    >>> site = Site(name='MPIA')
    >>> offs = lvmwcs.offset_px_to_azalt(800,550,900,650,site)
    >>> print(offs)
    """

    def __init__(self,wcs):
        """ A sub-class of the WCS astropy mappers.
        :param wcs: Information of the pixel-to-ra/dec transformation
        :type wcs: either a astropy.wcs.WCS right away or the HDU list
                 obtained by open()ing a FITS file which has the information
                 in the primary header.
        """

        if isinstance(wcs, astropy.wcs.WCS):
            self.wcs = wcs
        elif isinstance(wcs, astropy.io.fits.HDUList):
            # assume that the WCS is in the primary header data unit (PHDU)
            # Not tested on tile-compressed cases where the info might be in the second HDU
            phdu = wcs[0].header
            self.wcs = astropy.wcs.WCS(phdu)
        else:
            raise TypeError("invalid wcs data type")

    def __str__(self) :
        """
        Produce a human readable representation
        :return: same as the astropy version.
        :rtype: str
        """
        return str(self.wcs)

    def __truediv__(self, rhs) :
        """ Construct the operation of applying rhs then the inverse of this.
            Warning: this operation is not commutative; swapping
            the role of self and rhs gives the inverse result.
            Application: let self be a WCS which maps camera pixels via a WCS (for example
            obtained by astrometry) and
            let rhs be a different WCS which also maps camera pixels via a WCS (for example
            obtained by a model of the designed camera location), then the result is the
            affine transformation which maps the two pixel systems onto each other, building
            a "calibrated" location of the camera "as build" relative to the "as designed".
            The returned transformation applied to the rhs pixels will give the pixels of
            the self system. Note that such a mapping of pixels to pixels will also 
            include any errors of the pointing model assumed in the WCS(s), not just geometric
            distortions that follow from imperfect alignments of camera chips.
            A potential application could be a guider that has obtained for two
            images two WCS maps and wants to convert the "difference" in these two
            maps to a transformation between pixels.
        :return: an Affine2d object which correlates the two maps in their image coordinates.
        :rtype: lvmtipo.Affine2d
        """
        aff = Affine2d(wcs=self.wcs)
        aff_rhs = Affine2d(wcs=rhs.wcs)
        return aff / aff_rhs

    def offset_px_to_azalt(self, old_x, old_y, new_x, new_y, site, ambi=None, time=None):
        """
        Given an old pixel position (old_x,old_y) and a new pixel position (new_x,new_y)
        in a FITS coordinate system,
        compute the associated change of the reference pointing coordinate in
        azimuth and alitude of the horizontal coordinate system.
        This answers the question: if an object is at some (old_x,old_y) coordinate system,
        which offset is needed so it will appear at (new_x,new_y) after sending the offset?

        :param old_x: pixel first coordinate in FITS conventions
        :type old_x: int or float
        :param old_y: pixel second coordinate in FITS conventions
        :type old_y: int or float
        :param new_x: desired new first coordinate
        :type new_x: int or float
        :param new_y: desired new second coordinate
        :type new_y: int or float

        :param site: Observatory location
        :type site: lvmtipo.Site

        :param ambi: Ambient parameters used for refraction correction
        :type ambi: lvmtipo.Ambient

        :param time: time of the observation
        :type time: astropy.time.Time

        :return: a pair of two angles, first the change in azimuth (0=N, 90=E) then the 
             change in altitude (0=horizon)
        :rtype: astropy.coordinates.Angle pair
        """

        # compute where the old equatorial ra/dec are
        radec = self.wcs.all_pix2world([old_x],[old_y],1,ra_dec_order=True)
        object = astropy.coordinates.SkyCoord(ra=radec[0], dec=radec[1],unit="deg")
        old_targ = Target(object)

        # compute where the new  equatorial ra/dec are
        radec = self.wcs.all_pix2world([new_x],[new_y],1,ra_dec_order=True)
        object = astropy.coordinates.SkyCoord(ra=radec[0], dec=radec[1],unit="deg")
        new_targ = Target(object)

        if isinstance(time, astropy.time.Time):
            now = time
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()

        if isinstance(site, Site):
            if ambi is None:
                refr = Ambient(site=site)
            elif isinstance(ambi, Ambient):
                refr = ambi
            else:
                raise TypeError("invalid ambi data type")
            earthloc = site.toEarthLocation()
        else:
            raise TypeError("invalid site  data type")

        # compute where the old pixels are in alt/az 
        oldAa = old_targ.toHoriz(site,ambi=ambi,time=now)

        # compute where the new pixels are in alt/az 
        newAa = new_targ.toHoriz(site,ambi=ambi,time=now)

        # to be checked: is this the subtraction with the correct sign/order???
        # or does the PW mount want the opposite signs??
        # We assume here that the offset is the change of the pointing
        # direction of the mount. E.g. if we sent +1deg in alt then
        # the mount moves horizontally to the right. The effect of
        # the targets is the opposite, which means they appear more
        # to the left in the new coordinate system.
        delta_alt = oldAa.alt - newAa.alt

        delta_az = oldAa.az - newAa.az

        return delta_az, delta_alt

    def offset_px_to_radec(self, old_x, old_y, new_x, new_y):
        """
        Given an old pixel position (old_x,old_y) and a new pixel position (new_x,new_y),
        compute the associated change of the reference pointing coordinate in
        equatorial coordinates (ra/dec) as two differential angles.

        :param old_x: pixel first coordinate in FITS conventions
        :type old_x: int or float
        :param old_y: pixel second coordinate in FITS conventions
        :type old_y: int or float
        :param new_x: desired new first coordinate
        :type new_x: int or float
        :param new_y: desired new second coordinate
        :type new_y: int or float

        :return: a pair of two Angles, first the change in ra then the change in dec
        :rtype: astropy.coordinates.Angle pair
        """

        # compute where the old equatorial ra/dec are
        old_radec = self.wcs.all_pix2world([old_x],[old_y],1,ra_dec_order=True)
        old_object = astropy.coordinates.SkyCoord(ra=old_radec[0], dec=old_radec[1],unit="deg")

        # compute where the new  equatorial ra/dec are
        new_radec = self.wcs.all_pix2world([new_x],[new_y],1,ra_dec_order=True)
        new_object = astropy.coordinates.SkyCoord(ra=new_radec[0], dec=new_radec[1],unit="deg")

        # one could (easier, faster) just return old_radec[0]-new_radec[0]
        # and old_radec[1] - new_radec[1] here, converted to the Angle class...

        # To be checked: is this the subtraction with the correct sign/order???
        # or does the PW mount want the opposite signs??
        # This here constructs old_object minus new_object, so uses the same
        # sign convention as in offset_px_toazalt().

        return new_object.spherical_offsets_to(old_object)

class Affine2d():
    """ Affine transformations in 2d: planar rotations followed by shifts.
    """
    def __init__(self,rot_ang=0.0, scale=1.0, shift_x=0.0, shift_y=0.0, 
           mat=None, trans=None, wcs=None):
        """ A transformation in a plane (pixels, focalplane...)
        represented by a 2x2 numpy array for a shift-rotation and a
        2x1 numpy array for a translation. Note that the shift-rotation
        and scale are applied before the translation; these are not commutative.

        :param rot_ang: the rotation angle of the rotation in degrees
              or as an astropy.coordinates.Angle. A positive angle
              rotates objects clockwise.
        :type rot_ang: float or astropy.coordinates.Angle

        :param scale: the magnification factor for input values.
              This means if the transformation is applied to a pair
              of points of distance 1, they will have distance scale
              after the transformation.
              Note that negative scales are equivalent to introducing
              flips.
        :type scale: float

        :param shift_x: the shift along the first coordinate
             applied after the rotation. This is already in the up-scaled
             units. If e.g. scale=2 (representing 2 um) and one wants
             to move a point 2 um to the right, shift_x=2.
        :type shift_x: float

        :param shift_y: the shift along the second coordinate
             applied after the rotation. This is already in the up-scaled
             units. 
        :type shift_y: float

        :param mat: the rotation-shear matrix.
             If None, the rot_ang and scale parameters will be evaluated, 
             otherwise this is a substitue for a copy-constructor. It's a design
             flaw of python that function overloading is not supported.
        :type mat: 2x2 numpy.array

        :param trans: the translation vector in x and y.
             If None, shift_x and shift_y parameters will be evaluated, 
             otherwise works as a copy-constructor.
        :type trans: 2x1 numpy.array

        :param wcs: a wcs model that transforms pixels to world coordinates.
             If not None, all the other parameters are ignored. (It's a design
             error of python that function overloading is not supported.)
        :type wcs: astropy.wcs.WCS

        """
        if wcs is None :
            if trans is None :
                self.shift = numpy.array([shift_x, shift_y],numpy.double)
            elif isinstance(trans, numpy.ndarray) :
                self.shift = trans
            else :
                raise TypeError("trans not numpy array")
    
            if mat is None :
                if isinstance(rot_ang, (int,float)) :
                    ang=math.radians(rot_ang)
                elif isinstance(rot_ang, astropy.coordinates.Angle):
                    ang = rot_ang.rad
                else:
                    raise TypeError("invalid rot_ang data type")
        
                cosang= scale*math.cos(ang)
                sinang= scale*math.sin(ang)
        
                self.rot = numpy.empty(shape=(2,2))
                self.rot[0][0] = self.rot[1][1] = cosang 
                self.rot[0][1] = sinang
                self.rot[1][0] = -sinang
            elif isinstance(mat, numpy.ndarray) :
                self.rot = mat
            else :
                raise TypeError("mat not numpy array")
        else:
            if isinstance(wcs, astropy.wcs.WCS):
                # try to reconstruct the (approximate) linear
                # transformation by applying at least 3 test
                # transformations of the wcs.
                # To figure out the shift/translate part it suffices
                # to feed (0,0) into the wcs because the
                # self.rot part becomes ineffective then.
                p0_radec = wcs.all_pix2world([0],[0],1,ra_dec_order=True)
                # now p0_radec[0] is ra and p0_radec[1] is dec in degrees

                # print("p0 is",p0_radec)
                # pix2world returns arrays, so we have to extract 2-dim components here...
                self.shift = numpy.array([p0_radec[0][0], p0_radec[1][0]],numpy.double)

                # the matrix will be 2x2
                self.rot = numpy.empty(shape=(2,2))

                # Uses a sort of linearization throw to poke the WCS map.
                # Should be some positive integer, not too small to reduce
                # numerical jitter
                pix_throw = 100
                # get image of pixel (pix_throw,0) in ra/dec
                p1_radec = wcs.all_pix2world([pix_throw],[0],1,ra_dec_order=True)
                #print("p1 is",p1_radec)
                # first column of the matrix is essentially that image modulo the shift
                # p1_radec[0] -= p0_radec[0]
                # p1_radec[1] -= p0_radec[1]
                # print("shifted p1 is",p1_radec)
                self.rot[0][0] = (p1_radec[0][0] - p0_radec[0][0])/pix_throw
                self.rot[1][0] = (p1_radec[1][0] - p0_radec[1][0])/pix_throw

                # get image of pixel (0,pix_throw) in ra/dec
                p1_radec = wcs.all_pix2world([0],[pix_throw],1,ra_dec_order=True)
                # print("p1 is",p1_radec)
                # second column of the matrix is essentially that image modulo the shift
                # p1_radec[0] -= p0_radec[0]
                # p1_radec[1] -= p0_radec[1]
                # print("shifted p1 is",p1_radec)
                self.rot[0][1] = (p1_radec[0][0] - p0_radec[0][0])/pix_throw
                self.rot[1][1] = (p1_radec[1][0] - p0_radec[1][0])/pix_throw
            else:
                raise TypeError("wcs not an astropy.wcs.WCS")

        # end of __init__()

    def __str__(self) :
        """
        Produce a human readable representation
        :return: The 2x2 matrix and the 2x1 shift vector separated by blank
        :rtype: str
        """
        return str(self.rot) + " " + str(self.shift)

    def __inv__(self) :
        """ Construct the inverse transformation.
        :return: an Affine2d object which unapplies the transformation.
        :rtype: lvmtipo.Affine2d
        """
        rotinv = numpy.linalg.inv(self.rot)
        shifti = numpy.matmul(rotinv,self.shift)
        shifti *= -1
        return Affine2d(mat= rotinv, trans=shifti)

    def __mul__(self, rhs) :
        """ Construct the operation of first applying rhs then self.
            Warning: this operation is not commutative; the rhs is
            applied first and the right hand side in usual operator notation.
        :return: an Affine2d object which runs the two maps one after the other
        :rtype: lvmtipo.Affine2d
        """
        rotprod = numpy.matmul(self.rot,rhs.rot)
        shiftprod = numpy.matmul(self.rot,rhs.shift)
        shiftprod += self.shift
        return Affine2d(mat= rotprod, trans=shiftprod)

    def __truediv__(self, rhs) :
        """ Construct the operation of applying rhs then the inverse of this.
            Warning: this operation is not commutative; swapping
            the role of self and rhs gives the inverse result.
        :return: an Affine2d object which runs the two maps one after the other
        :rtype: lvmtipo.Affine2d
        """
        inv = self.__inv__()
        return inv * rhs

    def apply(self, xy_in) :
        """ Apply the transformation to a point
        :param xy_in: the x and y coordinate of the point.
        :type xy_in: numpy array 

        :return: the new 2 coordinates of the image of the point 
            after applying the transf.
        :rtype: numpy array 
        """
        roted = numpy.matmul(self.rot, xy_in)
        return numpy.add(roted,self.shift)

    def get_ang_scale(self) :
        """ extract angle and scale factor.
        If the transformation has no shear component and is just
        a pure scale-rotation and a shift, the results are exact.
        :return: two values for angle and scale factor.
                 If the initialization was derived from a WCS with pixels to degrees,
                 the scale factor has the units deg/px.
        :rtype: astropy.coordinates.Angle, float
        """
        # determinant (may be negative pending on flips....)
        # Because this is basically squared entries, we pull the sqroot to return
        det = abs(numpy.linalg.det(self.rot))
        # print("det",det)
        # entries of first colun are cosine and negated sign
        ang = math.atan2(- self.rot[1][0],self.rot[0][0])
        ang = astropy.coordinates.Angle(ang,astropy.units.radian)
        return ang, math.sqrt(det)
