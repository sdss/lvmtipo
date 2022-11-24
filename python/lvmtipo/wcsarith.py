# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-11-24
# @Filename: wcsarith.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for extended arithmetics of WCS mappers.
"""

import math
import astropy.coordinates
import astropy.units
import astropy.wcs

from lvmtipo.target import Target
from lvmtipo.site import Site
from lvmtipo.ambient import Ambient

__all__ = ['Wcsarith']


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
        """
        self.wcs = wcs

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
