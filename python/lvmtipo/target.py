# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-11-29
# @Filename: target.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for standard astropy equatorial coordinates with a richer interface
"""

import math
import astropy.coordinates
import astropy.time
import astropy.units

from .ambient import Ambient
from .site import Site

__all__ = ['Target']


class Target():
    """ sidereal astronomical target
    """

    def __init__(self, targ, pmra =0.0, pmdec=0.0):
        """ target coordinates
        :param targ: Position in equatorial coordinates
        :type targ: astropy.coordinates.Skycoord

        :param pmra: proper motion in right ascension in mas/yr
             The numerical value is equivalent to delta in alpha multiplied
             by cos(delta), as in most modern catalogues.
        :type pmra: float

        :param pmdec: proper motion in declination in mas/yr
        :type pmdec: float
        """

        if isinstance(targ, astropy.coordinates.SkyCoord):
            self.targ = targ
        else:
            raise TypeError("invalid data types")

        self.pmra = pmra
        self.pmdec = pmdec
        self.j2000 = astropy.time.Time('2000-01-01T00:00:00',format='isot',scale='utc')

    def push_prop(self,pmra=None, pmdec=None,time=None) :
        """ target coordinates
        :param pmra: proper motion in right ascension in mas/yr
             The numerical value is equivalent to delta in alpha multiplied
             by cos(delta), as in most modern catalogues.
             If None, the value at initialization time is used
        :type pmra: float

        :param pmdec: proper motion in declination in mas/yr
             If None, the value at initialization time is used
        :type pmdec: float

        :return: a new target moved to the new time but with zero proper motion
        :rtype: astropy.coordinates.Skycoord
        """
        if isinstance(time, astropy.time.Time):
            now = time
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()

        if pmra is None:
           pm_ra = self.pmra
        else:
           pm_ra = pmra

        if pmdec is None:
           pm_dec = self.pmdec
        else:
           pm_dec = pmdec

        if pm_ra != 0.0  or pm_dec != 0.0:
            deltayrs = now-self.j2000
            # deltayrs = deltayrs.to_value(format='jyear') # should work in astropy 5.1.1
            deltayrs = deltayrs.to_value(format='jd')/365.25
            # print('years',deltayrs)
            #measure both relative positions in arcseconds
            deltara = astropy.coordinates.Angle(pm_ra*deltayrs/1000.0, unit=astropy.units.arcsec)
            deltadec = astropy.coordinates.Angle(pm_dec*deltayrs/1000.0, unit=astropy.units.arcsec)

            # reinstall the cos(delta) factor for RA
            # because spherical_offsets_by() does not recognize it ?
            # No, apparently that's already done within spherical_offsets_by()!
            # deltara /= math.cos(self.targ.dec.radian)
            targnow = self.targ.spherical_offsets_by(deltara, deltadec)
        else:
            targnow = self.targ

        return targnow

    def toHoriz(self, site, ambi=None, time=None):
        """ convert from equatorial to horizontal coordinates
        :param site: Observatory location
        :type site: lvmtipo.Site

        :param ambi: Ambient parameters used for refraction correction
        :type ambi: lvmtipo.Ambient

        :param time: time of the observation
        :type time: astropy.time.Time

        :return: alt-az coordinates
        :rtype: type astropy.coordinates.AltAz
        """
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

        # print(earthloc)
        # print(astropy.units.Quantity(100.*refr.press,unit=astropy.units.Pa))
        # print(astropy.units.Quantity(refr.wlen,unit= astropy.units.um))
        # todo: render also proper motions (all 3 coords)
        # This is a blank form of Alt/aZ because the two angles are yet unknown
        # altaz = astropy.coordinates.builtin_frames.AltAz
        altaz = astropy.coordinates.AltAz(
            location=earthloc,
            obstime=now,
            pressure=astropy.units.Quantity(
                100.*refr.press, unit=astropy.units.Pa),
            temperature=astropy.units.Quantity(
                refr.temp, unit=astropy.units.deg_C),
            relative_humidity=refr.rhum,
            obswl=astropy.units.Quantity(refr.wlen, unit=astropy.units.um)
            )

        # move target by proper motions to now.
        targnow = self.push_prop(pmra = self.pmra, pmdec=self.pmdec,time=now)
        
        try:
            horiz = targnow.transform_to(altaz)

        except ValueError :
            # This is sometimes triggered by being offline or the
            # IERS data server being unreachable.
            # Try again with a sort of offline attempt of the IERS tables
            from astropy.utils.iers import conf
            conf.auto_download = False
            horiz = targnow.transform_to(altaz)

        return horiz

    def parallact(self, site, ambi, time):
        """
        Compute the parallactic angle for that target at that time.
        :param site: geographic location of the observatory
        :type site: lvmtipo.Site

        :param ambi: Ambient parameters used for refraction correction
        :type ambi: lvmtipo.Ambient

        :param time: time of the observation
        :type time: astropy.time.Time

        :return: The parallactic angle
        :rtype: astropy.coordinates.Angle
        """
        #define the zenith in the topocentric horizontal frame (alt,az)=(90,0)
        earthloc = site.toEarthLocation()
        zeni_hori = astropy.coordinates.AltAz(
            location=earthloc,
            obstime=time,
            pressure=astropy.units.Quantity(
                100.*ambi.press, unit=astropy.units.Pa),
            temperature=astropy.units.Quantity(
                ambi.temp, unit=astropy.units.deg_C),
            relative_humidity=ambi.rhum,
            obswl=astropy.units.Quantity(ambi.wlen, unit=astropy.units.um),
            az=astropy.coordinates.Angle(0,unit=astropy.units.degree),
            alt=astropy.coordinates.Angle(90,unit=astropy.units.degree))

        icrs_frame = astropy.coordinates.ICRS()
        # convert the zenith back to the ra/dec equatorial frame
        zeni= zeni_hori.transform_to(icrs_frame)
        # and the parallactic angle is the position angle of the zenith, coordinates.Angle
        pa = self.targ.position_angle(zeni)
        return pa
