# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-12-19
# @Filename: siderostat.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for siderostat field angles using homogeneous coordinates
"""

import math
import numpy
import datetime
import astropy.coordinates
import astropy.time
import astropy.units
import astropy.io.fits
import astropy.wcs

from lvmtipo.mirror import Mirror
from lvmtipo.target import Target
from lvmtipo.kmirror import Kmirror
from lvmtipo.fiber import Fiber
from lvmtipo.site import Site
from lvmtipo.ambient import Ambient
from lvmtipo.wcsarith import Wcsarith

__all__ = ['Siderostat']


class Siderostat():
    """ A siderostat of 2 mirrors
    ..todo.. add the roll angle as the 3rd angle besides azang and zenang
    """

    def __init__(self, zenang=90.0, azang=0.0, medSign=-1, m1m2dist=240.0, om1_off_ang = 118.969,
                 om2_off_ang = -169.752):
        """ A siderostat of 2 mirrors
        :param zenang: Zenith angle of the direction of the exit beam (degrees)
                   in the range 0..180. Default is the design value of the LVMT: horizontal
        :type zenang: float

        :param azang: Azimuth angle of the direction of the exit beam (degrees)
                   in the range -180..360 degrees, N=0, E=90.
                   Ought to be zero for the LCO LVMT where the FP is north of the
                   siderostat and 180 for the MPIA test setup where the FP is
                   south of the siderostat.. The default is the angle for LCO.
        :type azang: float

        :param medSign: Sign of the meridian flip design of the mechanics.
                       Must be either +1 or -1. Default is the LCO LVMT design as build (in newer
                       but not the older documentation).
        :type medSign: int

        :param m1m2dist: Distance between the centers of M1 and M2 in millimeter.
                       The default value is taken from
                       LVM-0098_Sky Baffle Design of 2022-04-18
                       by subracting the 84 and 60 cm distance of the
                       output pupils to M1 and M2.
        :type m1m2dist: float
        :param om2_off_ang: the offset angle which aligns PW motor angle of the M2 axis,
                       which is ax0 of the PW GUI, to the angles of the manuscript.
        :type om2_off_ang: float
        :param om1_off_ang: the offset angle which aligns PW motor angle of the M1 axis,
                       which is ax1 of the PW GUI, to the angles of the manuscript, degrees.
                      om1_off_ang and om2_off_ang are either angles in degrees or
                      both astropy.coordinates.Angle
        :type om1_off_ang: float
        """

        # the vector b[0..2] is the three cartesian coordinates
        # of the beam after leaving M2 in the topocentric horizontal system.
        # b[0] is the coordinate along E, b[1] along N and b[2] up.
        if isinstance(zenang, (int, float)) and isinstance(azang, (int, float)):
            self.b = numpy.zeros((3))
            self.b[0] = math.sin(math.radians(azang)) * \
                                 math.sin(math.radians(zenang))
            self.b[1] = math.cos(math.radians(azang)) * \
                                 math.sin(math.radians(zenang))
            self.b[2] = math.cos(math.radians(zenang))
        else:
            raise TypeError("invalid data types")

        self.m1m2len = m1m2dist

        if isinstance(medSign, int):
            if medSign in (1, -1) :
                self.sign = medSign
            else:
                raise ValueError("invalid medSign value")
        else:
            raise TypeError("invalid medSign data type")

        # axes orthogonal to beam. box points essentially to the zenith
        # and boy essentially to the East (at LCO) resp West (at MPIA)
        self.box = numpy.zeros((3))
        self.box[0] = 0.0
        self.box[1] = -self.b[2]
        self.box[2] = self.b[1]
        self.boy = numpy.cross(self.b, self.box)

        # The 3x3 B-matrix converts a (x,y,z) vector on the unit spehre
        # which is a direction from the observer to the star in alt-az (horizontal)
        # coordinates (x points East, y to the North and z to the zenith)
        # into a vector on the unit sphere where the two PW motor angles (and offsets)
        # play the role of the azimuth and polar angles.
        # See https://www.mpia.de/~mathar/public/mathar20201208.pdf
        self.B = numpy.array([[0,0,0],[0,0,0],[0,0,0]],dtype=numpy.double)
        bproj = math.hypot(self.b[1],self.b[2]) # sqrt(by^2+bz^2)
        self.B[0][0] = -self.sign*bproj
        self.B[0][1] = self.sign*self.b[0]*self.b[1]/bproj
        self.B[0][2] = self.sign*self.b[0]*self.b[2]/bproj
        self.B[1][0] = 0
        self.B[1][1] = -self.sign*self.b[2]/bproj
        self.B[1][2] = self.sign*self.b[1]/bproj
        self.B[2][0] = self.b[0]
        self.B[2][1] = self.b[1]
        self.B[2][2] = self.b[2]

        if isinstance(om1_off_ang, astropy.coordinates.Angle) :
            self.pw_ax_off = [om1_off_ang.radian, om2_off_ang.radian]
        else:
            self.pw_ax_off = [math.radians(om1_off_ang), math.radians(om2_off_ang)]

    def fieldAngle(self, site, target, ambi, time=None):
        """
        :param site: geographic ITRF location of the observatory
        :type site: lvmtipo.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates

        :param ambi: Ambient data relevant for refractive index
        :type ambi: lvmtipo.Ambient

        :param time: time of the observation /UTC; if None, the current time will be used.
        :type time:

        :return: field angle (direction to NCP) in radians
        :rtype: float
        """
        if isinstance(time, astropy.time.Time):
            now = time
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()

        # compute mirror positions
        horiz = target.toHoriz(site=site, ambi=ambi, time=now)
        # print(horiz)

        star = numpy.zeros((3))
        # same procedure as in the construction of b in the Sider ctor, but with 90-zenang=alt
        star[0] = math.sin(horiz.az.radian) * math.cos(horiz.alt.radian)
        star[1] = math.cos(horiz.az.radian) * math.cos(horiz.alt.radian)
        star[2] = math.sin(horiz.alt.radian)

        # unit vector from M2 to M1
        # we're not normalizing to self.m1m2len but keeping the vector
        # m2tom1 at length 1 to simplify the later subtractions to compute
        # normal vectors from other unit vector
        m2tom1 = numpy.cross(star, self.b)
        vlen = numpy.linalg.norm(m2tom1)
        m2tom1 /= self.sign*vlen

        # surface normal to M1 (not normalized to 1)
        m1norm = star - m2tom1
        # the orthogonal distance of the points of M1 to the origin
        # of coordinates are implied by the m1norm direction and
        # the fact that m2tom1 is on the surface. So the homogeneous
        # coordinate equation applied to m2tom1 should yield m2tom1 itself.
        # This requires (m2tom1 . m1norm -d) * n1morm=0 where dots are dot products.
        # the latter dot product actually requires a normalized m1norm
        vlen = numpy.linalg.norm(m1norm)
        m1norm /= vlen
        m1 = Mirror(m1norm, numpy.dot(m2tom1,m1norm))

        # surface normal to M2 (not normalized to 1)
        m2norm = self.b + m2tom1
        m2 = Mirror(m2norm, 0.0)

        # transformation matrix for the 2 reflections individually and in series
        m1trans = m1.toHomTrans()
        m2trans = m2.toHomTrans()
        trans = m2trans.multiply(m1trans)

        # for the field angle need a target that is just a little bit
        # more north (but not too little to avoid loss of precision)
	# 10 arcmin = 0.16 deg further to NCP
        targNcp = Target(target.targ.spherical_offsets_by(
                       astropy.coordinates.Angle("0deg"), astropy.coordinates.Angle("0.16deg")))

        # this here is the equivalent code to track a target with RA +0.16eg
        #    targNcp = Target(target.targ.spherical_offsets_by(
        #               astropy.coordinates.Angle("0.16deg"), astropy.coordinates.Angle("0deg")))
        # print("targNcp", to_north, targNcp.targ)
        # at this point using a positive value of the first paramter (0.16deg) gives
        # a target with larger RA value.
        # has
        horizNcp = targNcp.toHoriz(site=site, ambi=ambi, time=now)

        starNcp = numpy.zeros((3))
        # same procedure as in the construction of b in the Sider ctor,
        # but with 90 minus zenith angle=altitude
        starNcp[0] = math.sin(horizNcp.az.radian) * \
                              math.cos(horizNcp.alt.radian)
        starNcp[1] = math.cos(horizNcp.az.radian) * \
                              math.cos(horizNcp.alt.radian)
        starNcp[2] = math.sin(horizNcp.alt.radian)

        # image of targNcp while hitting M1
        m1img = trans.apply(m2tom1)

        # image of targNcp before arriving at M1
        star_off_m1 = m2tom1 + starNcp

        starimg = trans.apply(star_off_m1)

        # virtual direction of ray as seen from point after M2
        # no need to normalize this because the atan2 will do...
        # sign was wrong until 2022-11-19: we need to take the direction
        # from the on-axis star (m1img) to the off-axis start (starimg).
        starvirt = starimg - m1img

        # project in a plane orthogonal to  self.b
        cos_fang = numpy.dot(starvirt, self.box)
        sin_fang = numpy.dot(starvirt, self.boy)
        return math.atan2(sin_fang, cos_fang)

    def to_header(self, site, target, camera, ambi, k_mocon_steps, ag_cam_name,
           genrevx, genrevy, kmirr, time=None, pixsize=None, bin=None, wd=None, hd=None,
           flen=1839.8, dist_cam_edge=11.14771) :
        """ Convert the parameters to FITS WCS header cards.
        This traces the position angle implied by the siderostat and number of mirror
        reflections through the K-mirror (if applicable) and through the
        prisms if applicable to predict the key WCS keywords.
        .. warn:: this assumes that the FITS files are read out as
           in the MPIA test setup: all east/west/center cameras with the
           long edge up-down, no rotation by any 90 or 180 degrees in any
           python software that may spring up after Oct 2022 and no ad-hoc
           decisions to mount the cameras in different orientations.

        :param site: geographic ITRF location of the observatory
        :type site: lvmtipo.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates

        :param camera: the scales in the image, pixel and binning etc
            Can be None, but then the pixsize, bin, wd, hd and flen parameters must be set
        :type camera: SkymakerCamera

        :param ambi: Ambient data relevant for refractive index
        .. warn:: The pixel scales in the center camera are wrong as of 2022-11-04
                  because the configuration yml files do not contain the correct 4.5 um
        :type ambi: lvmtipo.Ambient

        :param k_mocon_steps: The location of the K-mirror at the time
             in units of Mocon steps. We avoid the use of bare angles because
             various different conventions are in use for this instrument,
             both with respect to sign and with respect to mechanical vs
             optical angles.
             If kmirr=None this parameter is not used.
        :type k_mocon_steps: int

        :param ag_cam_name: The telescope (skye, skyw, spec, sci) and camera
             location (age, agw, agc) or (east, west, center)
             of the FLIR camera that took the image. The string must be some
             concatenation of both names.
        :type ag_cam_name: string

        :param genrevx: True if the readout direction along x was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevx: bool

        :param genrevy: True if the readout direction along y was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevy: bool

        :param time: time of the observation /UTC; if None, the current time will be used.
        :type time:

        :param pixsize: pixel size in microns, not binned
              Can be None if camera has this information
        :type pixsize: float

        :param bin: binning factor in the image
              Can be None if camera has this information
        :type bin: int

        :param wd: width of the image in pixels (after binning)
              Can be None if camera has this information
              This is typically 1600 for the monochrome cameras and 3200 for
              the color camera (binning=1).
        :type wd: int

        :param hd: height of the image in pixels (after binning)
              Can be None if camera has this information
              This is typically 1100 for the monochrome cameras and 2200 for
              the color camera (binning=1).
        :type hd: int

        :param flen: focal lens of lens in mm.
              Can be None if camera has this information
        :type flen: float

        :param kmirr: K-mirror constants related to conversions of steps to angles
             May be None to indicate that there is no K-mirror (like for the spec telescope)
        :type kmirr: lvmtipo.Kmirror

        :param dist_cam_edge: the distance of the long east or west FLIR camera edge to
             the fiber bundle centre in millimeters
             Is 7.144+4.0 mm according to SDSS-V_0110 figure 6
             and 11.14471 according to figure 3-1 of LVMi-0081
        :type dist_cam_edge: float

        :return a set of FITS header cards with WCS keywords.
        :rtype astropy.io.fits.Header
        """

        # start with header cards of geographic coords of the observatory
        wcshdr = site.to_header()

        # compute the position angle (NCP) after the siderostat,
        # before the K-mirror, in radians.
        if isinstance(time, astropy.time.Time):
            now = time
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()
        ang = self.fieldAngle(site, target, ambi, time=now)
        # angtmp = self.fieldAngle(site, target, ambi, time=now,to_north=False)

        # at this point (after M2) the angular value for E targets
        # is 90 deg larger, ie. they are right from the targets for N targets
        # so the image is flipped, whch reflects that a lens flips images
        # not actually the lens but the act of turning around to look at the
        # picture on some downstream detector.
        # print("now",now,"ang N",math.degrees(ang),"deg")

        # number of mirror reflections so far, 2 of the siderostat
        n_refl = 2

        # Simplify the interface by also recognizing the descr of
        # the camera (which uses upper case for some letter).
        ag_cam_name = ag_cam_name.lower()

        # three LVMT benches have a K-mirror, one has not
        # Calculate the position angle after (any) K-mirror
        if kmirr is None:
            # this will be the case of the spec-telescope at LCO
            # This may be wrong for some MPIA setups
            # the angle passes unchanged
            pass
        else:
            # this will eventually be the case of the sci, skye, skyw telescopes at LCO
            # This may be wrong for some MPIA setups
            n_refl += 3

            # convert K-mirror orientation to radians
            kangle = kmirr.steps_to_radians(k_mocon_steps)
            # print("kmirr ",math.degrees(kangle))

            # flip the component of the position angle in the incidence
            # plane(s) of the K-mirror. So kangle+-90 deg is the
            # arithmetic mean of the position angle before and after the K-mirror.
            # (ang_in+ang_out)/2 = kangle+-pi/2
            # ang_out = 2*kangle+-pi - ang_in
            ang = math.pi + 2.*kangle - ang
            # angtmp = math.pi + 2.*kangle - angtmp

        # at this point (after K-mirror) the angular value for E targets
        # is 90 deg smaller, ie. they are left from the targets for N targets.
        # print("after kmirr",math.degrees(ang)," N deg")
        # print("after kmirr",math.degrees(angtmp)," E deg")
        tele='LVMT'
        if ag_cam_name.find('skye') >= 0 :
            tele +=' skye'
        elif ag_cam_name.find('skyw') >= 0 :
            tele +=' skyw'
        elif ag_cam_name.find('sci') >= 0 :
            tele +=' sci'
        elif ag_cam_name.find('spec') >= 0 :
            tele +=' spec'
        key = astropy.io.fits.Card("TELESCOP", tele, "Local Volume Mapper sci, skye, skyw, spec")
        wcshdr.append(key)

        # print(tele)
        # patch error in camera names (S/N assignments )before 2022-11-20.
        # This block of code should be removed to increase efficiency at some time in the future
        patchTime2 = astropy.time.Time('2022-11-20T00:00:00',format='isot',scale='tai')
        # positive if patched/corrected, i.e., after patchTime2
        sinceTime2 = now - patchTime2
        if sinceTime2.to_value('jd') < 0:
            # in the early configuration files the data streams of the east
            # and west cameras at MPIA were wrong/swapped
            # This occurs if the IP addresses (in configuration and static dhcp)
            # are wrong so the 16-bit-byte-streams of the east camera appear here
            # as the west camera and vice versa.
            swap_ew_cams = True
        else:
            swap_ew_cams = False

        is_agw = is_age = is_agc = False
        if ag_cam_name.find('west') >= 0 or ag_cam_name.find('agw') >= 0:
            if swap_ew_cams:
                is_age=True
            else :
                is_agw=True
        elif ag_cam_name.find('east') >= 0 or ag_cam_name.find('age') >= 0:
            if swap_ew_cams:
                is_agw=True
            else :
                is_age=True
        elif ag_cam_name.find('center') >= 0 or ag_cam_name.find('agc') >= 0:
            # this is last in the elif-chain because the name may have the substring 'agcam...'
            is_agc=True
        else:
            raise NameError("Unrecognized camera name " + ag_cam_name)

        # since the camera has already performed internally one
        # flip around the long axis assuming there is a single lens
        # and that the output stream should be in scan raster format we start
        # with flipCDy=true, which means the camera reads the pixels (assuming
        # the PoE connector is at the bottom) from bottom row to top row, left to right.
        # Note that almost always digrevy=true later on, such that we actully
        # may have read this top to bottom. flipCDy is relative to the scan raster format.
        flipCDy=True
        flipCDx=False

        # next recalculate the position angle passing by a prism
        # (if there is any) and recalculating it in the
        # camera orientations. Th ecamera orientations are all
        # with the long edge vertical. The new angles are those
        # where x is along the long edge and y is positive if pointing
        # away from the long side with the PoE connector.
        if is_agw :
            # prism incidence plane is horizontal. Flip that component
            ang *= -1.0
            n_refl += 1
            # assume camera is mounted vertically with the long lower edge
            # pointing south, away from the focal plane
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "right" or geographicall "to the north"
            # (passive rotation)
            ang -= math.pi/2.0
        elif is_age :
            # prism incidence plane is horizontal. Flip that component
            ang *= -1.0
            n_refl += 1
            # assume camera is mounted vertically with the long lower edge
            # pointing south (basically rotated 180 deg around the horizontal
            # axis relative to the west camera) away from the FP.
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "right" or geographically "to the south"
            # (passive rotation)
            ang += math.pi/2.0
        else:
            # center camera no prism
            # assume camera is mounted vertically with the long lower edge
            # pointing east ( basically rotating the east camera 90 deg
            # around the vertical axis and then 180 deg around the horizontal axis.
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "left" or "to the west"
            # (passive rotation)
            ang += math.pi/2.0

        if is_agc :
            tele += '.agc'
        elif is_agw :
            tele += '.agw'
        elif is_age :
            tele += '.age'
        else:
            tele += '.'
        key = astropy.io.fits.Card("CAMERA", tele, " AG camera age, agw, agc on the bench")
        wcshdr.append(key)
        # print(tele,"on chip",math.degrees(ang),"deg")

        # up to here ang is refering to the short-edge-up orientation
        # of the camera if that were read out by streaming the araviscam
        # bytes in standard order (genrevx=gerevy=false). Because we need
        # one flip up-down (first row becomes last row) to flush the stream
        # into the standard FITS order where the lowest row is the first to be
        # read, the optically relevant value for that flip is opposite to the
        # nominal value used with the araviscam interface

        # we're left with 4 combinations of genrev[xy] individually true or false
        if genrevx:
            if genrevy:
                # both flipped (i.e., a 180 deg rotation)
                # ang += math.pi
                flipCDx = not flipCDx
                flipCDy = not flipCDy
            else:
                # x flipped, y not flipped
                # ang = -ang
                flipCDx = not flipCDx
        else:
            if genrevy:
                # x not flipped but y flipped
                # Note that this is the operation as above for the kmirror at angle 0
                # a y-flip is the same as an x-flip plus a 180 deg rotation
                # ang = math.pi - ang
                flipCDy = not flipCDy
        # print(tele,"with flips",math.degrees(ang),"deg")
        # print(tele,"flipCDx",flipCDx,"flipCDy",flipCDy)

        # add the main pointing keywords
        key = astropy.io.fits.Card("CUNIT1", "deg", "WCS units along axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CUNIT2", "deg", "WCS units along axis 2")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE1", "RA---TAN", "WCS type axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE2", "DEC--TAN", "WCS type axis 2")
        wcshdr.append(key)

        hexrep = target.targ.ra.to_string(unit=astropy.units.hour,precision=2)
        key = astropy.io.fits.Card(
            "CRVAL1", target.targ.ra.degree, "[deg] RA ref pixel " + hexrep)
        wcshdr.append(key)
        hexrep = target.targ.dec.to_string(unit=astropy.units.degree,precision=2)
        key = astropy.io.fits.Card(
            "CRVAL2", target.targ.dec.degree, "[deg] Dec ref pixel " + hexrep)
        wcshdr.append(key)

        # convert disgtcamEdgeCtr from mm to micron and then to pixels
        # This is the distance from the long edge that is narrowest to
        # the fiber center after projection into the focal plane.
        if camera is None:
            dist_cam_edge *= 1000.0 / (pixsize*bin)
        else:
            dist_cam_edge *= 1000.0 / (camera.pixsize*camera.binning[1])

        # distance from the middle of the east or west camera to fiber ctner in pixels
        if camera is None:
            # note that hd is already in binned pixels (as in the FITS NAXIS)
            dist_cam_mid = hd / 2 + dist_cam_edge
        else:
            dist_cam_mid = camera.detector_size.hd / \
                2 / camera.binning[1] + dist_cam_edge

        # where is the center of image away in the pixel coordinate system
        # of the camera? For the x-position this is in the middle of the
        # camera along the long axis because all cameras are installed up-down
        # For age and agw the wd parameters are 1600 and the hd 1100 (not binned).
        if camera is None:
            crpix1 = wd / 2 + 0.5
            crpix2 = hd / 2 + 0.5
        else :
            crpix1 = camera.detector_size.wd / 2 / camera.binning[0] + 0.5
            crpix2 = camera.detector_size.hd / 2 / camera.binning[1] + 0.5

        if is_agc :
            # ra/dec in the center of the image
            pass
        elif is_agw :
            # the direction from the prism center to the fiber center
            # is walkign east, which is in the camera walking to the lower edge
            # in raster scan order (flipCDy=false) this means increasing pixy
            if flipCDy:
                crpix2 -= dist_cam_mid
            else:
                crpix2 += dist_cam_mid
        elif is_age :
            # the direction from the prism center to the fiber center
            # is walkign west, which is in the camera walking to the lower edge
            # in raster scan order (flipCDy=false) this means increasing pixy
            if flipCDy:
                crpix2 -= dist_cam_mid
            else:
                crpix2 += dist_cam_mid

        key=astropy.io.fits.Card( "CRPIX1", crpix1, "[px] point cntr along axis 1")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CRPIX2", crpix2, "[px] point cntr along axis 2")
        wcshdr.append(key)

        # the two main values in the CD matrix, irrespective of sign
        if camera is None:
            # plate scale, microns in the focal plane per radian on the sky
            # image scale at focus is 8.92 um/arcsec = 8.92um/(pi/180/3600 rad) =1.8e^6 um/rad
            pscal=1000.*flen
            # plate scale, pixels in the focal plane per radian on the sky
            pscal /= pixsize
            # plate scale, pixels in the focal plane per degrees on the sky
            pscal = math.radians(pscal)
            cosperpix = math.cos(ang)/pscal
            sinperpix = math.sin(ang)/pscal
        else:
            cosperpix = camera.degperpix*math.cos(ang)
            sinperpix = camera.degperpix*math.sin(ang)
        # suppose the parity is positive and ang is small, then  in matrix notation
        # ( x)     ( sin ang, -cos ang ) (delta)
        #       =
        # ( y)     ( cos ang, sin ang  ) (alpha)
        # and the CD is the inverse
        # ( delta)     ( sin ang, cos ang ) (x)
        #       =
        # ( alpha)     ( -cos ang, sin ang  ) (y)
        # but with alpha the first coordinate:
        # ( alpha)     ( -cos ang, sin ang ) (x)
        #       =
        # ( delta)     ( sin ang, cos ang  ) (y)

        xy_to_ad=[[-cosperpix,sinperpix],[sinperpix,cosperpix]]
        if flipCDx :
            # direction of x pixels flipped
            xy_to_ad[0][0] *= -1
            xy_to_ad[1][0] *= -1
        if flipCDy :
            # direction of y pixels flipped
            xy_to_ad[0][1] *= -1
            xy_to_ad[1][1] *= -1

        # negative CD determinant for upright images (imgparity=True), because
        # (ra/dec) are a left-handed coordinate system; positive for flipped (imparity=False)

        # Assign an image parity defined by the even/odd number of mirror reflections
        # this is the parity one expects in the FITS image, irrespective
        # of whether readout serialization has changed directions with genrev[xy].
        imgparity = True if ( (n_refl % 2) == 0 ) else False
        # ... and there is the lens which also flips once, but the FLIR has
        # already taken this into account by permuting the pixels in the byte stream.
        # imgparity = True if ( (n_refl % 2) != 0 ) else False

        if not imgparity:
            # direction of alpha (i.e. the image) is flipped
            xy_to_ad[0][0] *= -1
            xy_to_ad[0][1] *= -1

        key=astropy.io.fits.Card( "CD1_1", xy_to_ad[0][0], "[deg/px] WCS matrix diagonal x->ra")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD1_2", xy_to_ad[0][1], "[deg/px] WCS matrix outer diagonal y->ra")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_1", xy_to_ad[1][0], "[deg/px] WCS matrix outer diagonal x->delta")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_2", xy_to_ad[1][1], "[deg/px] WCS matrix diagonal y-> delta")
        wcshdr.append(key)

        # convert the equatorial to horizontal coordinates
        horiz = target.toHoriz(site=site, ambi=ambi, time=now)
        key=astropy.io.fits.Card(
            "ALT", horiz.alt.deg, "[deg] object altitude above horizon")
        wcshdr.append(key)
        key=astropy.io.fits.Card(
            "AZ", horiz.az.deg, "[deg] object azimuth N=0, E=90")
        wcshdr.append(key)

        # compute the parallactic angle
        pa = target.parallact(site, ambi, now)
        key=astropy.io.fits.Card("PARANG", pa.deg, "[deg] parallactic angle")
        wcshdr.append(key)

        # append another WCS 'F' which maps pixels to microns in the focal plane
        self.to_header_wcs_f(wcshdr,is_agw, is_age, genrevx, genrevy,
            pixsize, wd, hd, bin=bin, dist_cam_edge=dist_cam_edge)

        # append another WCS 'A' which maps pixels to the alt-az coordinates
        self.to_header_wcs_a(wcshdr, numpy.array(xy_to_ad,numpy.double), pa, horiz, crpix1,crpix2)

        # update DATE
        fwritetime = datetime.datetime.utcnow()
        key=astropy.io.fits.Card(
            "DATE", fwritetime.strftime("%Y-%m-%dT%H:%M:%S"), "UTC " + str(type(self)) + " keys modification")
        wcshdr.append(key)

        return wcshdr

    def to_header_wcs_a(self, wcshdr, xy_to_ad, pa, horiz,crpix1,crpix2) :
        """ Collect the parameters for a WCS system for alt/az angles.
        The strategy is to mix the parallactic angle into the CDi_j matrix
        of the primary ra/dec WCS such that aspects of flips etc are already known from this.

        :param wcshdr: the list of header cards to be augmented here
        :type wcshdr: astropy.io.fits.Header

        :param xy_to_ad: The CDi_j rotation for the converstion to ra/dec angles
        :type xy_to_ad: numpy.array

        :param pa: The parallactic angle
        :type pa: astropy.coordinates.Angle

        :param horiz: the alt/az of the pointing center
        :type horiz: astropy.coordinates.AltAz

        :param crpix1: pixel along x (FITS coordnates) of the horiz point
        :type crpix1: float

        :param crpix1: pixel along y (FITS coordnates) of the horiz point
        :type crpix1: float
        """

        # select suffix 'A' for this WCS, valid values are 'A'to 'Z' and
        # this must not be the same as in the to_header_wcs_f() add-on
        suff='A'

        # add the main pointing keywords
        key = astropy.io.fits.Card("CUNIT1"+suff, "deg", "WCS units along axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CUNIT2"+suff, "deg", "WCS units along axis 2")
        wcshdr.append(key)
        # key = astropy.io.fits.Card("CTYPE1"+suff, "AZ---TAN", "WCS type axis 1")
        key = astropy.io.fits.Card("CTYPE1"+suff, "HOLN-TAN", "WCS type axis 1")
        wcshdr.append(key)
        # key = astropy.io.fits.Card("CTYPE2"+suff, "ALT--TAN", "WCS type axis 2")
        key = astropy.io.fits.Card("CTYPE2"+suff, "HOLT--TAN", "WCS type axis 2")
        wcshdr.append(key)

        key=astropy.io.fits.Card( "CRPIX1"+suff, crpix1, "[px] point cntr along axis 1")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CRPIX2"+suff, crpix2, "[px] point cntr along axis 2")
        wcshdr.append(key)

        # put the origin of the focal plane coordinate system at the center of
        # the center camera of fiber bundle
        key = astropy.io.fits.Card( "CRVAL1"+suff, horiz.az.deg, "[deg] azimuth at ref pixel")
        wcshdr.append(key)
        key = astropy.io.fits.Card( "CRVAL2"+suff, horiz.alt.deg, "[deg] altitude at ref pixel")
        wcshdr.append(key)

        # if xy_to_ad is the existing transformation to Ra/dec we need
        # to find the (left) matrix ad_to_aa depending on pa which rotations
        # this further. If pa is positive and the images are flipped or not flipped
        # (az)    (-cos pa sin pa)   (ra)
        #       =                  *
        # (alt)   ( sin pa cos pa)  (dec)
        cospa = math.cos(pa.radian)
        sinpa = math.sin(pa.radian)
        ad_to_aa = numpy.array([[-cospa,sinpa],[sinpa,cospa]],numpy.double)
        xy_to_aa = numpy.matmul(ad_to_aa, xy_to_ad)

        key=astropy.io.fits.Card( "CD1_1"+suff, xy_to_aa[0][0], "[deg/px] WCS matrix diagonal x->az")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD1_2"+suff, xy_to_aa[0][1], "[deg/px] WCS matrix outer diagonal y->az")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_1"+suff, xy_to_aa[1][0], "[deg/px] WCS matrix outer diagonal x->alt")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_2"+suff, xy_to_aa[1][1], "[deg/px] WCS matrix diagonal y-> alt")
        wcshdr.append(key)

        #key=astropy.io.fits.Card( "LONPOLE"+suff, 50., "[deg]")
        #wcshdr.append(key)
        #key=astropy.io.fits.Card( "LATPOLE"+suff, 30., "[deg]")
        #wcshdr.append(key)

    def to_header_wcs_f(self, wcshdr, is_agw, is_age,
           genrevx, genrevy, pixsize, wd, hd, bin =1,
           dist_cam_edge=11.14771) :
        """ Collect the parameters for a second WCS system for focal plane coordinates.
        The primary list of parameters collected in to_header_wcs() convert pixels
        to RA/DEC coordinates, whereas this set here converts pixels to
        micrometers in a comman focal plane of the guide cameras and fiber heads.
        The names for this system are ending on 'F' (short for focal plane).
        The main purpose is to have a quick online conversion from guide camera
        pixels to fiber head positions.

        :param wcshdr: the list of FITS header cards to be augmented here.
        :type wcshdr: SkymakerCamera

        :param is_agw: true if this describes pixels on a west AG camera
        :type is_agw: bool

        :param is_age: true if this describes pixels on a east AG camera
             If neither is_agw nor is_age are True we consider this a center camera.
        :type is_age: bool

        :param genrevx: True if the readout direction along x was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevx: bool

        :param genrevy: True if the readout direction along y was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevy: bool

        :param pixsize: pixel size in microns, not binned
        :type pixsize: float

        :param bin: binning factor in the image
        :type bin: int

        :param wd: width of the image in pixels (after binning)
              This is typically 1600 for the monochrome cameras and 3200 for
              the color camera (binning=1).
        :type wd: int

        :param hd: height of the image in pixels (after binning)
              This is typically 1100 for the monochrome cameras and 2200 for
              the color camera (binning=1).
        :type hd: int

        :param dist_cam_edge: the distance of the long east or west FLIR camera edge to
             the fiber bundle centre in millimeters
             Is 7.144+4.0 mm according to SDSS-V_0110 figure 6
             and 11.14471 according to figure 3-1 of LVMi-0081
        :type dist_cam_edge: float
        """

        # suffix for the WCS keywords. Must be on of the letters A to Z
        suff = 'F'

        if is_agw or is_age:
            is_agc = False
        else:
            is_agc = True

        # since the camera has already performed internally one
        # flip around the long axis assuming there is a single lens
        # and that the output stream should be in scan raster format we start
        # with flipCDy=true, which means the camera reads the pixels (assuming
        # the PoE connector is at the bottom) from bottom row to top row, left to right.
        # Note that almost always digrevy=true later on, such that we actully
        # may have read this top to bottom. flipCDy is relative to the scan raster format.
        flipCDy=True
        flipCDx=False

        # 4 combinations of genrev[xy] individually true or false
        if genrevx:
            if genrevy:
                # both flipped (i.e., a 180 deg rotation)
                flipCDx = not flipCDx
                flipCDy = not flipCDy
            else:
                # x flipped, y not flipped
                flipCDx = not flipCDx
        else:
            if genrevy:
                # x not flipped but y flipped
                # Note that this is the operation as above for the kmirror at angle 0
                # a y-flip is the same as an x-flip plus a 180 deg rotation
                flipCDy = not flipCDy

        #if is_age or is_agw:
        #    n_refl = 1
        #elif is_agc :
        #    n_refl = 0
        #else :
        #    raise ValueError("one of age, age, agc must be True ")

        # add the main pointing keywords
        key = astropy.io.fits.Card("CUNIT1"+suff, "um", "WCS units along axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CUNIT2"+suff, "um", "WCS units along axis 2")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE1"+suff, "FOCP-X", "WCS type axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE2"+suff, "FOCP-Y", "WCS type axis 2")
        wcshdr.append(key)

        # put the origin of the focal plane coordinate system at the center of
        # the center camera of fiber bundle
        key = astropy.io.fits.Card( "CRVAL1"+suff, 0, "[um] FP x for ref pixel")
        wcshdr.append(key)
        key = astropy.io.fits.Card( "CRVAL2"+suff, 0, "[um] FP y for ref pixel")
        wcshdr.append(key)

        # relevant pixel size for our conversion are the binned ones
        pixsize *= bin

        # convert disgtcamEdgeCtr from mm to micron and then to pixels
        # This is the distance from the long edge that is narrowest to
        # the fiber center after projection into the focal plane.
        dist_cam_edge *= 1000.0 / pixsize

        # distance from the middle of the east or west camera to fiber ctner in pixels
        # note that hd is already in binned pixels (as in the FITS NAXIS)
        dist_cam_mid = hd / 2 + dist_cam_edge

        # where is the center of image away in the pixel coordinate system
        # of the camera? For the x-position this is in the middle of the
        # camera along the long axis because all cameras are installed up-down
        # For age and agw the wd parameters are 1600 and the hd 1100 (not binned).
        crpix1 = wd / 2 + 0.5
        crpix2 = hd / 2 + 0.5

        if is_agc :
            # x pixel points up and y-pixel points left
            xy_to_ad=[[0,-pixsize],[pixsize,0]]
        elif is_agw :
            # the direction from the prism center to the fiber center
            # is walkign east, which is in the camera walking to the lower edge
            # in raster scan order (flipCDy=false) this means increasing pixy
            if flipCDy:
                crpix2 -= dist_cam_mid
            else:
                crpix2 += dist_cam_mid
            # x pixel points down and y-pixel points into focal plane and
            # considering the prism reflection left in the focal plane
            xy_to_ad=[[0,-pixsize],[-pixsize,0]]
        elif is_age :
            # the direction from the prism center to the fiber center
            # is walkign west, which is in the camera walking to the lower edge
            # in raster scan order (flipCDy=false) this means increasing pixy
            if flipCDy:
                crpix2 -= dist_cam_mid
            else:
                crpix2 += dist_cam_mid
            # x pixel points up and y-pixel points into focal plane and
            # considering the prism reflection right in the focal plane
            xy_to_ad=[[0,pixsize],[pixsize,0]]
        else :
            raise ValueError("one of age, age, agc must be True ")

        key=astropy.io.fits.Card( "CRPIX1"+suff, crpix1, "[px] point cntr along axis 1")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CRPIX2"+suff, crpix2, "[px] point cntr along axis 2")
        wcshdr.append(key)

        if flipCDx :
            # direction of x pixels flipped
            xy_to_ad[0][0] *= -1
            xy_to_ad[1][0] *= -1
        if flipCDy :
            # direction of y pixels flipped
            xy_to_ad[0][1] *= -1
            xy_to_ad[1][1] *= -1

        key=astropy.io.fits.Card( "CD1_1"+suff, xy_to_ad[0][0], "[um/px] WCS matrix diagonal x->umHor")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD1_2"+suff, xy_to_ad[0][1], "[um/px] WCS matrix outer diagonal y->umHor")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_1"+suff, xy_to_ad[1][0], "[um/px] WCS matrix outer diagonal x->umVert")
        wcshdr.append(key)
        key=astropy.io.fits.Card( "CD2_2"+suff, xy_to_ad[1][1], "[um/px] WCS matrix diagonal y->umVert")
        wcshdr.append(key)

    def to_wcs(self, site, target, camera, ambi, k_mocon_steps, ag_cam_name,
           genrevx, genrevy, kmirr, time=None,pixsize=None,bin=None,wd=None,hd=None,flen=1839.8,
           dist_cam_edge=11.14771):
        """ Convert the parameters to astropy WCS information.
        This is a simple convenience wrapper around to_header().

        :param site: geographic ITRF location of the observatory
        :type site: lvmtipo.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates

        :param camera: the scales in the image, pixel and binning etc
        :type camera: SkymakerCamera

        :param ambi: Ambient data relevant for refractive index
        .. warn:: The pixel scales in the center camera are wrong as of 2022-11-04
                  because the configuration yml files do not contain the correct 4.5 um
        :type ambi: lvmtipo.Ambient

        :param k_mocon_steps: The location of the K-mirror at the time
             in units of Mocon steps. We avoid the use of bare angles because
             various different conventions are in use for this instrument,
             both with respect to sign and with respect to mechanical vs
             optical angles.
             If kmirr=None this argument is not used.
        :type k_mocon_steps: int

        :param ag_cam_name: The telescope (skye, skyw, spec, sci) and camera
             location (age, agw, agc) or (east, west, center)
             of the FLIR camera that took the image. The string must be some
             concatenation of both names.
        :type ag_cam_name: string

        :param genrevx: True if the readout direction along x was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevx: bool

        :param genrevy: True if the readout direction along y was reversed
             in the araviscam interface before creating the FITS file.
        :type genrevy: bool

        :param time: time of the observation /UTC; if None, the current time will be used.
        :type time:

        :param pixsize: pixel size in microns, not binned
              Can be None if camera has this information
        :type pixsize: float

        :param bin: binning factor in the image
              Can be None if camera has this information
        :type bin: int

        :param wd: width of the image in pixels (after binning)
              Can be None if camera has this information
        :type wd: int

        :param hd: height of the image in pixels (after binning)
              Can be None if camera has this information
        :type hd: int

        :param flen: focal lens of lens in mm.
              Can be None if camera has this information
        :type flen: float

        :param kmirr: K-mirror constants related to conversions of steps to angles
             May be None to indicate that there is no K-mirror (like for the spec telescope)
        :type kmirr: lvmtipo.Kmirror

        :param dist_cam_edge: the distance of the long east or west FLIR camera edge to
             the fiber bundle centre in millimeters
        :type dist_cam_edge: float

        :return: a set of astropy WCS keywords.
        :rtype: astropy.wcs.WCS
        """

        hdr= self.to_header(site, target, camera, ambi, k_mocon_steps, ag_cam_name,
           genrevx, genrevy, kmirr, time=time,pixsize=pixsize,bin=bin,wd=wd,hd=hd,flen=flen,
           dist_cam_edge=dist_cam_edge )

        return astropy.wcs.WCS(hdr)

    def update_fits(self, fits_in, fits_out):
        """
        Update the WCS header cards of an existing FITS file with
        computed a-priori parameters  found in the header.
        The header parameters OBSERVAT, DATE-OBS, GENREV[XY], BINX,
            KMIRDROT, CAMNAME etc are fed into the computed position angles for
            the fits_in image (in the primary header), and fits_out
            is created with the same image but WCS keywords replaced according to the
            theory of the LVMT mirror trains.
        This is essentially a checker for existing FITS files so the result can
          be compared with astrometry.net solutions, and a demonstrator of
          how first esimators to run astrometry.net  could be found.

        ..warn..: this works only assuming that some standard conventions in the
            header data (as applying in Nov 2022) are met.
        :param fits_in: existing fits file name on local disks
        :type fits_in: str

        :param fits_out: file name of a non-existing FITS file on local disks to be created
            To prevent creation of the FITS file this can be None, so effectively
            the returned WCS structure is computed.
        :type fits_out: str

        :return: a set of astropy WCS keywords computed from the other header keywords.
        :rtype: astropy.wcs.WCS

        ..todo..: should probably put in the __main__ of the fieldrotation sample

        Usage example: (first in bash rm tst[WCE].fits in the shell, then thisprog.py, then
          ds9 -mosaic tst[WCE].fits)
        #!/usr/bin/env python3

        from lvmtipo.siderostat import Siderostat
        import astropy.io.fits

        fits_in="/home/mathar/lvm.sci.agcam.west_00001112.fits"
        fits_out="tstW.fits"
        sid=Siderostat()
        sid.update_fits(fits_in,fits_out)

        fits_in="/home/mathar/lvm.sci.agcam.east_00001112.fits"
        fits_out="tstE.fits"
        sid=Siderostat()
        sid.update_fits(fits_in,fits_out)

        fits_in="/home/mathar/lvm.spec.agcam.center_00000731.fits"
        fits_out="tstC.fits"
        sid=Siderostat()
        sid.update_fits(fits_in,fits_out)
        """

        # grab header of PHDU of fits_in
        hdu = astropy.io.fits.open(fits_in)
        hdr = hdu[0].header
        # extract all the wcs-relevant data

        # this is MPIA or LCO
        observat = hdr['OBSERVAT'].strip()
        site = Site(name=observat)

        # ra and dec in degrees
        ra = hdr['RA']
        dec = hdr['DEC']
        # hack 1: correct the absolutely wrong values of 2022-10-26 for Alcyrone
        # ra=(3+47.0/60.0+29.19/3600.)*15 ;
        # dec= 24+6/60.0+14.8/3600. ;

        # horizontal and vertical avariscam flips in the image
        genrevx = hdr['GENREVX']
        genrevy = hdr['GENREVY']

        # K-mirror (middle mirror) position angle in degrees
        # The KMIRDROT uses a convnetion where this is apparently
        # the mechanical angle with a sign such that -135 indicates the home position.
        kmang = -hdr['KMIRDROT']

        # DATE-OBS is in TAI time scales (37 seconds off UTC)
        time = astropy.time.Time(hdr['DATE-OBS'],format='isot',scale='utc')

        # patch error in FITS headers before 2022-10-30.
        # This block of code should be removed to increase efficiency at some time in the future
        patchTime1 = astropy.time.Time('2022-10-30T00:00:00',format='isot',scale='tai')
        # positive if patched/corrected, i.e., after patchTime2
        sinceTime1 = time - patchTime1
        if sinceTime1.to_value('jd') < 0:
            genrevy = True

        # binning factor 1,2,4
        bin = hdr['BINX']
        # FLIR camera type (to deduce the pixel size in microns)
        flirtyp = hdr['CAMTYPE']
        ag_cam_name = hdr['CAMNAME'].strip()

        # some less important data to build a refractivity model
        celsius = hdr['BENTEMP']
        relhum = hdr['BENHUM']
        pres = hdr['BENPRESS']

        # binned size of the image
        wd=hdr['NAXIS1']
        hd=hdr['NAXIS2']

        # try to guess a telescope bench from the file name
        tele=''
        if 'sci' in fits_in :
            tele += 'sci'
        elif 'skyw' in fits_in :
            tele += 'skyw'
        elif 'skye' in fits_in :
            tele += 'skye'
        elif 'spec' in fits_in :
            tele += 'spec'

        # consistency check: if the telescope is not spec this
        # value needs to be <= 137 which is a rough estimate for the home switch offset
        # For the spec telescope there is no k-mirror and this angle does not matter
        if 'spec' not in tele:
            if abs(kmang) > 140:
                # for the sci, skye and skyw we need an angle
                # but this seems not to be correct here...
                raise ValueError("unsupported K-mirror angle " + str(kmang) + " on " + tele)
            else:
                kmirr = Kmirror()
                k_mocon_steps=kmirr.radians_to_steps(math.radians(kmang))
        else:
            kmirr = None
            k_mocon_steps=None


        if relhum >= 0 and pres >0 and celsius > -20 :
            ambi=Ambient(press=pres, temp=celsius, rhum=relhum)
        else:
            ambi=Ambient()
        obj = astropy.coordinates.SkyCoord(ra=ra, dec=dec,unit="deg")
        target=Target(obj)
        # print(target)

        # Because the arcsecs per pixel kw in the FITS files is usually wrong
        # for the center camera, we derive our own patched values here....
        if '70S7C' in flirtyp :
            pixsize=4.5
        else:
            pixsize=9.0

        wcshdr = self.to_header(site, target, None, ambi, k_mocon_steps, tele + " " + ag_cam_name,
                genrevx, genrevy, kmirr, time, pixsize, bin=bin, wd=wd, hd=hd)
        # print(wcshdr)

        if fits_out is not None :
            # In this case we want to create a new not yet existing FITS file
            # copy the data/image to a new HDU
            hduout = astropy.io.fits.PrimaryHDU(hdu[0].data)

            # copy also the header.
            # Need copy() here, otherwise the update() further down will refuse to replace
            # the old CD values which are in the old FITS file.
            hduout.header = hdu[0].header.copy()
            hduout.header.update(wcshdr,unique=True)
            # todo: remove/replace CHECKSUM and edit DATE, currently astropy
            # will (errnously) not refer to UTC in the checksum/datasum headers
            # but to the local timezone of the computer
            hduout.writeto(fits_out,checksum=True)

        wcs= astropy.wcs.WCS(hdr)
        # print(wcs)
        return wcs

    def centrTarg(self, site, target, ambi, fib, flen=1839.8, time=None):
        """
        Compute a target (in the usual astrometic sense) by assuming
        the astronomer wants to place the target of the argument list on
        the head of the fiber in the argument list.
        This covers the generic task to let the telescope acquire a
        target that is supposed to land on one of the P1-x or P2-x fiber
        heads of the spectophotometric channel.


        This involves basically calculating the current direction
        to the NCP at the nominal target in the focal plane, calculating
        the position of the fiber head in the focal plane, and using
        the center of the fiber bundle as a distance and position
        angle relativ to the nominal target...
        :param site: location of the observatory
        :type site: lvmtipo.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates

        :param ambi: Ambient data relevant for refractive index
        :type ambi:

        :param fib: One of the fibers in one of the 4 tables
        :type fib: Fiber

        :param flen: Focal length (defining th eplate scale) in mm
                    This should be the same as in the YAML file of the cameras...
        :type flen: float

        :param time: time of the observation /UTC; if None, the current time will be used.
        :type time:
        :return: a new target which is in the fiber bundle center when target is at the fiber
        """

        # direction to NCP from target in the FP orientation (=0 if NCP=up)
        fang=self.fieldAngle(site, target, ambi, time)

        # direction of the fiber relative to the bundle center
        labang=fib.labAngle()
        # coordinates of fiber in microns relative to bundle center
        xfib, yfib=fib.xyFocalPlane()

        # plate scale, microns in the focal plane per radian on the sky
        # image scale at focus is 8.92 um/arcsec = 8.92um/(pi/180/3600 rad) =1.8e^6 um/rad
        pscal=1000.*flen

        # distance center to fiber in radians on the sky
        throwrad=math.hypot(xfib, yfib)/pscal

        # The lab angle of the fiber location to the fiber bundle center
        # (again up=0) is essentially 180 degress plus the labang
        # (quasi by inversion of the role of the coordinate systems of the two)
        # The angle from "up" to NCP is fang, and the angle from "up" to center
        # is 180+labang, so the angle from NCP (clockwise) to center is
        # 180+laban-fang. The position angle on the sky is defined c.c.w.
        # Whether "clockwise" means a position angle
        # (for flipped images) or the negative position angle (for non-flipped)
        # depends on a K-mirror being present or not. The non-flipped
        # images are on the spectrophotom. (names P*) the flipped on
        # the science and background (name S*, A*, B*)
        posang=fang - labang - math.pi
        if fib.name[0:1] != 'P':
            posang *= -1

        # we do not care her to wrap this into [-180..180] or [0..360] or whatever

        # print("field angle target " + str(math.degrees(fang)))
        # print("lab angle fiber " + str(math.degrees(labang)))
        # print("position angle ctr " + str(math.degrees(posang)))
        # print("throw " + str(3600*math.degrees(throwrad)) + "arcsec")

        # apply offset to the given target
        targ_ctr =Target(target.targ.directional_offset_by(
                       posang * astropy.units.rad, throwrad*astropy.units.rad))

        return targ_ctr

    def mpiaMocon(self, site, target, ambi, degNCP=0.0, deltaTime=45.0, polyN=20,
                  time=None, stepsAtStart=135*180*100,
                  homeIsWest=False, homeOffset=135.0, stepsPerturn=360*180*100):
        """
        Compute the polynomial coefficients to rotate the K-mirror for
        a total of polyN*deltaTime seconds in the future with the MPIA
        MoCon, starting at 'time'. The result is a 2dim list of lists in the
        format [[time0,vel0,pol0,acce0,jer0],[time1,vel1,],[],...]

        :param site: location of the observatory
        :type site: lvmtipo.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates.SkyCoord

        :param ambi: Ambient data relevant for refractive index
        :type ambi:

        :param degNCP: The angle in degrees where the NCP (direction of +delta)
              in the field should be fixed in the focal plane, basically
              a fiber selector.
              A value of zero means the +delta is up in the laboratory,
              the value +90 means the +delta direction is right (horizontally E)
              If this is a fiber, the direction from the center of the fiber bundle
              to that fiber is implied.
        :type degNCP: float or lvmtipo.Fiber

        :param deltaTime: Time covered by a single polynomial in seconds
          The default is 45 seconds, which is small enough to keep
          the target aligned with the 0.23 mrad requirement  using only
          first order polynomials. Note that using much smaller time
          intervals may lead to increased download times of trajectories.
        :type deltaTime: float

        :param polyN: Number of polynomials to be constructed.
          The default is 20. The product of polyN and deltaTime
          should at least be as long as the exposure time of the next exposure
          on that optical table/fiber bundle/camera. So defaults of 20
          and defaults of 45 seconds cover the 15 minutes of what is supposed
          to be some standard of the LVM (South). Note that the trajectory
          will stop after that total time; the motor can also be forced to
          stop earlier (which is not in the scope of this documentation or software.)
        :type polyN: int

        :param time: start time of the derotation /UTC, the time when the first
          polynomial in the trajectory should start; if None, the current time will be used.
        :type time:

        :param stepsAtStart: A prefered angle in units of Mocon motor steps for the start
          of the trajectory. The default is to start close to where the K-mirror is "up".
          Background: the parameter degNCP leaves two choices of the motor
          angle to map the sky to the focal plane, separated by 180 degrees,
          because (at least ignoring K-mirror wobbles) the sky rotation is
          twice the mechanical rotation. This parameter allows incremental
          chaining of the trajectories to provide a hint where the previous
          polynomial of the trajectory ended such that the algorithm here can
          connect to that end of the trajectory as smoothly as possible. If
          that option would not exist, an automated choice might "jump" by
          180 degrees (mechanically) to another prefered position which might
          lead to 360 degree rotations of the field around the pointing center
          during exposures.
        :type stepsAtaStart: float

        :param homeIsWest: True if the western of the two limit switches of
              the K-mirror is the home switch, false if the eastern limit
              switch is the home switch. Here "east/west" are the topocentric
              direction at LCO. Because the test setup at MPIA is rotated by
              180 degrees, these meanings are the opposite at the MPIA.
              Default is what's be in fact the cabling at MPIA in Feb 2022 and Nov 2022.
        :type homeIsWest: bool

        :param homeOffset: The angular difference between the K-mirror position
              at home and if its M2 is up, measured in degrees. This value is always positive
              and definitely must be calibrated before this function can be used.
              Default is an estimate from the engineering design, where the
              hall sensor (defining home) is slightly "inside" the mechanical switch.
              The maximum usable range (mechanically) is roughly twice that value,
              because the two limit switches are approximately symmetrical at
              the west and east.
        :type homeOffset: float

        :param stepsPerturn: The number of steps to move the K-mirror
              by 360 degrees. According to information of Lars Mohr of 2021-11-25 we
              have 100 steps per degree, 180 microsteps per step, 360 degrees per turn.
              The product of these 3 numbers defines the default.
        :type stepsPerturn: int

        :return: The list of list of integer values for the MPIA motion controller
              external profile. This is the bare parameter set for the
              SetExternalProfileData command (221). All entries are
              scaled with the applicable powers of 2^16 or 2^32.
              If polyN<=0, that array
              is empty, obviously. If the velocity parameter of many
              consecutive polynomials stays the same within its integer
              representation, you are specifying too short deltaTime values
              without having any benefit (but just longer computation times,
              longer download times to the MoCon, and less predictable
              trajectory start times....)
        .. warning::
          The program does not check that the targets are reachable, which
          means whether they are in a zenith angle < 60 deg or above the horizon
          at that epoch or date.
        """

        # in case the angle is implicitly coded as a fiber, take the fiber
        # position as the desired angle in the focal plane. From here on
        # posang_fp is in radians.
        if isinstance(degNCP, Fiber) :
            posang_fp = degNCP.labAngle()
        elif isinstance(degNCP, astropy.coordinates.Angle) :
            posang_fp = degNCP.radian
        else:
            posang_fp = math.radians(degNCP)

        moc=[]
        if polyN > 0:
            if isinstance(time, astropy.time.Time):
                now=time
            elif isinstance(time, str):
                now=astropy.time.Time(time, format='isot', scale='utc')
            elif time is None:
                now=astropy.time.Time.now()

            tdiff=astropy.time.TimeDelta(deltaTime*astropy.units.second)

            # collect array of bare field angles in radians
            rads=[]
            # number of steps to switch branch of the arctan (avoid +-180 deg wraps)
            degsteps=0
            for poly in range(polyN+1):
                # print(now)
                ang=self.fieldAngle(site, target, ambi, time=now)
                ang += degsteps*2.0*math.pi

                if poly > 0:
                    if ang > rads[poly-1] + math.pi:
                        degsteps -= 1
                        ang -= 2.0*math.pi
                    elif ang < rads[poly-1] - math.pi:
                        degsteps += 1
                        ang += 2.0*math.pi

                rads.append(ang)
                # advance clock to the start of next polynomial
                now += tdiff

            # The K-mirror flips (with its 3-mirror design) the component
            # of incoming images parallel to the (common) incidence plane
            # of the 3 mirrors. The component perpendicular to the incidence
            # plane stays where it is. So the K-mirror angle component
            # perpendicular to the incidence plane bisects the angles of
            # the incoming and outgoing images. Since we want to keep degNCP
            # of the outgoing beam fixed, the bisecting angle (r+degNCP)/2
            # defines where the component perpendicular to the incidence
            # plane is supposed to be. Then another 90 deg rotation
            # defines where the K-mirror "up-down" incidence plane should be,
            # which is the (mechanical) plane of the K-mirror motor.
            rads=[(math.pi + r + posang_fp)/2. for r in rads]

            # Use an optional jump of 180 deg (that's optically 360 deg)
            # to keep trajectory near the angle of stepsAtStart.
            kmirr = Kmirror(home_is_west = homeIsWest, home_offset = homeOffset,
                steps_per_turn=stepsPerturn)
            radsAtStart = kmirr.steps_to_radians(stepsAtStart)

            # if we are too far away from the desired initial position,
            # add or subtract 180 degrees .
            if rads[0] < radsAtStart - 0.5*math.pi:
                rads=[r + math.pi for r in rads]
            elif rads[0] > radsAtStart + 0.5*math.pi:
                rads=[r - math.pi for r in rads]

            # convert all angles from radians to steps
            rads=[ kmirr.radians_to_steps(r) for r in rads]

            # 1 cycle = 614.4 microsecs, see section 9.3 of MoCon User's Guide
            cycsteps=deltaTime/614.4e-6
            for poly in range(polyN):
                # Scale velocity and acceleration with 2^16=65536, yerk with 2^32.
                # We do not use acceleration and yerk (almost 0 for LVMT)
                pos=round(rads[poly])

                # rads[poly+1]-rads[poly]  is velocity in units of counts
                # per deltaTime. Divide by deltaTime to get counts per second
                # and multiply by 614.4e-6 to get counts per cycle.
                vel=round(65536*(rads[poly+1]-rads[poly])/cycsteps)

                # Note the order: duration, velocity is before position....
                traj=[round(cycsteps), vel, pos, 0, 0]
                moc.append(traj)

            # last entry with time 0 to stop the motor
            # (otherwise it would cycle and restart/rewind at entry 0)
            moc.append([0, 0, round(rads[-1]), 0, 0])

        return moc

    def _off_norm2(self, offom1, offom2, motang, signom1= -1, signom2= -1):
        """ Compute a deviation of calculated and measured PW motor angle.
        This is basically the function which defines what is to be minimized
        while searching for best fitting PW motor angle offsets in find_pw_ax_off().

        :param offom1: offset (zero point) for M1 motor axis (radians)
        :param offom2: offset (zero point) for M2 motor axis (radians)
        :param motang: [0]=motor angle ax0 [1]=motor angle ax1 [2]=azimuth [3]=altitude (degrees)
        :param signom1: +1 turns in direction of manuscript, else -1 opposite
        :param signom2: +1 turns in direction of manuscript, else -1 opposite
              Heuristics shows that these must be -1 as of 2022-12-12 tests.
              Which means we assume that the PW motor angles are measuring
              c.c.w.
        """
        # extract motor angles and assumed offsets
        # omega_2 of the manuscript is mot0 in PW and omega_1 is mot1 in PW
        om1 = math.radians(motang[1])-offom1
        om2 = math.radians(motang[0])-offom2
        om1 *= signom1
        om2 *= signom2

        # left hand side equ (39), vector of the motor angles
        # as measured with assumed offset
        om_mot = numpy.zeros((3))
        om_mot[0] = math.cos(om2)*math.sin(om1)
        om_mot[1] = math.sin(om2)*math.sin(om1)
        om_mot[2] = math.cos(om1)

        # extract azimuth and zenith distance from arguments
        az = math.radians(motang[2])

        # extract zenith distance, 90 deg comlemeent of altitude
        Z = math.radians(90-motang[3])

        # right hand side equ (39), vector of the sky coordinates.
        aa = numpy.zeros((3))
        aa[0] = math.sin(az)*math.sin(Z)
        aa[1] = math.cos(az)*math.sin(Z)
        aa[2] = math.cos(Z)
        # computed/predicted left hand side by 3x3 matrix multiplication
        om = numpy.dot(self.B, aa)
        # return a distance between the computed and measured Cartesian directions
        return numpy.linalg.norm(om_mot- om)

    def _off_sum(self, offom1, offom2, motang):
        """ Accumulate over all measured entries in motang the squared
            distance between predicted and computed motor angle (vector)
        :param offom1: offset (zero point) for first motor axis (radians)
        :type offom1: float
        :param offom2: offset (zero point) for second motor axis (radians)
        :type offom2: float
        :param motang: list of 4-tuples
                       [0]=PW angle ax0 [1]=PW angle ax1 [2]=azimuth [3]=altitude (degrees)
        :return: a figure of merit to be minimised to find the best offom1 and offom2
        :rtype: float
        """
        off = 0.0
        for i in range(len(motang)):
            off += self._off_norm2(offom1,offom2, motang[i])
        return off

    def find_pw_ax_off(self, motang) :
        """ Find a pair of motor offsets of the PW motors.
        The task of this calibration routine is to find a pair of offset angles
        to be fed into the __init__ function which matches PW motor axis
        angles to the right-handed angles defined in the paper of the imperfect K-mirrors.

        :param motang:

        This is a list of PW motor angles (ax0 and ax1) and associated
        azimuth and altitude angles, all 4 values in degrees.
        The 4 value should be read off the PW GUI at the same time;
        it is not required that any actual star is observed and this
        should be done with zero offsets in the interface and tracking
        disabled (so the motor angles and horizontal coordinates are fixed.)
        The list has as many 4-tuples as motor and sky coordinates noted down;
        For fitting purposes of the motor angle offsets there ought to be
        at least 2 4-tuples in that list of lists.
        Example test case here is of the 2022-12-12 observations at the MPIA:
        motang = [ [101.5, 54.29,176+54/60.0+34/3600.0, 64+40/60.0+49/3600.0],
            [88.23, 19.44, 308+52/60.0+37.0/3600.0, 74+42/60.0+19.3/3600.0],
            [163.60, 28.87,89+55/60.0+3.8/3600.0, 26+40/60.0+51.5/3600.0] ]
        :type motang: list of 4-tuples of numbers

        :return: two angles in degrees which can be used as om_1_off_ang
                 and om2_off_ang to have a calibrated routine for generating PW motor
                 angles from alt-az angles.
        """
        # a pessimistic starting estimate of the best fit: each value misses by 2*pi
        # the aim of this fitting routine is to get a pair of offsets
        # such that this value becomes almost zero (up to jitter....)
        min_sqr = 6.*len(motang)

        # start with random assumption that both offsets are zero.
        # better_offi are the values found in 4 iterations, and
        # best_offi are the angles/values that are so far globally the best
        better_off1 = better_off2 = 0
        best_off1 = best_off2 = 0

        # the search range of the grid of testing offsets is 360 degrees
        # in both angles (as wide as possible)
        off_range = 2.*math.pi

        # In 4 iterations start with a coarse grid and zoom into the
        # region around the best angles found so far to find a global minimum
        for itr in range(4):
            # subdivide the search range of the angles into
            # 30 samples
            off_step = off_range/30.0
            # in an outer loop for the offset of omega1 and
            # an inner loop for the offset of omega2 get a deviation
            # from the observed value at each subsampled point by
            # calling _off_sum.
            for idx1 in range(30):
                offom1 = best_off1 + (idx1-15)*off_step
                for idx2 in range(30):
                    offom2 = best_off2 + (idx2-15)*off_step
                    this_sqr = self._off_sum(offom1, offom2, motang)
                    if this_sqr < min_sqr :
                        # if this pair of offsets produced a better fit,
                        # remember which angles these were and what their
                        # deviation from the hitherto best fit was.
                        better_off1 = offom1
                        better_off2 = offom2
                        min_sqr = this_sqr

            # new pivotal middel of the next iteration are the
            # angles that reached the minimum deviation
            best_off1 = better_off1
            best_off2 = better_off2
            # print(best_off1,best_off2,min_sqr) # debugging to print progress
            # shrink the search range of the two angles by a factor of 20
            # (so they are slightly larger than the factor of30 to catch
            # cases where the minima are on boundaries of the subintervals)
            off_range /= 20.0
        return math.degrees(best_off1), math.degrees(best_off2)

    def alt_az_to_pw_motang(self,horiz):
        """ Compute PW motor angles given altitude and azimuth of the sky.
        The (x,y,z) components of a unit vector that is defined by using the two
        PW motor angles in a polar coordinate system is given by multiplying
        the 3x3 matrix self.B by the unit vector that is defined by the (x,y,z)
        components of the topocentric horizontal coordinate system of the
        alt/az direction to the star.

        :param horiz: the alt/az of the target
        :type horiz: astropy.coordinates.AltAz

        :return: a pair of PW motor angles, ax0, ax1 in degrees
        :rtype: float
        """

        # generate a vector representation of the pointing
        # exactly as in _off_norm2()

        # extract azimuth and zenith distance in radians
        az = horiz.az.radian
        Z = math.pi/2 - horiz.alt.radian

        # right hand side equ (39), vector of the sky coordinates.
        aa = numpy.zeros((3))
        aa[0] = math.sin(az)*math.sin(Z)
        aa[1] = math.cos(az)*math.sin(Z)
        aa[2] = math.cos(Z)
        # computed/predicted left hand side by 3x3 matrix multiplication
        om = numpy.dot(self.B, aa)
        # convert the vector (cos(omega2)sin(omega1), sin(omega2)sin(omega1)..
        # back to ax0 and ax1 with the offsets
        om1 = math.acos(om[2])
        om[0] /= math.sin(om1)
        om[1] /= math.sin(om1)
        om2 = math.atan2(om[1],om[0])

        # this reverses the two signs (handedness) of the motor axes
        # compatible with _off_norm2(). This here operates what's there in reverse order.
        om1 *= -1
        om2 *= -1
        om1 += self.pw_ax_off[0]
        om2 += self.pw_ax_off[1]
        # normalize both angles to 0..360 deg
        om1 = (om1 + 2*math.pi) % (2*math.pi)
        om2 = (om2 + 2*math.pi) % (2*math.pi)
        return math.degrees(om2), math.degrees(om1)

    def guide_wcs_to_pw_motang(self,wcs_ref, wcs_current, site, ambi=None, time=None):
        """ Convert a WCS pair of the guider to a pair of PW motor angle offsets
        Example:
        >>> fits_inr="/home/mathar/lvm.spec.agcam.center_00000051.fits"
        >>> f_inr = astropy.io.fits.open(fits_inr)
        >>> wcs_inr =  astropy.wcs.WCS(f_inr[0].header)
        >>> fits_inc= "/home/mathar/lvm.spec.agcam.center_00000052.fits"
        >>> f_inc = astropy.io.fits.open(fits_inc)
        >>> wcs_inc =  astropy.wcs.WCS(f_inc[0].header)
        >>> sid=Siderostat(azang=180.0)
        >>> site = Site(name='MPIA')
        >>> ax = sid.guide_wcs_to_pw_motang(wcs_inr,wcs_inc,site)
        >>> print(ax)

        :param wcs_ref: the WCS map where the stars should be.
           This is typically obtained by astrometry at the start of the exposure.
        :type wcs_ref: astropy.wcs.WCS

        :param wcs_current: the WCS map where the stars currently are.
           This is typically obtained in a loop while exposing.
        :type wcs_current: astropy.wcs.WCS

        :param site: geographic ITRF location of the observatory
        :type site: lvmtipo.Site

        :param ambi: Ambient data relevant for refractive index
        :type ambi: lvmtipo.Ambient

        :param time: time of wcs_current either as astropy.time.Time
             or a string in ISO format/UTC; if None, the current time will be used.
        :type time:

        :return: a 2-uple of PW differential motor angles (offsets) for ax0, ax1 in degrees.
            These ought to pull the angles such that the forthcoming WCS
            matches the wcs_ref.
        :rtype: float

        """

        # first step is to get the affine transformation that
        # transforms pixels from one WCS to the other.
        aff_ref = Wcsarith(wcs_ref)
        aff_current = Wcsarith(wcs_current)
        # the affine transwformation that transforms pixels in the current 
        # WCS to pixels in the reference WCS. pixref= W^(-1)ref * W(current)*pixcur
        aff = aff_current / aff_ref
        # print(str(aff)) # debugging

        # we ignore the rotational aspect and use only the shift/translate
        # result aff.shift[0,1]. The rotational aspect in aff.rot would only play a role
        # if one would try to update the K-mirror velocities....
        # The following is basically the same as Wcsarith.offset_px_to_azalt().

        # compute where the reference equatorial ra/dec are
        radec = wcs_ref.all_pix2world([0.],[0.],1,ra_dec_order=True)
        object = astropy.coordinates.SkyCoord(ra=radec[0], dec=radec[1],unit="deg")
        targ_ref = Target(object)

        # compute where the current  equatorial ra/dec are
        # Note that this is using wcs_ref again, not wcs_current.
        radec = wcs_ref.all_pix2world([aff.shift[0]],[aff.shift[1]],1,ra_dec_order=True)
        object = astropy.coordinates.SkyCoord(ra=radec[0], dec=radec[1],unit="deg")
        targ_current = Target(object)

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

        # compute where the reference pixels are in alt/az 
        aa_ref = targ_ref.toHoriz(site,ambi=ambi,time=now)

        # compute where the current pixels are in alt/az 
        aa_current = targ_current.toHoriz(site,ambi=ambi,time=now)

        #compute where the reference pixels are in PWI angles
        pw_ax_ref = self.alt_az_to_pw_motang(aa_ref)

        #compute where the current pixels are in PWI angles
        pw_ax_current = self.alt_az_to_pw_motang(aa_current)
        # print(pw_ax_ref, pw_ax_current) # debugging

        # to be checked: is this the subtraction with the correct sign/order???
        # or do we need to subtract ref from current?
        delta_ax0 = pw_ax_ref[0] - pw_ax_current[0]
        delta_ax1 = pw_ax_ref[1] - pw_ax_current[1] 
        return delta_ax0, delta_ax1

    def nearby_target(self, t, m2=[0, 0, 0]):
        '''
        Compute azimuth and zenith angle of a nearby target following
        https://svn.mpia.de/trac/gulli/lvmt/attachment/wiki/software/horCoordLocTarg.pdf
        :param t: is a vector of the three cartesian coordinates of the target (in the observatory)
                in units of millimeters. t[0] is along East, t[1] along North and t[2] up.
                The origin of the
                coordinates is the midpoint between the middle tables.
        :param m2: is the position of the center of M2 in the same coordinate
                system as t. For computations which involve all 4
                benches call this function in a loop with variable m2[0].
                The distance of M2 to the concrete floor is 1.0 according to
                Req-pier-6 value of the LVM-MPIA-PROC-0002_ProcurementSpec-Siderostat/LVM-MPIA-PROC-0002_Siderostat
        :return: a triple [z,A,sdist] with two angles z and A
                in radians: zenith distance and
                azimuth (north over east) of the direction of t as seen from M1, and
                the distance from M1 to the target in millimeters.

        A prototype of tabulating the alt-az-angles for a series
        of calibration screens:

        roof = CalibScreen()
        for cidx in range(40):
                # loop over 40 positions of the calibration screen west to east, millimeters
                xt = -2000.+cidx*100.0
                # 3D vector of where that would be under the sloping roof
                t = roof.cart_coord(xt)
                print("%f " % xt, end='')
                for tidx in range(4):
                   # tidx =0 for telescope skyw up to xidx=3 for telescope skye
                   # The horizontal distance between the M2 pairs is 170 cm according to
                   # SDSS-V_0111_LVMi_ICD_Telescope_to_Enclosure_Oct21_Rev-revGB.docx
                   xm2 = -2550+1700.0*tidx
                   # position of that M2 on the reference system
                   M2 = [xm2 ,0,0]
                   zA = nearby_target(t,m2=M2)
                   print(" %f %f %f" %
                         (math.degrees(zA[0]),math.degrees(zA[1]),zA[2]), end='')
                print()
        '''

        # vector T measured from M2 to t with 3 Cartesian coordinates
        m2_to_t=numpy.array([t[0]-m2[0], t[1]-m2[1], t[2]-m2[2]], numpy.double)

        # second auxiliary base vector
        bCrossT=numpy.cross(self.b, m2_to_t)
        # length second auxiliary base vector
        lenbCrossT=numpy.linalg.norm(bCrossT)

        # third auxiliary base vector
        bCrossTcrossb=numpy.cross(bCrossT, self.b)

        # coefficient of m along third base vector
        alpha3=(self.m1m2len/lenbCrossT)**2

        # coefficient of m along 2nd base vector
        # intermediate value m^2 - alpha3^2 * length square of 3rd vector
        alpha2=self.m1m2len ** 2 - \
            (alpha3 * numpy.linalg.norm(bCrossTcrossb))**2
        alpha2=math.sqrt(alpha2)/lenbCrossT

        # construct vector m from M2 to M1 from known alpha-coefficients
        m=numpy.add(
            numpy.multiply(alpha2, bCrossT),
            numpy.multiply(alpha3, bCrossTcrossb)
            )


        # compute vector s from M1 to target
        s=numpy.subtract(m2_to_t, m)
        slen=numpy.linalg.norm(s)

        # normalize s to unit length
        s=s/slen

        # zenith angle from 3rd component of normalized s-vector
        z=math.acos(s[2])

        # azimuth from ratio of first and second component of s
        # Note that this is correct: we do NOT use atan2(s[1],s[0]) here.
        A=math.atan2(s[0], s[1])

        # (perhaps return math.pi-z at the end instead, the altitude?)
        # (perhaps return a astropy sky object instead?)
        return [z, A, slen]


class CalibScreen():
    """ A model of the place of the calibration screen under the LVMT roof.
        Default parameters digitized from Fig 2-1 of
        SDSS-V_0111_LVMi_ICD_Telescope_to_Enclosure_Oct21_Rev-revGB.docx
    """
    def __init__(self, height=2184, slope=0.29157, tabl_off=460):
        """ This is a v-shaped geometry that helps to relate
        positions at the ceiling above the 4 telescopes to the M2 mirrors of
        the telescopes.
        :param height: The maximum height of the roof (center) above the plane of
                  the M2 mirrors in mm. Because the M2 mirrors are 1 m above the
                  floor, the default value of 2.184m means that this is 3.184 m
                  above the telescope platform floor.
        :type height: float

        :param slope: The inclination of the roof versus the horizontal
                  in the usual dheight/deast differential sense. The arctan
                  of this gives 0.283 rad = 16.25 deg.
        :type slope: float

        :param tabl_off: There is a sort of natural east-west coordinate origin
                  in the middle between the two M2 of the sci and spec tables
                  (that is 1700/2 = 850 mm away from both). If we use this spot as the
                  zero-value of the x-coordinate, this tablOff is the distance
                  (positive along East) of the top of the roof from that center.
                  The default is that the top is 460 mm away from the origin,
                  which means 850-460mm = 360 mm away from one and 850+460=1310 mm
                  away from the other of the two middle telescopes.
        :type tabl_off: float
        """
        self.height=height
        self.slope=slope
        self.east_off=tabl_off

    def cart_coord(self, x_east_pos):
        """
        The cartesian coordinates of a point that is in the
        the vertical plane of the M2 mirros and x_east_pos millimeters
        to the East of the mid-point between the middle two telescopes.
        :param x_east_pos: the offset relative to the mid-point between
               the middle two telescopes in millimeters.
        :type x_east_pos: float
        :return: a vector of cartesian east, north and up cordinates in millimeters
        """

        # relate x coordinate to the location of the middle of the roof, east = positive
        # assume roof is symmetric east to west
        # So this value is zero if the position is at the maximum
        # height of the inner side of the roof.
        x_rel_roof=abs(x_east_pos-self.east_off)

        # equation of z coordinate is abscisa intersection with 2.184 m above
        # M2 and slope dz/dy
        zcoo=self.height - self.slope*x_rel_roof
        return [x_east_pos, 0, zcoo]
