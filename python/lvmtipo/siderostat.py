# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-11-04
# @Filename: siderostat.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for siderostat field angles using homogeneous coordinates
"""

import math
import numpy
import astropy.coordinates
import astropy.time
import astropy.units

from lvmtipo.mirror import Mirror
from lvmtipo.target import Target
from lvmtipo.kmirror import Kmirror

__all__ = ['Siderostat']


class Siderostat():
    """ A siderostat of 2 mirrors
    """

    def __init__(self, zenang=90.0, azang=0.0, medSign=-1, m1m2dist=240.0):
        """ A siderostat of 2 mirrors
        :param zenang: Zenith angle of the direction of the exit beam (degrees)
                   in the range 0..180. Default is the design value of the LCO LVMT.
        :type zenang: float
        :param azang: Azimuth angle of the direction of the exit beam (degrees)
                   in the range -180..360 degrees, N=0, E=90.
                   Ought to be zero for the LCO LVMT where the FP is north of the
                   siderostat and 180 for the MPIA test setup where the FP is
                   south of the siderostat..
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

        # axes orthogonal to beam
        self.box = numpy.zeros((3))
        self.box[0] = 0.0
        self.box[1] = -self.b[2]
        self.box[2] = self.b[1]
        self.boy = numpy.cross(self.b, self.box)

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
        m2tom1 = numpy.cross(star, self.b)
        vlen = numpy.linalg.norm(m2tom1)
        m2tom1 /= self.sign*vlen

        # surface normal to M1
        m1norm = star - m2tom1
        vlen = numpy.linalg.norm(m1norm)
        m1norm /= vlen
        m1 = Mirror(m1norm, 1.0)

        # surface normal to M2
        m2norm = self.b + m2tom1
        vlen = numpy.linalg.norm(m2norm)
        m2norm /= vlen
        m2 = Mirror(m2norm, 0.0)

        # transformation matrix for the 2 reflections
        m1trans = m1.toHomTrans()
        m2trans = m2.toHomTrans()
        trans = m2trans.multiply(m1trans)

        # for the field angle need a target that is just a little bit
        # more north (but not too little to avoid loss of precision)
	# 10 arcmin = 0.16 deg further to NCP
        targNcp = Target(target.targ.spherical_offsets_by(
                       astropy.coordinates.Angle("0deg"), astropy.coordinates.Angle("0.16deg")))
        horizNcp = targNcp.toHoriz(site=site, ambi=ambi, time=now)

        starNcp = numpy.zeros((3))
        # same procedure as in the construction of b in the Sider ctor,
        # but with 90-zenith angle=altitude
        starNcp[0] = math.sin(horizNcp.az.radian) * \
                              math.cos(horizNcp.alt.radian)
        starNcp[1] = math.cos(horizNcp.az.radian) * \
                              math.cos(horizNcp.alt.radian)
        starNcp[2] = math.sin(horizNcp.alt.radian)

        # image of targNcp while hitting M1
        m1img = trans.apply(m2tom1)
        # image of targNcp before arriving at M1
        starOffm1 = m2tom1 + starNcp
        starimg = trans.apply(starOffm1)

        # virtual direction of ray as seen from point after M2
        # no need to normalize this because the atan2 will do...
        starvirt = m1img - starimg

        # project in a plane orthogonal to  self.b
        cosFang = numpy.dot(starvirt, self.box)
        sinFang = numpy.dot(starvirt, self.boy)
        return math.atan2(sinFang, cosFang)

    def to_header(self, site, target, camera, ambi, k_mocon_steps, ag_cam_name,
           genrevx, genrevy, kmirr, time=None,
           dist_cam_edge=11.14771):
        """ Convert the parameters to FITS WCS header cards.
        This traces the position angle implied by the siderostat and number of mirror
        reflections through the K-mirror if applicable and through the
        prisms if applicable to predict the key WCS keywords. The main difference
        to the implementation in lvmcam/python/lvmcam/models/wcs.py
        is that the variable genvrev[xy] are included in this computation.
        .. warn:: this assumes that the FITS files are read out as
           in the MPIA test setup: all east/west/center cameras with the
           long edge up-down, no rotation by any 90 or 180 degrees in any
           python software that may spring up after Oct 2022.

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
             For one of the three benches (hidden in ag_cam_name) this is
             a dummy argument because there is no K-mirror.
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

        :param kmirr: K-mirror constants related to conversions of steps to angles
        :type kmirr: lvmtipo.Kmirror

        :param dist_cam_edge: the distance of the long east or west FLIR camera edge to
             the fiber bundle centre in millimeters
             Is 7.144+4.0 mm according to SDSS-V_0110 figure 6
             and 11.14471 according to figure 3-1 of LVMi-0081
        :type dist_cam_edge: float

        :return a set of FITS header cards with WCS keywords.
        :rtype astropy.io.fits.Header
        """
        # start with an empty set of header cards
        wcshdr = astropy.io.fits.Header()

        # compute the position angle (NCP) after the siderostat
        # and before the K-mirror in radians.
        if isinstance(time, astropy.time.Time):
            now = time
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()
        ang = self.fieldAngle(site, target, ambi, time=now)

        # number of mirror reflections so far (2 of the siderostat)
        n_refl = 2

        # Simplify the interface by also recognizing the descr of
        # the camera (which uses upper case for some letter).
        ag_cam_name = ag_cam_name.lower()
        # three LVMT benches have a K-mirror, one has not
        if ag_cam_name.find('skye') >= 0 or ag_cam_name.find('skyw') >= 0 or ag_cam_name.find('sci') >= 0:
            n_refl += 3

            # convert K-mirror orientation to radians
            kangle = kmirr.steps_to_radians(k_mocon_steps)

            # flip the component of the position angle in the incidence
            # plane(s) of the K-mirror. So kangle+-90 deg is the
            # arithmetic mean of the position angle before and after the K-mirror.
            # (ang_in+ang_out)/2 = kangle+-pi/2
            # ang_out = 2*kangle+-pi - ang_in
            ang = math.pi + 2.*kangle - ang
        elif ag_cam_name.find('spec') >= 0:
            pass
        else:
            raise NameError("Unrecognized camera name " + ag_cam_name)

        if ag_cam_name.find('center') >= 0 or ag_cam_name.find('agc') >= 0:
            # no prism
            # assume camera is mounted vertically with the long lower edge
            # pointing east ( basically rotating the east camera 90 deg
            # around the vertical axis and then 180 deg around the horizontal axis.
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "left" or "to the west"
            # (passive rotation)
            ang += math.pi/2.0
        elif ag_cam_name.find('west') >= 0 or ag_cam_name.find('agw') >= 0:
            # prism incidence plane is horizontal. Flip that component
            ang *= -1.0
            n_refl += 1
            # assume camera is mounted vertically with the long lower edge
            # pointing south, away from the focal plane
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "right" or geographicall "to the north"
            # (passive rotation)
            ang -= math.pi/2.0
        elif ag_cam_name.find('east') >= 0 or ag_cam_name.find('age') >= 0:
            # prism incidence plane is horizontal. Flip that component
            ang *= -1.0
            n_refl += 1
            # assume camera is mounted vertically with the long lower edge
            # pointing north (basically rotated 180 deg around the vertical
            # axis relative to the west camera) toward the focal plane
            # We want to express an angle where the reference is "up" to an
            # angle where the reference "right" or geographically "to the south"
            # (passive rotation)
            ang += math.pi/2.0
        else:
            raise NameError("Unrecognized camera name " + ag_cam_name)

        # assign an image parity defined by the number of mirror reflections
        if (n_refl % 2) == 0:
            imgparity = True
        else:
            imgparity = False

        # up to here ang is refering to the short-edge-up orientation
        # of the camera if that were read out by streaming the araviscam
        # bytes in standard order (genrevx=gerevy=false). Because we need
        # one flip up-down (first row becomes last row) to flush the stream
        # into the standard FITS order where the lowest row is the first to be
        # read, the optically relevant value for that flip is opposite to the
        # nominal value used for the araviscam itnerface
        genrevy = not genrevy

        # we're left with 4 combinations of genrev[xy]
        if genrevx:
            if genrevy:
                # both flipped (i.e., a 180 deg rotation)
                ang += math.pi
            else:
                # x flipped
                ang = -ang
                imgparity = not imgparity
        else:
            if genrevy:
                # x not flipped but y flipped
                # Note that this is the operation as above for the kmirror at angle 0
                # a y-flip is the same as an x-flip plus a 180 deg rotation
                ang = math.pi - ang
                imgparity = not imgparity

        # add the main pointing keywords
        key = astropy.io.fits.Card("CUNIT1", "deg", "WCS units along axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CUNIT2", "deg", "WCS units along axis 2")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE1", "RA---TAN", "WCS type axis 1")
        wcshdr.append(key)
        key = astropy.io.fits.Card("CTYPE2", "DEC--TAN", "WCS type axis 2")
        wcshdr.append(key)
        key = astropy.io.fits.Card(
            "CRVAL1", target.ra.degree, "[deg] RA at reference pixel")
        wcshdr.append(key)
        key = astropy.io.fits.Card(
            "CRVAL2", target.dec.degree, "[deg] Dec at reference pixel")
        wcshdr.append(key)

        # convert disgtcamEdgeCtr from mm to pixels
        # this is the distance from the long edge that is narrowest to
        # the fiber center after projection into the focal plane.
        dist_cam_edge *= 1000.0 / (camera.pixsize*camera.binning[1])
        # distance from the middle of the east or west camera to fiber ctner in pixels
        dist_cam_mid = camera.detector_size.hd / \
            2 / camera.binning[1] + dist_cam_edge

        # where is the center of image away in the pixel coordinate system
        # of the camera? For the x-position this is in the middel of the
        # camera along the long axis because all cameras are installed up-down
        crpix1 = camera.detector_size.wd / 2 / camera.binning[0] + 0.5
        crpix2 = camera.detector_size.hd / 2 / camera.binning[1] + 0.5
        if ag_cam_name.find('center') >= 0 or ag_cam_name.find('agc') >= 0 :
            # ra/dec in the center of the image
            pass
        elif ag_cam_name.find('west') >= 0 or ag_cam_name.find('agw') >= 0 :
            # the direction from the prism center to the fiber center
            # is walkign east, which is in the camera walking to the lower edge
            if genrevy:
                crpix2 += dist_cam_mid
            else:
                crpix2 -= dist_cam_mid
        elif ag_cam_name.find('east') >= 0 or ag_cam_name.find('age') >= 0 :
            # the direction from the prism center to the fiber center
            # is walkign west, which is in the camera walking to the upper edge
            if genrevy:
                crpix2 -= dist_cam_mid
            else:
                crpix2 += dist_cam_mid

        key=astropy.io.fits.Card(
            "CRPIX1", crpix1, "[px] point cntr along axis 1")
        wcshdr.append(key)
        key=astropy.io.fits.Card(
            "CRPIX2", crpix2, "[px] point cntr along axis 2")
        wcshdr.append(key)

        # the two main values in the CD matrix, irrespective of sign
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

        # negative CD determinant for uprigt images (imgparity=True), because
        # (ra/dec) are a left-handed coordinate system; positive for flipped (imparity=False)
        if imgparity:
            key=astropy.io.fits.Card(
                "CD1_1", -cosperpix, "[deg/px] WCS matrix diagonal")
            wcshdr.append(key)
            key=astropy.io.fits.Card(
                "CD1_2", sinperpix, "[deg/px] WCS matrix outer diagonal")
            wcshdr.append(key)
        else:
            # in this case the alpha coordinate values are in the opposite direction
            # so CD1_1 and CD1_2 are negated compared to above
            key=astropy.io.fits.Card(
                "CD1_1", cosperpix, "[deg/px] WCS matrix diagonal")
            wcshdr.append(key)
            key=astropy.io.fits.Card(
                "CD1_2", -sinperpix, "[deg/px] WCS matrix outer diagonal")
            wcshdr.append(key)

        key=astropy.io.fits.Card(
            "CD2_1", sinperpix, "[deg/px] WCS matrix outer diagonal")
        wcshdr.append(key)
        key=astropy.io.fits.Card(
            "CD2_2", cosperpix, "[deg/px] WCS matrix diagonal")
        wcshdr.append(key)

        return wcshdr


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
        :type site: fieldrotation.Site
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
        :return a new target which is in the fiber bundle center when target is at the fiber
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
        :type site: fieldrotation.Site

        :param target: sidereal target in ra/dec
        :type target: astropy.coordinates.SkyCoord

        :param ambi: Ambient data relevant for refractive index
        :type ambi:

        :param degNCP: The angle in degrees where the NCP (direction of +delta)
              in the field should be fixed in the focal plane, basically
              a fiber selector.
              A value of zero means the +delta is up in the laboratory,
              the value +90 means the +delta direction is right (horizontally E)
        :type degNCP: float

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
            rads=[(math.pi + r + math.radians(degNCP))/2. for r in rads]

            # Use an optional jump of 180 deg (that's optically 360 deg)
            # to keep trajectory near the angle of stepsAtStart.
            kmirr = Kmirror(home_is_west = homeIsWest, home_offset = homeOffset, steps_per_turn=stepsPerturn)
            radsAtStart= kmirr.steps_to_radians( stepsAtStart)

            if rads[0] < radsAtStart - 0.5*math.pi:
                rads=[r + math.pi for r in rads]
            elif rads[0] > radsAtStart + 0.5*math.pi:
                rads=[r - math.pi for r in rads]

            # convert all angles from radians to steps
            rads=[r * stepsPerturn/(2.0*math.pi) for r in rads]

            homeOffsetSteps=homeOffset * stepsPerturn / 360.0
            # convert all counts measured from 0=up to counts
            # relative to the home/reference position of the MoCon
            if homeIsWest:
                rads=[r + homeOffsetSteps for r in rads]
            else:
                rads=[homeOffsetSteps - r for r in rads]


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
        T=numpy.array([t[0]-m2[0], t[1]-m2[1], t[2]-m2[2]], numpy.double)

        # second auxiliary base vector
        bCrossT=numpy.cross(self.b, T)
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
        s=numpy.subtract(T, m)
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
    def __init__(self, height=2184, slope=0.29157, tablOffs=460):
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

        :param tablOffs: There is a sort of natural east-west coordinate origin
                  in the middle between the two M2 of the sci and spec tables
                  (that is 1700/2 = 850 mm away from both). If we use this spot as the
                  zero-value of the x-coordinate, this tablOff is the distance
                  (positive along East) of the top of the roof from that center.
                  The default is that the top is 460 mm away from the origin,
                  which means 850-460mm = 360 mm away from one and 850+460=1310 mm
                  away from the other of the two middle telescopes.
        :type tablOffs: float
        """
        self.height=height
        self.slope=slope
        self.eastOffs=tablOffs

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
        x_rel_roof=abs(x_east_pos-self.eastOffs)

        # equation of z coordinate is abscisa intersection with 2.184 m above
        # M2 and slope dz/dy
        zcoo=self.height - self.slope*x_rel_roof
        return [x_east_pos, 0, zcoo]
