# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: siderostat.py
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

from lvmtipo.mirror import Mirror
from lvmtipo.target import Target

__all__ = ['Siderostat']


class Siderostat():
    """ A siderostat of 2 mirrors
    """
    def __init__(self, zenang=90.0, azang=0.0, medSign=1) :
        """ A siderostat of 2 mirrors
        :param zenang Zenith angle of the direction of the exit beam (degrees)
                   in the range 0..180. Default is the design value of the LCO LVMT.
        :type zenang float
        :param azang Azimuth angle of the direction of the exit beam (degrees)
                   in the range -180..360 degrees, N=0, E=90.
                   Default is the design value of the LCO LVMT.
        :type azang float
        :param medSign Sign of the meridian flip design of the mechanics.
                       Must be either +1 or -1. Default is the LCO LVMT design (in most
                       but not all of the documentation).
        :type medSign int
        """

        if isinstance(zenang, (int,float) ) and isinstance(azang, (int,float)) :
            self.b = numpy.zeros((3))
            self.b[0] = math.sin( math.radians(azang)) * math.sin( math.radians(zenang))
            self.b[1] = math.cos( math.radians(azang)) * math.sin( math.radians(zenang))
            self.b[2] = math.cos( math.radians(zenang))
        else :
            raise TypeError("invalid data types")

        if isinstance(medSign,int) :
            if medSign == 1 or medSign == -1:
                self.sign = medSign
            else:
                raise ValueError("invalid medSign value")
        else :
            raise TypeError("invalid medSign data type")

        # axes orthogonal to beam
        self.box = numpy.zeros((3))
        self.box[0] = 0.0
        self.box[1] = -self.b[2]
        self.box[2] = self.b[1]
        self.boy = numpy.cross(self.b,self.box)
        

    def fieldAngle(self, site, target, ambi, wlen=0.5, time=None) :
        """ 
        :param site location of the observatory
        :type site fieldrotation.Site
        :param target sidereal target in ra/dec
        :type target astropy.coordinates
        :param ambi Ambient data relevant for refractive index
        :type ambi
        :param wlen wavelength of observation in microns.
             Default is 500 nm, the canonical light in the visible and an
             estimator of the mean band widht of the LVM.
        :type wlen float
        :param time time of the observation /UTC; if None, the current time will be used.
        :type time
        :return field angle (direction to NCP) in radians
        """
        if isinstance(time, astropy.time.Time) :
            now = time 
        elif isinstance(time, str):
            now = astropy.time.Time(time, format='isot', scale='utc')
        elif time is None:
            now = astropy.time.Time.now()

        # copute mirror positions 
        horiz = target.toHoriz(site=site, ambi=ambi, time = now) 
        # print(horiz)

        star = numpy.zeros((3))
        # same procedure as in the construction of b in the Sider ctor, but with 90-zenang=alt
        star[0] = math.sin( horiz.az.radian) * math.cos( horiz.alt.radian )
        star[1] = math.cos( horiz.az.radian ) * math.cos( horiz.alt.radian)
        star[2] = math.sin( horiz.alt.radian)

        # unit vector from M2 to M1
        m2tom1 = numpy.cross(star,self.b)
        len = numpy.linalg.norm(m2tom1)
        m2tom1 /= self.sign*len

        # surface normal to M1
        m1norm = star - m2tom1
        len = numpy.linalg.norm(m1norm)
        m1norm /= len
        m1 = Mirror(m1norm,1.0)

        # surface normal to M2
        m2norm = self.b + m2tom1
        len = numpy.linalg.norm(m2norm)
        m2norm /= len
        m2 = Mirror(m2norm,0.0)

        # transformation matrix for the 2 reflections
        m1trans = m1.toHomTrans()
        m2trans = m2.toHomTrans()
        trans = m2trans.multiply(m1trans)

        # for the field angle need a target that is just a little bit
        # more north (but not too little to avoid loss of precision)
	# 10 arcmin = 0.16 deg further to NCP
        targNcp = Target(target.targ.spherical_offsets_by(
                       astropy.coordinates.Angle("0deg"), astropy.coordinates.Angle("0.16deg")))
        horizNcp = targNcp.toHoriz(site=site, ambi=ambi, time = now) 

        starNcp = numpy.zeros((3))
        # same procedure as in the construction of b in the Sider ctor, 
        # but with 90-zenith angle=altitude
        starNcp[0] = math.sin( horizNcp.az.radian) * math.cos( horizNcp.alt.radian )
        starNcp[1] = math.cos( horizNcp.az.radian ) * math.cos( horizNcp.alt.radian)
        starNcp[2] = math.sin( horizNcp.alt.radian)

        # image of targNcp while hitting M1
        m1img = trans.apply(m2tom1)
        # image of targNcp before arriving at M1
        starOffm1 = m2tom1 + starNcp
        starimg = trans.apply(starOffm1)

        # virtual direction of ray as seen from point after M2
        # no need to normalize this because the atan2 will do...
        starvirt = m1img - starimg 

        # project in a plane orthogonal to  self.b
        cosFang = numpy.dot(starvirt,self.box)
        sinFang = numpy.dot(starvirt,self.boy)
        return math.atan2(sinFang, cosFang)

    def centrTarg(self, site, target, ambi, fib, flen = 1839.8, wlen=0.5, time=None) :
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
        :param site location of the observatory
        :type site fieldrotation.Site
        :param target sidereal target in ra/dec
        :type target astropy.coordinates
        :param ambi Ambient data relevant for refractive index
        :type ambi
        :param fib One of the fibers in one of the 4 tables
        :type Fiber
        :param flen Focal length (defining th eplate scale) in mm
                    This should be the same as in the YAML file of the cameras...
        :type float
        :param wlen wavelenght of observation in microns
        :type wlen float
        :param time time of the observation /UTC; if None, the current time will be used.
        :type time
        :return a new target which is in the fiber bundle center when target is at the fiber
        """

        # direction to NCP from target in the FP orientation (=0 if NCP=up)
        fang  = self.fieldAngle(site, target, ambi, wlen, time)

        # direction of the fiber relative to the bundle center
        labang = fib.labAngle()
        # coordinates of fiber in microns relative to bundle center
        xfib, yfib = fib.xyFocalPlane() 

        # plate scale, microns in the focal plane per radian on the sky
        # image scale at focus is 8.92 um/arcsec = 8.92um/(pi/180/3600 rad) =1.8e^6 um/rad
        pscal = 1000.*flen

        # distance center to fiber in radians on the sky
        throwrad = math.hypot( xfib, yfib)/pscal

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
        posang = fang -labang -math.pi
        if fib.name[0:1] != 'P' :
            posang *= -1

        # we do not care her to wrap this into [-180..180] or [0..360] or whatever

        # print("field angle target " + str(math.degrees(fang)))
        # print("lab angle fiber " + str(math.degrees(labang)))
        # print("position angle ctr " + str(math.degrees(posang)))
        # print("throw " + str(3600*math.degrees(throwrad)) + "arcsec")

        # apply offset to the given target
        targCtr = Target(target.targ.directional_offset_by(
                       posang *astropy.units.rad, throwrad*astropy.units.rad))

        return targCtr

    def mpiaMocon(self, site, target, ambi, degNCP=0.0, deltaTime =45.0, polyN=20, 
                  wlen=0.5, time=None, stepsPerturn = 6480000) :
        """ 
        Compute the polynomial coefficients to rotate the K-mirror for
        a total of polyN*deltaTime seconds in the future with the MPIA
        MoCon, starting at 'time'. The result is a 2dim list of lists in the
        format [[index0,time0,vel0,pol0,acce0,jer0],[index1,time1,vel1,],[],...]

        :param site location of the observatory
        :type site fieldrotation.Site

        :param target sidereal target in ra/dec
        :type target astropy.coordinates

        :param ambi Ambient data relevant for refractive index
        :type ambi

        :param degNCP The angle in degrees where the NCP (direction of +delta)
              in the field should be fixed in the focal plane, basically 
              a fiber selector. 
              A value of zero means the +delta is up in the laboratory,
              the value +90 means the +delta direction is right (horizontally E)
        :type float

        :param deltaTime Time covered by a single polynomial in seconds
          The default is 45 seconds, which is small enough to keep
          the target aligned with the 0.23 mrad requirement  using only
          first order polynomials. Note that using much smaller time
          intervals may lead to increased download times of trajectories.
        :type float

        :param polyN Number of polynomials to be constructed.
          The default is 20. The product of polyN and deltaTime
          should at least be as long as the exposure time of the next exposure
          on that optical table/fiber bundle/camera. So defaults of 20
          and defaults of 45 seconds cover the 15 minutes of what is supposed
          to be some standard of the LVM (South). Note that the trajectory
          will stop after that total time; the motor can also be forced to
          stop earlier (which is not in the scope of this documentation or software.)
        :type int

        :param wlen wavelength of observation in microns
        :type wlen float

        :param time start time of the derotation /UTC; if None, the current time will be used.
        :type time

        :param stepsPerturn The number of steps to move the K-mirror
              by 360 degrees. According to information of Lars Mohr of 2021-11-25 we
              have 100 steps per degree, 18000 microsteps per degree, which
              defines the default.
        :type int

        :return The list of list of integer values for the MPIA motion controller
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
        Notes:
          The program does not check that the targets are reachable, which
          means whether they are in a zenith angle < 60 deg or above the horizon
          at that epoch or date.
        """
        moc = []
        if polyN > 0:
            if isinstance(time, astropy.time.Time) :
                now = time 
            elif isinstance(time, str):
                now = astropy.time.Time(time, format='isot', scale='utc')
            elif time is None:
                now = astropy.time.Time.now()

            tdiff = astropy.time.TimeDelta(deltaTime*astropy.units.second)

            # collect array of bare field angles in radians
            rads = []
            # number of steps to switch branch of the arctan (avoid +-180 deg wraps)
            degsteps =0 
            for poly in range(polyN+1):
                # print(now)
                ang = self.fieldAngle(site, target, ambi, wlen=wlen, time=now)
                ang += degsteps*2.0*math.pi

                if poly > 0 :
                    if ang > rads[poly-1] + math.pi :
                        degsteps -= 1
                        ang -= 2.0*math.pi 
                    elif ang < rads[poly-1] - math.pi :
                        degsteps += 1
                        ang += 2.0*math.pi 

                rads.append(ang)
                # advance clock to the start of next polynomial
                now += tdiff
 
            # same calcuation as in the araviscam.singleFram() function
            # action of mirror is (cos(2*m), sin(2*m);sin(2*m), -cos(2*m))
            # where m is the mechanical mirror angle (0 if M2 is up), and
            # if f is the direction of NCP (+delta) before the K-mirror,
            # the x and y projections are (sin(f),cos(f)) and
            # the direction is ( sin(2m+f), -cos(2m+f)) after the K-mirror
            # = (sin(Pi-2m-f), cos(pi-2m-f)). Setting pi-2m-f= decNcp gives
            # angle of 2m=pi-f-decNcp (optical, before division thru 2).
            # The sign flip r->-r means we are DErotating, and the division
            # through 2 is the usual optical-to-mechanical rotation angle factor.
            # This applies factors -1/2 (sign to derotate, 1/2 for opt-mech-transf)
            rads = [ (math.pi -r -math.radians(degNCP))/2.  for r in rads]

            # Use an arbitrary jump of 180 deg (that's optically 360 deg)
            # to keep trajectory near the angle of 0 (stay away from
            # the stops/limit switches at +-137 deg)
            if rads[0] < -0.5*math.pi:
                rads = [r + math.pi for r in rads]
            elif rads[0] > 0.5*math.pi:
                rads = [r - math.pi for r in rads]

            # convert all angles from radians to counts
            rads = [r *stepsPerturn/(2.0*math.pi) for r in rads]

            # 1 cycle = 614.4 microsecs, see section 9.3 of MoCon User's Guide
            cycsteps = deltaTime/614.4e-6
            for poly in range(polyN):
                # Scale velocity and acceleration with 2^16=65536, yerk with 2^32.
                # We do not use acceleration and yerk (almost 0 for LVMT)
                pos = round(rads[poly])

                # rads[poly+1]-rads[poly]  is velocity in units of counts
                # per deltaTime. Divide by deltaTime to get counts per second
                # and multiply by 614.4e-6 to get counts per cycle.
                vel = round(65536*(rads[poly+1]-rads[poly])/cycsteps)
                # Note the order: index, duration, velocity is before position....
                traj = [round(cycsteps), vel, pos, 0,0]
                moc.append(traj)

            # last entry with time 0 to stop the motor
            # (otherwise it would cycle and restart/rewind at entry 0)
            moc.append([0,0,round(rads[-1]),0,0])
 
        return moc

