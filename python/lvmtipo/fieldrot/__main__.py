#!/usr/bin/env python3

"""
# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: homtrans.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

Python3 class for siderostat field angles using homogeneous coordinates
"""

import sys
import math
import numpy
import astropy.coordinates
import astropy.time
import astropy.units

from lvmtipo.site import Site
from lvmtipo.siderostat import Siderostat
from lvmtipo.fiber import Fiber
from lvmtipo.target import Target


def main():
    """ Example application demonstrating the interface.
    .. todo demonstrate use of proper motions 
    """
    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument("-v", '--verbose', action='store_true',
    #                     help="print some notes to stdout")

    # right ascension in degrees
    parser.add_argument("-r", '--ra', help="RA J2000 in degrees or in xxhxxmxxs format")

    # declination in degrees
    parser.add_argument("-d", '--dec', help="DEC J2000 in degrees or in +-xxdxxmxxs format")

    # shortcut for site coordinates: observatory
    parser.add_argument("-s", '--site', default="LCO", help="LCO or MPIA or APO or KHU")

    # optional name of a fiber
    parser.add_argument("-f", '--fiber', help="fiber name like P2-1 or P1-11")

    args = parser.parse_args()

    # check ranges and combine ra/dec into a astropy SkyCoord
    if args.ra is not None and args.dec is not None :
        if args.ra.find("h") < 0 :
            # apparently simple floating point representation
            targ = astropy.coordinates.SkyCoord(ra=float(args.ra), dec=float(args.dec),unit="deg")
        else :
            targ = astropy.coordinates.SkyCoord(args.ra + " " + args.dec)
    else :
        targ = None

    # step 1: define where the observatory is (on Earth)
    geoloc = Site(name = args.site)
    # print(geoloc)

    # step 2: define where the output beam of the siderostat points to
    sid = Siderostat()
    # print(sid)

    # step 3: define where the sidereostat is pointing on the sky
    point = Target(targ)
    print("target is ",targ)

    # calculate the field angle (in radians)
    rads = sid.fieldAngle(geoloc, point, None)
    print("field angle " + str(math.degrees(rads)) + " deg")

    # if a P[12]-[1..12] fiber head was specified, calculate
    # also the virtual target in the fiber bundle center.
    if args.fiber is not None and targ is not None :
        fib=Fiber(args.fiber)
        # print("lab angle " + str(math.degrees(fib.labAngle())) + " deg")
        ctrTarg = sid.centrTarg(geoloc, point, None, fib)
        print(ctrTarg.targ)

if __name__ == "__main__":
    main()

