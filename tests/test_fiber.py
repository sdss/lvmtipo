# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2021-11.21
# @Filename: fiber.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import sys
import math

import pytest

from lvmtipo.fiber import Fiber


def test_S():
    """ Test locations of fibers in the range S1-1 to S3-600
    """
    fib = Fiber("S4-1")
    # self.assertRaises(NameError, fib.labAngle)
    with pytest.raises(NameError):
        fib.labAngle()

    fib = Fiber("S1-1")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx(0.0)

    fib = Fiber("S1-2")
    degs = math.degrees(fib.labAngle())
    assert degs== pytest.approx( 0.0)

    fib = Fiber("S1-24")
    degs = math.degrees(fib.labAngle())
    assert degs== pytest.approx( 0.0)

    fib = Fiber("S1-25")
    degs = math.degrees(fib.labAngle())
    assert degs > 2.0

    fib = Fiber("S1-48")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 60.0)

    fib = Fiber("S1-49")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 90.0)

    fib = Fiber("S1-600")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 60.0)

    fib = Fiber("S1-601")
#    self.assertRaises(NameError == pytest.approx( fib.labAngle)
    with pytest.raises(NameError):
        fib.labAngle()

    fib = Fiber("S2-1")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 120.0)

    fib = Fiber("S2-24")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 120.0)

    fib = Fiber("S2-25")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 120.0, 2.0)

    fib = Fiber("S2-48")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 180.0)

    fib = Fiber("S2-49")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( -150.0)

    fib = Fiber("S2-50")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 180.0)

    fib = Fiber("S2-600")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 180.0)

    #fib = Fiber("S3-1")
    #degs = math.degrees(fib.labAngle())
    #assert degs == pytest.approx( -120.0)

    fib = Fiber("S3-25")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( -120.0)

    fib = Fiber("S3-601")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( -60.0)


def test_P():
    """ Test locations of fibers in the range P1-1 to P2-12
    """
    fib = Fiber("P3-1")
#    self.assertRaises(NameError == pytest.approx( fib.labAngle)
    with pytest.raises(NameError):
        fib.labAngle()

    fib = Fiber("P2-0")
#    self.assertRaises(NameError == pytest.approx( fib.labAngle)
    with pytest.raises(NameError):
        fib.labAngle()

    fib = Fiber("P1-13")
#    self.assertRaises(NameError == pytest.approx( fib.labAngle)
    with pytest.raises(NameError):
        fib.labAngle()

    # Note there are rounded angles in LVMi-0086...
    fib = Fiber("P1-1")
    degs = math.degrees(fib.labAngle())
    print(degs)
    print(pytest.approx(15.61))
    assert degs == pytest.approx(15.61, 0.02)

    fib = Fiber("P1-4")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 104.39, 0.02)
    
    # Note angles in LVMi-0086... are 0..360 and here -180..180
    fib = Fiber("P1-10")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( -75.61, 0.02)

    fib = Fiber("P2-1")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 14.56, 0.02)

def test_A():
    """ Test locations of fibers in the range A1-1 to B2-20
    """
    fib = Fiber("A1-18")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 90.0)

    fib = Fiber("B1-18")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 90.0)

    fib = Fiber("A2-20")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( 180.0)

    fib = Fiber("B3-21")
    degs = math.degrees(fib.labAngle())
    assert degs == pytest.approx( -60.0)

