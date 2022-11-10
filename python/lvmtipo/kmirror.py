# -*- coding: utf-8 -*-
#
# @Author: Richard J. Mathar <mathar@mpia.de>
# @Date: 2022-11-07
# @Filename: kmirror.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

"""
Python3 class for K-mirror constants and angle convention
"""

import math

__all__ = ['Kmirror']

class Kmirror():
    """ Functions to represent K-mirror angles as Mocon steps
    """

    def __init__(self, home_is_west = False, home_offset = 135.0, steps_per_turn=360*180*100):
        """
        :param home_is_west: True if the western of the two limit switches of
               the K-mirror is the home switch, false if the eastern limit
               switch is the home switch. Here "east/west" are the topocentric
               direction at LCO. Because the test setup at MPIA is rotated by
               180 degrees, these meanings are the opposite at the MPIA.
               Default is what's be in fact the cabling at MPIA in Feb 2022 and Nov 2022.
        :type home_is_west: bool

        :param home_offset: The angular difference between the K-mirror position
               at home and if its M2 is up, measured in degrees. This value is always positive
               and definitely must be calibrated before this function can be used.
               Default is an estimate from the engineering design, where the
               hall sensor (defining home) is slightly "inside" the mechanical switch.
               The maximum usable range (mechanically) is roughly twice that value,
               because the two limit switches are approximately symmetrical at
               the west and east.
        :type home_offset: float

        :param steps_per_turn: The number of steps to move the K-mirror
              by 360 degrees. According to information of Lars Mohr of 2021-11-25 we
              have 100 steps per degree, 180 microsteps per step, 360 degrees per turn.
              The product of these 3 numbers defines the default.
        :type steps_per_turn: int
        """
        self.home_is_west  = home_is_west
        self.home_offset  = home_offset
        self.steps_per_deg  = steps_per_turn/360

    def steps_to_radians(self, steps) :
        """
        Convert Mocon steps to position angle of the middle mirror in radians.

        :param steps: The Mocon steps
        :type steps: int

        :return: The position angle in radians equivalent to the steps.
              Note that some TaN functions in the LVM package may return the oppositely signed angle.
        :rtype: float
        """
        if self.home_is_west :
            # angle = -135 + steps*degreesperstep
            ang = -self.home_offset + steps / self.steps_per_deg
        else:
            # angle = 135 - steps*degreesperstep
            ang = self.home_offset - steps / self.steps_per_deg
        return math.radians(ang)

    def radians_to_steps(self, rads) :
        """
        Convert radians of position angle to Mocon steps
        Inverse function to steps_to_radians().

        :param rads: position angle in radians. 0 if middle mirror up.
        :type rads: float

        :return: The steps  equivalent to the position angle
              Note that some TaN functions in the LVM package may return the oppositely signed angle.
        :rtype: int
        """
        if self.home_is_west :
            # angle = -135 + steps*degreesperstep
            # steps = (angle +135)/degreesperstep
            stp = (math.degrees(rads) + self.home_offset)*self.steps_per_deg
        else:
            # angle = 135 - steps*degreesperstep
            # steps = (135- angle)/degreeperstep
            stp = (self.home_offset - math.degrees(rads))*self.steps_per_deg
        return int(stp+0.5)
