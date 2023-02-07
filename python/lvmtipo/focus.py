# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: lvmtipo/focus.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np

tempdata = np.array([
    [8.98, 35.412],
    [12.99, 38],
    [17.7, 41.0],
])

def temp2focus( temperature:float, polynomial:float = 2):
    """Gathering focus based on temperature"""

    d = tempdata.T
    foc_est = np.polyfit(d[0], d[1], polynomial)

    return np.polyval(foc_est, temperature)
