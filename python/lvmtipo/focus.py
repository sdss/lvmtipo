# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: lvmtipo/temp2focus.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import numpy as np

# TODO: put real values from the fits header into here.
# these values are from Heidelberg from science and to spec 2 has been added
temp_data = {
    "lvm.sci.foc":
        np.array([
            [8.98, 35.412],
            [12.99, 38],
            [17.7, 41.0],
        ]),
    "lvm.skye.foc":
        np.array([
            [8.98, 35.412],
            [12.99, 38],
            [17.7, 41.0],
        ]),
    "lvm.skyw.foc":
        np.array([
            [8.98, 35.412],
            [12.99, 38],
            [17.7, 41.0],
        ]),
    "lvm.spec.foc":
        np.array([
            [8.98, 37.412],
            [12.99, 40],
            [17.7, 43.0],
        ]),
}

def temp2focus(telescope:str, temperature:float,  polynomial:float = 2):
    """Gathering focus based on temperature"""

    # get the polynomial coefficients
    d = temp_data[telescope].T
    foc_est = np.polyfit(d[0], d[1], polynomial)

    return np.polyval(foc_est, temperature)
