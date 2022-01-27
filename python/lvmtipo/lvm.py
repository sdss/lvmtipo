# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: lvm.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import json
from types import SimpleNamespace

lvm_tree = '{ "SCI":  \
              { "FOC":   "lvm.sci.foc", \
                "KM":    "lvm.sci.km", \
                "PWI":   "lvm.sci.pwi" \
              }, \
             "SKYW":  \
              { "FOC":   "lvm.skyw.foc", \
                "KM":    "lvm.skyw.km", \
                "PWI":   "lvm.skyw.pwi" \
              }, \
             "SKYE":  \
              { "FOC":   "lvm.skye.foc", \
                "KM":    "lvm.skye.km", \
                "PWI":   "lvm.skye.pwi" \
              }, \
             "SPEC":  \
              { "FOC":   "lvm.spec.foc", \
                "FISEL": "lvm.spec.fibsel", \
                "PWI":   "lvm.spec.pwi" \
              } \
            }'

LVM = json.loads(lvm_tree, object_hook=lambda d: SimpleNamespace(**d))


