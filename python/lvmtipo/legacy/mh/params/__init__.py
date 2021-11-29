# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: params/__init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import importlib

def load(name:str):
    return importlib.import_module(name)
