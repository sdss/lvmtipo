# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: __init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from sdsstools import get_config, get_logger, get_package_version

# pip package name
NAME = "sdss-lvmtipo"

# Loads config. config name is the package name.
# config = get_config('lvmtan')

# Inits the logging system as NAME. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(path).  Filename can be different
# than NAME.
log = get_logger(NAME)


# package name should be pip package name
__version__ = get_package_version(path=__file__, package_name=NAME)


