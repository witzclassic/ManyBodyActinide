#!/usr/bin/env python

# kfwrapper.py - Debug code for the kf code
#
# Imports: kf
# Requires: dmpkf ADF utility (valid license)
#
# Copyright (C) 2009 Brown University Physics (Prof. J.B. Marston)
# Author: Steve Horowitz
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Brown owns the intellectual property right for this file and reserves the
# right to distribute it under a license other than LGPL
# <======================================================================>

import kf

# <======================================================================>
# A wrapper for the ADF kf routine which dumps everything read from the adf tape 
# ... file to stdout
# <======================================================================>
class kffileWrapper:

    def __init__(self, filename):
        self.kf = kf.kffile(filename)
        print 'Using ', filename

    def read( self, section, variable ):
        result = self.kf.read(section, variable)
        print section, variable, result
        return result
