#!/usr/bin/env python

# analyzeLLS - utility function used for post diag analysis
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
import ManyBodyState as mbs
import numpy as np
import cPickle as pickle

f = open('mbDimension.pkl', 'r')
mbDim = pickle.load(f)
f.close()

f = open('LowLyingStates.mtx', 'r')
header = f.readline()
entries, gsEnergy = header.split()
vList = np.zeros(mbDim, complex)

for index, line in enumerate(f):
    real, imag = line.split()
    vList[index] = complex(float(real), float(imag))

f.close()
v = np.array(vList)
mbs.analyzeOccupations(v)
