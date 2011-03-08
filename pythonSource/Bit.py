#!/usr/bin/env python

# <======================================================================>
# Bit - A Bit String Class BitStr
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
# BitStr class used for bit manipulations, not used in this package because
# it is too slow, the bit stuff was moved in line.
# <======================================================================>

class BitStr:
    # Default to a 32 bit string
    # but python is flexible
    def __init__(self, len = 32, value=0x00000000):
        self.value = value
        self.len = len

    def clear(self, pos):
        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr clear) Pos'+str(pos)+'out of range '+str(self.len)

        # --------------------------------------------------------------------

        self.value = self.value & ~(0x1 << pos)

    def set(self, pos):

        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr set) Pos'+str(pos)+'out of range '+str(self.len)
        # --------------------------------------------------------------------

        self.value = self.value | (0x1 << pos)

    def get(self, pos):
        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr get) Pos'+str(pos)+'out of range '+str(self.len)
        if self.value & (0x1 << pos):
            return True
        else:
            return False

    def display(self):
        # Slow, so rewrite if we need performance
        bstr = '0b'
        bOccList = []
        for b in range(self.len,0,-1):
            if self.get(b-1):
                bstr += '1'
                bOccList.append(b-1)
            else:
                bstr += '0'
        return bstr+str(bOccList)

    def __getitem__(self, index):
        return self.get(index)

