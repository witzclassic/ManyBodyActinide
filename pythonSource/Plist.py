#!/usr/bin/env python

# Plist.py
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
# 
# Propety list processing including classes:
# => PropertyList - wrapper for plistlib routines
# => AoEnum - enumeration for the atomic orbitals
# => AoBasisType - cart or spherical enum
# => SoCouplingEnum - enum for spin-orbit
# => AtomicOrbitalPlist - a property list that includes atomic orbitals (all in our case)
# => atomic specific property lists which saves me copying it into every project directory
#
# <======================================================================>
import plistlib
import AO

# <======================================================================>
# PropertyList base class
# Uses plistlib ... read the python docs
# Fields:
# - fname: file name including any suffix
# - pl: the python property list 
# Methods:
# - readPlist
# - setPlist
# - writePlist
# <======================================================================>
class PropertyList:
    # ------------------------------------------------------------
    def __init__(self,plistFname='Ham.plist'):
        self.fname = plistFname
    
    # ------------------------------------------------------------
    def readPlist(self):
        self.pl = plistlib.readPlist(self.fname)
        
    # ------------------------------------------------------------
    def setPlist(self,pl):
        self.pl = pl
        
    # ------------------------------------------------------------
    def writePlist(self):
        plistlib.writePlist(self.pl, self.fname)

# <======================================================================>
# Enumeration for the atomic orbital list
# 'H', 1, 'P:x', 1, 'A'
# Means Hydrogen, first atom, orbital (n=2) P:x, spin up
# The n=2 comes from the orbitalNo which is the first n that
# that a p orbital can be in.  For example ...
# 'H', 1, 'S', 1, 'A" means orbital (n=1) S
# It is just because ADF does it this way.
# 
# atom - text name of the atom, must match ADF fragment name
# fragNo - the ADF fragNo .. really the atom number indexed within the atom name
# orbital - the text name of the ADF SFO
# orbitalNo - the ADF index for n .. see above
# spin - 'A' or 'B'
# l - angular momentum
# lz - if cartesian basis, >0 is the real component, <0 the imaginary
# ... if spherical, lz azimuth value
# <======================================================================>

class AoEnum:
    atom, fragNo, orbital, orbitalID, spin, l, lz = range(7)

# Enumeration to tell what basis set we're using
class AoBasisType:
    cartesian, spherical,j = range(3)

# Enumeration for coupling constant tuplets
# l that this applies to, Xi and Alpha constants
class SoCouplingEnum:
    l, soXi, soAlpha = range(3)

# <======================================================================>
# Atomic Orbital Property List
# .. interprets the property list file as a atomic orbital property list
# Inherits Property List
# Fields: 
# Methods:
# - readPlist ... read in the file
# The other routines interpret the pl data and return the values:
# - fragNo
# - atom
# - orbsAreEqual
# - spinList
# - basisSetType
# - slaterCondonParameters
# - soCouplingConstants
# - SOalpha
# <======================================================================>
class AtomicOrbitalPlist(PropertyList):

    # ------------------------------------------------------------
    def __init__(self, fname='Ham.plist'):
        PropertyList.__init__(self, fname)

    # ------------------------------------------------------------
    def readPlist(self):
        PropertyList.readPlist(self)
        self.aoList = AO.AtomicOrbitalList(self)

    # ------------------------------------------------------------
    def fragNo(self, index):
        aoList = self.pl['aoList']
        ao = aoList[index]
        return ao[AoEnum.fragNo] 
            
    # ------------------------------------------------------------
    def atom(self, index):
        aoList = self.pl['aoList']
        ao = aoList[index]
        return ao[AoEnum.atom]

    # ------------------------------------------------------------
    def orbsAreEqual(self, index1, index2):
        aoList = self.pl['aoList']
        ao1 = aoList[index1]
        ao2 = aoList[index2]

        if (ao1[AoEnum.atom] == ao2[AoEnum.atom]
            and ao1[AoEnum.fragNo] == ao2[AoEnum.fragNo]
            and ao1[AoEnum.orbital] == ao2[AoEnum.orbital]
            and ao1[AoEnum.orbitalID] == ao2[AoEnum.orbitalID]):

            return True

        else:
            return False

    # ------------------------------------------------------------
    def spinList(self):
        spinList = []
        aoList = self.pl['aoList']
        for ao in aoList:
            if ao[AoEnum.spin] == 'A':
                spinList.append(1)
            else:
                spinList.append(0)
        return spinList
    
    # ------------------------------------------------------------
    def basisSetType(self):
        return self.pl['aoBasisType']

    # ------------------------------------------------------------
    def slaterCondonParameters(self):
        return self.pl['slaterCondonParameters']

    # ------------------------------------------------------------
    # SO Coupling constants 
    def soCouplingConstants(self):
        return self.pl['soCouplingConstants'];

    # ------------------------------------------------------------
    def SOalpha(self):
        return self.pl['SOalpha']

# class AnO2Plist(AtomicOrbitalPlist):

#     def __init__(self, An, scl, soCC, fname='Ham.plist'):

#         AtomicOrbitalPlist.__init__(self, fname)

#         self.pl = dict(
#             slaterCondonParameters=scl,
#             soCouplingConstants=soCC,
#             aoBasisType=AoBasisType.cartesian,
#             # The order of this is important for performance reasons
#             aoList = [
#                 ('O', 2, 'P:x', 1, 'A', 1, 1),
#                 ('O', 2, 'P:y', 1, 'A', 1, -1),
#                 ('O', 2, 'P:z', 1, 'A', 1, 0),
#                 ('O', 3, 'P:x', 1, 'A', 1, 1),
#                 ('O', 3, 'P:y', 1, 'A', 1, -1),
#                 ('O', 3, 'P:z', 1, 'A', 1, 0),
#                 ('O', 2, 'P:x', 1, 'B', 1, 1),
#                 ('O', 2, 'P:y', 1, 'B', 1, -1),
#                 ('O', 2, 'P:z', 1, 'B', 1, 0),
#                 ('O', 3, 'P:x', 1, 'B', 1, 1),
#                 ('O', 3, 'P:y', 1, 'B', 1, -1),
#                 ('O', 3, 'P:z', 1, 'B', 1, 0),
#                 (An, 1, 'F:z3', 1, 'A', 3, 0),
#                 (An, 1, 'F:x', 1, 'A', 3, 1),
#                 (An, 1, 'F:y', 1, 'A', 3, -1),
#                 (An, 1, 'F:z', 1, 'A', 3, 2),
#                 (An, 1, 'F:xyz', 1, 'A', 3, -2),
#                 (An, 1, 'F:z2x', 1, 'A', 3, 3),
#                 (An, 1, 'F:z2y', 1, 'A', 3, -3),
#                 (An, 1, 'F:z3', 1, 'B', 3, 0),
#                 (An, 1, 'F:x', 1, 'B', 3, 1),
#                 (An, 1, 'F:y', 1, 'B', 3, -1),
#                 (An, 1, 'F:z', 1, 'B', 3, 2),
#                 (An, 1, 'F:xyz', 1, 'B', 3, -2),
#                 (An, 1, 'F:z2x', 1, 'B', 3, 3),
#                 (An, 1, 'F:z2y', 1, 'B', 3, -3)
#                 ]
#             )
#         self.writePlist()
        

# class AnPlist(AtomicOrbitalPlist):
#     def __init__(self, An, scl, soCC, aoBT=AoBasisType.cartesian, fname='Ham.plist'):
# #    def __init__(self, An, scl, soCC, aoBT=AoBasisType.spherical, fname='Ham.plist'):

#         AtomicOrbitalPlist.__init__(self, fname)
        
#         print 'Building Plist for ', An
#         print 'Slater Condon Parameters:'
#         for sc in scl:
#             print scl
#         print 'Spin Orbit'
#         for cc in soCC:
#             print cc

#         self.pl = dict(
#             soCouplingConstants=soCC,
#             slaterCondonParameters=scl,
#             aoBasisType=aoBT,
#             aoList = [
#                 (An, 1, 'F:z3', 1, 'A', 3, 0), # 0
#                 (An, 1, 'F:x', 1, 'A', 3, 1),  # 1
#                 (An, 1, 'F:y', 1, 'A', 3, -1), # 2
#                 (An, 1, 'F:z', 1, 'A', 3, 2),  # 3
#                 (An, 1, 'F:xyz', 1, 'A', 3, -2), # 4
#                 (An, 1, 'F:z2x', 1, 'A', 3, 3), # 5
#                 (An, 1, 'F:z2y', 1, 'A', 3, -3), # 6
#                 (An, 1, 'F:z3', 1, 'B', 3, 0), # 7
#                 (An, 1, 'F:x', 1, 'B', 3, 1), # 8
#                 (An, 1, 'F:y', 1, 'B', 3, -1), # 9
#                 (An, 1, 'F:z', 1, 'B', 3, 2), # 10
#                 (An, 1, 'F:xyz', 1, 'B', 3, -2), # 11
#                 (An, 1, 'F:z2x', 1, 'B', 3, 3), # 12
#                 (An, 1, 'F:z2y', 1, 'B', 3, -3) # 13
#                 ]
#             )
#         self.writePlist()
    

# <======================================================================>
# Actinide Property List
# .. interprets the property list file as a actinide orbital property list
# 
# Inherits Atomic Orbital Property List
# Used for An or AnO_2 type molecules
# Fields: 
# Methods:
#
# <======================================================================>
class AnPlist(AtomicOrbitalPlist):
    """
    AnPlist - Actinide molecule property list
    """

    # ------------------------------------------------------------
    def __init__(self, An, scl, soCC, withO=False, withD=False, fname='Ham.plist', basisType=AoBasisType.cartesian):

        AtomicOrbitalPlist.__init__(self, fname)

        print withO, withD

        if withD and withO:
            aol = [
p                ('O', 2, 'P:x', 1, 'A', 1, 1),
                ('O', 2, 'P:y', 1, 'A', 1, -1),
                ('O', 2, 'P:z', 1, 'A', 1, 0),
                ('O', 3, 'P:x', 1, 'A', 1, 1),
                ('O', 3, 'P:y', 1, 'A', 1, -1),
                ('O', 3, 'P:z', 1, 'A', 1, 0),
                ('O', 2, 'P:x', 1, 'B', 1, 1),
                ('O', 2, 'P:y', 1, 'B', 1, -1),
                ('O', 2, 'P:z', 1, 'B', 1, 0),
                ('O', 3, 'P:x', 1, 'B', 1, 1),
                ('O', 3, 'P:y', 1, 'B', 1, -1),
                ('O', 3, 'P:z', 1, 'B', 1, 0),
                (An, 1, 'F:z3', 1, 'A', 3, 0),
                (An, 1, 'F:x', 1, 'A', 3, 1),
                (An, 1, 'F:y', 1, 'A', 3, -1),
                (An, 1, 'F:z', 1, 'A', 3, 2),
                (An, 1, 'F:xyz', 1, 'A', 3, -2),
                (An, 1, 'F:z2x', 1, 'A', 3, 3),
                (An, 1, 'F:z2y', 1, 'A', 3, -3),
                (An, 1, 'F:z3', 1, 'B', 3, 0),
                (An, 1, 'F:x', 1, 'B', 3, 1),
                (An, 1, 'F:y', 1, 'B', 3, -1),
                (An, 1, 'F:z', 1, 'B', 3, 2),
                (An, 1, 'F:xyz', 1, 'B', 3, -2),
                (An, 1, 'F:z2x', 1, 'B', 3, 3),
                (An, 1, 'F:z2y', 1, 'B', 3, -3),
                (An, 1, 'D:z2', 1, 'A', 2, 0),
                (An, 1, 'D:xz', 1, 'A', 2, 1),
                (An, 1, 'D:yz', 1, 'A', 2, -1),
                (An, 1, 'D:xy', 1, 'A', 2, 2),
                (An, 1, 'D:x2-y2', 1, 'A', 2, -2),
                (An, 1, 'D:z2', 1, 'B', 2, 0),
                (An, 1, 'D:xz', 1, 'B', 2, 1),
                (An, 1, 'D:yz', 1, 'B', 2, -1),
                (An, 1, 'D:xy', 1, 'B', 2, 2),
                (An, 1, 'D:x2-y2', 1, 'B', 2, -2)
                ]
        elif withD:
            aol = [
                (An, 1, 'F:z3', 1, 'A', 3, 0),
                (An, 1, 'F:x', 1, 'A', 3, 1),
                (An, 1, 'F:y', 1, 'A', 3, -1),
                (An, 1, 'F:z', 1, 'A', 3, 2),
                (An, 1, 'F:xyz', 1, 'A', 3, -2),
                (An, 1, 'F:z2x', 1, 'A', 3, 3),
                (An, 1, 'F:z2y', 1, 'A', 3, -3),
                (An, 1, 'F:z3', 1, 'B', 3, 0),
                (An, 1, 'F:x', 1, 'B', 3, 1),
                (An, 1, 'F:y', 1, 'B', 3, -1),
                (An, 1, 'F:z', 1, 'B', 3, 2),
                (An, 1, 'F:xyz', 1, 'B', 3, -2),
                (An, 1, 'F:z2x', 1, 'B', 3, 3),
                (An, 1, 'F:z2y', 1, 'B', 3, -3),
                (An, 1, 'D:z2', 1, 'A', 2, 0),
                (An, 1, 'D:xz', 1, 'A', 2, 1),
                (An, 1, 'D:yz', 1, 'A', 2, -1),
                (An, 1, 'D:xy', 1, 'A', 2, 2),
                (An, 1, 'D:x2-y2', 1, 'A', 2, -2),
                (An, 1, 'D:z2', 1, 'B', 2, 0),
                (An, 1, 'D:xz', 1, 'B', 2, 1),
                (An, 1, 'D:yz', 1, 'B', 2, -1),
                (An, 1, 'D:xy', 1, 'B', 2, 2),
                (An, 1, 'D:x2-y2', 1, 'B', 2, -2)
                ]
        elif withO:
            aol = [
                ('O', 2, 'P:x', 1, 'A', 1, 1),
                ('O', 2, 'P:y', 1, 'A', 1, -1),
                ('O', 2, 'P:z', 1, 'A', 1, 0),
                ('O', 3, 'P:x', 1, 'A', 1, 1),
                ('O', 3, 'P:y', 1, 'A', 1, -1),
                ('O', 3, 'P:z', 1, 'A', 1, 0),
                ('O', 2, 'P:x', 1, 'B', 1, 1),
                ('O', 2, 'P:y', 1, 'B', 1, -1),
                ('O', 2, 'P:z', 1, 'B', 1, 0),
                ('O', 3, 'P:x', 1, 'B', 1, 1),
                ('O', 3, 'P:y', 1, 'B', 1, -1),
                ('O', 3, 'P:z', 1, 'B', 1, 0),
                (An, 1, 'F:z3', 1, 'A', 3, 0),
                (An, 1, 'F:x', 1, 'A', 3, 1),
                (An, 1, 'F:y', 1, 'A', 3, -1),
                (An, 1, 'F:z', 1, 'A', 3, 2),
                (An, 1, 'F:xyz', 1, 'A', 3, -2),
                (An, 1, 'F:z2x', 1, 'A', 3, 3),
                (An, 1, 'F:z2y', 1, 'A', 3, -3),
                (An, 1, 'F:z3', 1, 'B', 3, 0),
                (An, 1, 'F:x', 1, 'B', 3, 1),
                (An, 1, 'F:y', 1, 'B', 3, -1),
                (An, 1, 'F:z', 1, 'B', 3, 2),
                (An, 1, 'F:xyz', 1, 'B', 3, -2),
                (An, 1, 'F:z2x', 1, 'B', 3, 3),
                (An, 1, 'F:z2y', 1, 'B', 3, -3)
                ]
        else:
            aol = [
                (An, 1, 'F:z3', 1, 'A', 3, 0),
                (An, 1, 'F:x', 1, 'A', 3, 1),
                (An, 1, 'F:y', 1, 'A', 3, -1),
                (An, 1, 'F:z', 1, 'A', 3, 2),
                (An, 1, 'F:xyz', 1, 'A', 3, -2),
                (An, 1, 'F:z2x', 1, 'A', 3, 3),
                (An, 1, 'F:z2y', 1, 'A', 3, -3),
                (An, 1, 'F:z3', 1, 'B', 3, 0),
                (An, 1, 'F:x', 1, 'B', 3, 1),
                (An, 1, 'F:y', 1, 'B', 3, -1),
                (An, 1, 'F:z', 1, 'B', 3, 2),
                (An, 1, 'F:xyz', 1, 'B', 3, -2),
                (An, 1, 'F:z2x', 1, 'B', 3, 3),
                (An, 1, 'F:z2y', 1, 'B', 3, -3)
                ]
        self.pl = dict(
            slaterCondonParameters=scl,
            soCouplingConstants=soCC,
            aoBasisType=basisType,
            aoList = aol
            )
        self.writePlist()
                
                
