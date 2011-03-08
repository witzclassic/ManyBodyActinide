#!/usr/bin/env python

# AO - Atomic Orbital Classes
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
# Atomic Orbital Classes used extensively for expanding an orbital basis
# into another orbital basis, or for parsing thru the atomic orbitals.
#
# <======================================================================>
# 
# AngMoExp - supports expansion of one orb basis to another
# AtomicOrbital - base class orbID, l and aoIndex
# JAtomicOrbital - base class for j,jz
# nonJAtomicOrbtial - base class for lz, sz
# SphericalAtomicOrbital - typical hydrogenic basis 
# CartesianAtomicOrbital - e.g. F:xz F:yz
# AtomicOrbitalList - takes care of creating a list of the appropriate atom type
# 
# <======================================================================>

import CG
import Plist
import cmath

debugFlagAO = False

# <======================================================================>
# Class AngMoExp - contains AtomicOrbital (ao) and a coefficient c
# ... allows expansion of one atomic orbital into series of another
# Fields: 
# - ao: atomic orbital
# - c: coefficient
# Methods:
# - dump: print out contents of the class to stdout
# - equal: compares two instances
# <======================================================================>
class AngMoExp():
    """
    Angular Momentum Expansion
    """

    # ------------------------------------------------------------
    def __init__(self, ao, c):
        self.ao = ao
        self.c = c

    # ------------------------------------------------------------
    def dump(self):
        self.ao.dump()
        print self.c

    # ------------------------------------------------------------
    def equal(self, rhs):
        self.dump()
        rhs.dump()
        if self.j != rhs.j or self.jz != rhs.jz:
            return False
        else:
            return True

# <======================================================================>
# Class AtomicOrbital - prototypical atomic orbital quantum numbers
# Virtual Class .... 
# Fields:
# - atom
# - fragNo: index into frag list
# - orb: n quantum ID
# - orbID:
# - l: quantum number
# - aoIndex: 
# <======================================================================>
class AtomicOrbital():
    """
    Atomic Orbital base class
    """

    # ------------------------------------------------------------
    def __init__(self, atom, fragNo, orb, orbID, l, aoIndex):

        self.atom = atom
        self.fragNo = fragNo
        self.orb = orb

        self.orbID = orbID  
        self.l = l

        self.aoIndex = aoIndex


# <======================================================================>
# JAtomicOrbital - atomic orbital in the j, jz basis
# Inherits: Atomic Orbital
# Fields: j and jz atomic numbers
# Methods: 
# - dump: print out contents of the class to stdout
# - jExpGen: a generator iterator
# <======================================================================>
class JAtomicOrbital(AtomicOrbital):
    """
    Atomic Orbital with j, jz quantum numbers
    """

    # ------------------------------------------------------------
    def __init__(self, atom, fragNo, orb, orbID, j, jz, l, aoIndex):
        AtomicOrbital.__init__(self, atom, fragNo, orb, orbID, l, aoIndex)
        self.j = j
        self.jz = jz

    # ------------------------------------------------------------
    def dump(self):
        print 'Index', self.aoIndex, 'n', self.fragNo, self.orb, self.orbID, 'l', self.l, 'j', self.j, 'jz', self.jz
        
    # ------------------------------------------------------------
    def jExpGen(self):
        yield self

# <======================================================================>
# nonJAtomicOrbital 
# Inherits: Atomic Orbital
# Fields: lz and sz atomic numbers
# Methods: 
# - dump: print out contents of the class to stdout
# - _eq_: works for ao1==ao2
# - _ne_: works for ao1!=ao2
# <======================================================================>
class nonJAtomicOrbital(AtomicOrbital):
    """
    Atomic Orbital with lz and sz quantum numbers
    """

    # ------------------------------------------------------------
    def __init__(self, atom, fragNo, orb, orbID, l,lz,sz,aoIndex):
        AtomicOrbital.__init__(self, atom, fragNo, orb, orbID, l, aoIndex)
        self.lz = lz
        self.sz = sz

    # --------------------------------------------------            
    # equal? built in for ao1 == ao2 ?
    # --------------------------------------------------            
    def __eq__(self, rhs):
        if self.orbID != rhs.orbID \
                or self.l != rhs.l \
                or self.lz != rhs.lz \
                or self.sz != rhs.sz:
            return False
        else:
            return True
        
    # --------------------------------------------------            
    # not equal? built in for ao1 != ao2 ?
    # --------------------------------------------------            
    def __ne__(self, rhs):
        return not self.__eq__(rhs)
        
    # --------------------------------------------------            
    # dump the class to the terminal
    # --------------------------------------------------            
    def dump(self):
        print 'Index', self.aoIndex, self.atom, 'n=', self.fragNo, self.orb, self.orbID, 'l=', self.l, 'lz=', self.lz, 'sz=', self.sz
        
# <======================================================================>
# Spherical Atomic Orbitals
#     Typical Hydrogenic atomic orbitals in orbID, l, lz, sz basis
# Inherits nonJAtomicOrbitals
# Fields: None
# Methods:
# - dump: print out contents of the class to stdout
# - sphericalExpGen: generator iterates through spherical atomic orbs
# - switchExpGen: generic generator that iterates in other representation, in this case
#                 the cartesian atomic orbitals.
# - jExpGen: generator expands in j, jz basis
# <======================================================================>
class SphericalAtomicOrbital(nonJAtomicOrbital):
    """
    Spherical ATomic Orbital
    """

    # ------------------------------------------------------------
    def __init__(self, atom, fragNo, orb, orbID, l,lz,sz,aoIndex):
        nonJAtomicOrbital.__init__(self, atom, fragNo, orb, orbID, l,lz,sz,aoIndex)

    # --------------------------------------------------            
    # dump the class to the terminal
    # --------------------------------------------------            
    def dump(self):
        print 'Spherical Atomic Orbital: '
        nonJAtomicOrbital.dump(self)

    # --------------------------------------------------            
    # sphericalExpGen a generator fn to iterate through the expansion
    # of this ao into spherical atomic orbitals ... in this case, it does 
    # nothing.
    # 
    # Example usage: for sao in ao.sphericalExpGen():
    # --------------------------------------------------            
    def sphericalExpGen(self):
        sAOExp = AngMoExp(self, 1.0)
        yield sAOExp
        return

    # --------------------------------------------------            
    # switchExpGen a generator fn to iterate through the expansion
    # of this ao into the opposite basis 
    # ... in this case into the Cartesian Basis
    # 
    # For example, a p orbital is expanded as:
    #
    # lz = 0  => 1.0 P:z
    # lz = 1  => 1/sqrt(2) [P:x + i P:y]
    # lz = -1 => 1/sqrt(2) [P:x - i P:y]
    #
    # Example usage: for xao in ao.switchExpGen():
    # --------------------------------------------------            
    def switchExpGen(self):
        if self.lz == 0:
            cAOExp = AngMoExp(CartesianAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0)
            yield cAOExp
        elif self.lz > 0:
            cAOExp = AngMoExp(CartesianAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield cAOExp
            cAOExp = AngMoExp(CartesianAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), 1.0j/cmath.sqrt(2.0))
            yield cAOExp
        elif self.lz < 0:
            cAOExp = AngMoExp(CartesianAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield cAOExp
            cAOExp = AngMoExp(CartesianAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), -1.0j/cmath.sqrt(2.0))
            yield cAOExp

            
    # --------------------------------------------------            
    # jExpGen a generator fn to iterate through the expansion
    # of this ao into the j,jz basis
    # Example usage: for jao in ao.jExpGen():
    # --------------------------------------------------            
    def jExpGen(self):
        jp = self.l + 0.5
        jm = self.l - 0.5
        jz = self.lz + self.sz

        # Build the expansion list of this atomic orbital in the j/jz basis
        # ClebschGordan (j1, j2, j, jz1, jz2, jz)
        # where |j1 - j2| <= j <= |j1 + j2| and jz = jz1 + jz2 
        coeff = CG.ClebschGordan(self.l, 0.5, jp, self.lz, self.sz, jz)
        # print self.lz, self.sz, jp, jz, coeff
        if abs(coeff) > 0:
            yield AngMoExp(JAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, jp, jz, self.l, self.aoIndex), coeff)
        coeff = CG.ClebschGordan(self.l, 0.5, jm, self.lz, self.sz, jz)
        # print self.lz, self.sz, jm, jz, coeff
        if abs(coeff) > 0:
            yield AngMoExp(JAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, jm, jz, self.l, self.aoIndex), coeff)

# <======================================================================>
# Cartesian Atomic Orbitals
#
# Like P:x, P:y and P:z
# contains a list of spherical orbitals .. 
# 
# Spherical to Cartesian:
# [1,0] = P:z
# [1,1] = 1/sqrt(2) (P.x + i P.y)
# [1,-1] = 1/sqrt(2) (P.x - i P.y)
#
# Cartesian to Spherical
# P:z = [1,0]
# P:x = 1/sqrt(2) ( [1,1] + [1,-1] )
# P:y = -i/sqrt(2) ( [1,1] - [1,-1] )
# Inherits nonJAtomicOrbitals
# Fields: None
# Methods:
# - dump: print out contents of the class to stdout
# - sphericalExpGen: generator iterates through spherical atomic orbs
# - switchExpGen: generic generator that iterates in other representation, in this case
#                 the spherical atomic orbitals.
# - jExpGen: generator expands in j, jz basis
# <======================================================================>
class CartesianAtomicOrbital(nonJAtomicOrbital):
    """
    Cartesian Atomic Orbital
    """

    # --------------------------------------------------            
    def __init__(self, atom, fragNo, orb, orbID, l, lz, sz, aoIndex):
        nonJAtomicOrbital.__init__(self, atom, fragNo, orb, orbID, l, lz, sz, aoIndex)

    # --------------------------------------------------            
    # dump the class to the terminal
    # --------------------------------------------------            
    def dump(self):
        print 'Carteisan Atomic Orbital: '
        nonJAtomicOrbital.dump(self)

    # --------------------------------------------------            
    # switchExpGen a generator fn to iterate through the expansion
    # of this ao into the opposite basis 
    # ... in this case into the Spherical Basis
    # 
    # Example usage: for xao in ao.switchExpGen():
    # --------------------------------------------------            
    def switchExpGen(self):

        if self.lz == 0:
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0)
            yield sAOExp

        # e.g. P.x = 1/sqrt(2) ( [1] + [-1] )
        elif self.lz > 0: 
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield sAOExp
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield sAOExp

        # e.g. P.y = -i/sqrt(2) ( [1] - [-1] )
        else:
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), -1.0j/cmath.sqrt(2.0))
            yield sAOExp
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0j/cmath.sqrt(2.0))
            yield sAOExp

    # --------------------------------------------------            
    # sphericalExpGen a generator fn to iterate through the expansion
    # of this ao into spherical atomic orbitals ... in this case, 
    # the spherical basis.
    #
    # Example usage: for sao in ao.sphericalExpGen():
    # --------------------------------------------------            
    def sphericalExpGen(self):

        if self.lz == 0:
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0)
            yield sAOExp

        # e.g. P.x = 1/sqrt(2) ( [1] + [-1] )
        elif self.lz > 0: 
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield sAOExp
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), 1.0/cmath.sqrt(2.0))
            yield sAOExp

        # e.g. P.y = -i/sqrt(2) ( [1] - [-1] )
        else:
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, -self.lz, self.sz, self.aoIndex), -1.0j/cmath.sqrt(2.0))
            yield sAOExp
            sAOExp = AngMoExp(SphericalAtomicOrbital(self.atom, self.fragNo, self.orb, self.orbID, self.l, self.lz, self.sz, self.aoIndex), 1.0j/cmath.sqrt(2.0))
            yield sAOExp


    def jExpGen(self):
        for sAO in self.sphericalExpGen():
            for jAO in sAO.ao.jExpGen():
                jAO.c = jAO.c * sAO.c
                yield jAO
        return
                
# <======================================================================>
# AtomicOrbitalList - a list of atomic orbitals
# Fields: list: list of atomic orbitals
# Methods: 
# - validSFO: Checks if SFO if one we want
# - get: returns atomic orbital matching requested attributes
# - aoByLGen: generator iterates by quantum id l, returns atomic orbital
# - saoByLGen: generator iterates by quantum id l, returns spherical atomic orbital
# __getitem__: builtin ao = aoList[index]
# __iter__: builtin for ao in aoList:
# reverse: same as iter but in reverse order.
# <======================================================================>
class AtomicOrbitalList:

    def __init__(self, pl):
        aoList = pl.pl['aoList']
        self.list = []

        for aoIndex, aoTuple in enumerate(aoList):

            # n + l because ADF is stupid
            # for example, ADF starts p orbitals at n=1 (not 2), etc.
            # I probably should bury this in the ADF specific code ... TO DO list
            n = aoTuple[Plist.AoEnum.orbitalID] + aoTuple[Plist.AoEnum.l]

            if aoTuple[Plist.AoEnum.spin] == 'A':
                sz = 0.5
            else:
                sz = -0.5
            if (pl.basisSetType() == Plist.AoBasisType.cartesian):
                self.list.append( \
                    CartesianAtomicOrbital( \
                        aoTuple[Plist.AoEnum.atom], \
                            aoTuple[Plist.AoEnum.fragNo], \
                            aoTuple[Plist.AoEnum.orbital], \
                            aoTuple[Plist.AoEnum.orbitalID], \
                            aoTuple[Plist.AoEnum.l], \
                            aoTuple[Plist.AoEnum.lz], \
                            sz, \
                            aoIndex))
            elif (pl.basisSetType() == Plist.AoBasisType.spherical):
                self.list.append( \
                    SphericalAtomicOrbital( \
                        aoTuple[Plist.AoEnum.atom], \
                            aoTuple[Plist.AoEnum.fragNo], \
                            aoTuple[Plist.AoEnum.orbital], \
                            aoTuple[Plist.AoEnum.orbitalID], \
                            aoTuple[Plist.AoEnum.l], \
                            aoTuple[Plist.AoEnum.lz], \
                            sz, \
                            aoIndex))

        self.nStates = aoIndex+1
        
        if debugFlagAO:
            for aoIndex, ao in enumerate(self.list):
                print 'Atom', aoIndex
                ao.dump()

    # --------------------------------------------------
    # get a specific atomic orbital ... valid SFO (ADF)
    # ... same as get below, but does not include redundant info of l,lz
    # --------------------------------------------------
    def validSFO(self, atom, fragNo, orb, orbID, spin):

        if spin == 'A':
            s = 0.5
        else:
            s = -0.5
            
        for ao in self.list:

            if (ao.atom == atom and 
                ao.fragNo == fragNo and
                ao.orb == orb and 
                ao.orbID == orbID and 
                ao.sz == s):

                return ao

        return None

        
    # --------------------------------------------------
    # get a specific atomic orbital
    # --------------------------------------------------
    def get(self, atom, fragNo, orb, orbID, l, lz, sz):

        for ao in self.list:
            if (ao.atom == atom and 
                ao.fragNo == fragNo and
                ao.orb == orb and 
                ao.orbID == orbID and 
                ao.l == l and
                ao.lz == lz and
                ao.sz == sz):

                return ao

        return None

    # --------------------------------------------------
    # usage: for ao in aoList.aoByLGen(desired l)
    # loops through aoList only returning aos that have appropriate l
    # --------------------------------------------------
    def aoByLGen(self, l):
        for ao in self.list:
            if ao.l == l:
                yield ao
        
    def saoByLGen(self, L):
        for ao in self.list:
            if ao.l == l:
                yield ao.sphericalExpGen()

    # --------------------------------------------------
    # usage: ao = aoList[index]
    # --------------------------------------------------
    def __getitem__(self, index):
        return self.list[index]

    # --------------------------------------------------
    # usage: for ao in aoList:
    # --------------------------------------------------
    def __iter__(self):
        return iter(self.list)
    
    # --------------------------------------------------
    # same as iter but in reverse
    # --------------------------------------------------
    def reverse(self):
        return self.list.reverse()
