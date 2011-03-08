#!/usr/bin/env python

# ADF.py - Class files to build Kohn-Sham data structures from ADF TAPE15 and 
#          TAPE21 output files
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
# Classes included:
# - SFO - essentially the atomic orbital, contained by Frag
# - Frag - information associated with each fragment file including an
# ... array of SFOs
# - KSOrbital - data associated with the KS molecular orbitals
# ... including the energy, occupation and expansion coefficients
# - irrep
# <======================================================================>
#
# Convention, index is zero based, ID is one based
#
# <======================================================================>

import kf      # key file is an ADF supplied interface to dmpkf utility

debugFlagADF = False

AU2eV = 27.2114   # 1 AU (hartree) is 27.211e eV

# <======================================================================>
# Class SFO - Symmetrized Fragment Orbital (really, the atomic orbital)
# Fields: 
# - sym: symbolic text name of orbital 
# - norb: number of orbitals associated with this SFO
# - id: unique numeric id
# Methods: 
# - dump: print out contents of the class to stdout
# <======================================================================>
class SFO:
    """
    Symmetrized Fragment Orbital
    """

    # ------------------------------------------------------------
    def __init__(self, sym, norb):
        self.sym = sym   # symbolic name like F:x
        self.norb = norb # number of orbs assoc with this SFO
        self.id = 0

        if debugFlagADF:
            self.dump()

    # ------------------------------------------------------------
    def dump(self):
        print 'SFO ===> ', self.id, self.sym, self.norb

# <======================================================================>
# Class Frag - Fragment file info, includes SFOs of this fragment
# Fields:
# - atomtype: text name of atomic element
# - sfoList: a python list of SFO's
# Methods:
# - dump: print out contents of the class to stdout
# <======================================================================>
class Frag:
    """
    Fragment
    """

    # ------------------------------------------------------------
    def __init__(self, atomtype, symlab, norb):

        self.atomtype = atomtype  # text name of atom
        self.sfoList = []
        for sym, lnorb in zip(symlab, norb):
            self.sfoList.append(SFO(sym, lnorb))

        if debugFlagADF :
            self.dump()

    # ------------------------------------------------------------
    def dump(self):
        print '<================= Frag  ===========================>'
        print self.atomtype
        for sfo in self.sfoList:
            sfo.dump()
        print '<===================================================>'

        

# <======================================================================>
# Class Irrep - Info generic to all KS orbitals of this symmetry
# ... really, just the SFO coefficents 
# Fields:
# These names are odd, but they match the names of the fields in the 
# ADF key files:
# - SFO_proto: overall SFO coefficient 
# - nofrs: the SFO index within the fragment
# - nofr:  The fragment number associated with the SFO
# - nrcoef: number of coeff associated for symmetry
# Methods:
# - dump: print out contents of the class to stdout
# <======================================================================>
class Irrep:
    """
    Irreducable Representation
    """
    
    # ------------------------------------------------------------
    def __init__(self, irrepName, SFO_proto, nofr, nofrs, nrcoef):
        self.irrepName = irrepName  # symmetry name

        # SFO association with this irrep
        # array of info each the same size 
        self.SFO_proto = SFO_proto # the SFO overall coeff.
        self.nofrs = nofrs # the fragment SFO index with the frag
        self.nofr = nofr # the fragment number of the SFO
        self.nrcoef = nrcoef # number of coeff assoc for symmetry

        if debugFlagADF :
            self.dump()

    # ------------------------------------------------------------
    def dump(self):
        print '<================= Irrep  ===============================>'
        print self.irrepName
        print self.SFO_proto
        print self.nofr
        print self.nofrs
        print self.nrcoef


# <======================================================================>
# Class KSOrbital - Kohn-Sham Molecular Orbitals
# Fields:
# - irrepID: ID of irrep
# - instance: unique ID
# - spin: A or B
# - escale: energy eigenvalue in eV
# - eigenLow: 
# - occ: occupation of this orbital
# Methods: 
# - dump: print out contents of the class to stdout
# - link: a debug routine that links all of the structures so we can verify
#         the extracted data against the data in the ADF logfile.
# <======================================================================> 
class KSOrbital:

    # ------------------------------------------------------------
    def __init__(self, irrepID, instance, spin, escale, eigenLow, occ):
        self.irrepID = irrepID
        self.instance = instance
        self.spin = spin
        self.escale = escale*AU2eV
        self.eigenLow = eigenLow
        self.occ = occ
        if debugFlagADF :
            self.dump()
        
    # ------------------------------------------------------------
    def dump(self):
        print '<=================KS Orbital ===========================>'
        print self.irrepID, '.', self.instance, '.', self.spin
        print 'Energy (eV): ', self.escale
        print 'Occupancy ', self.occ
        
    # ------------------------------------------------------------
    # Basically a sanity check that the code is working - check against logfile
    # ------------------------------------------------------------
    def link(self, irrepList, fragList, sfoGlobalID):

        irrep = irrepList[self.irrepID - 1]
        
        print '<=================KS Orbital ===========================>'
        print irrep.irrepName, '.', self.instance, '.', self.spin
        print 'Energy (eV): ', self.escale
        print 'Occupancy ', self.occ

        coeffIndex = 0
        nrcoefIndex = 0
        count = 0
        
        for fragID, sfoFragID, SFO_proto in zip(irrep.nofr,
                                                 irrep.nofrs,
                                                 irrep.SFO_proto):

            # ADF doesnt print out each coefficient if they are the same
            # move to next coefficients or reset back to reuse previous ones
            if count == 0:
                count = irrep.nrcoef[nrcoefIndex]
                nrcoefIndex = nrcoefIndex + 1
            else: 
                coeffIndex = saveCoeffIndex

            count = count - 1

            # ADF indexes each fragment file ...
            # but it doesn't index each SFO
            f = fragList[fragID - 1]
            sfo = f.sfoList[sfoFragID - 1]

            # set the global SFO ID
            if sfo.id == 0:
                incrementIDFlag = True
                sfo.id = sfoGlobalID
            else:
                incrementIDFlag = False

            saveCoeffIndex = coeffIndex
            # norb - multiple SFOs per SFO, as it were
            for norb in range(sfo.norb):
                coeff = SFO_proto * self.eigenLow[coeffIndex]
                print 'SFO :',sfo.id,':', f.atomtype, '.', fragID, '[', sfo.sym, norb+1, ']',SFO_proto, '*', self.eigenLow[coeffIndex], '=', coeff
                coeffIndex = coeffIndex + 1
                if incrementIDFlag:
                    sfoGlobalID = sfoGlobalID + 2  # 2 IDs, one for each spin
                
        print  '<=======================================================>'
        return sfoGlobalID
        

