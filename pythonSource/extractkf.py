#!/usr/bin/env python

# extractkf.py - Implementation of a Python interface to extract data
#         from an ADF TAPE21 file.  This is an implementation that
#         uses the KF utilities to do the actual work, so you need a
#         working ADF installation. 
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

# Executive Brief: Parse the key files generated from the ADF output files
# to extract information about the KS orbitals, their energy levels, and their 
# expansion into atomic (Lowdin in this case) orbitals (the coefficients).
# That said, the ADF terminology and key file structure is poorly documented and 
# hard to understand ... it takes a lot of comparing between the log file and the key
# file to fully understand.
#
# <======================================================================>
# Extract the KS MOs from the key files (TAPE21 and TAPE15) output
#         from ADF calculations.  Each MOs is associated with the
#         orbital energy and electronic occupancy. Each MO is expanded
#         in terms of Lowdin orbitals associated with a unique SFO.
#         Extract the MO to SFO associations and any overall
#         coefficients due to symmetry (like identical atoms).
# <======================================================================>
# More Detail:  This is confusing, no doubt, so read carefully.  
# Also, try to read the ADF user guide, especially the key file part.
#
# Each KS MO is an eigenfunction of the density functional
#         calculation. Each KS MO is associated with a molecular
#         symmetry irrep which can be expanded into STO basis
#         functions, or, more conveniently, into atomic orbitals.  
# SUMMARY: KS MO (irrep name, id, energy eigenvalue, expansion
#         coeffs)
#
# Each KS MO irrep, regardless of eigenfunction, is generically
#         associated with the underlying orbitals. For example, there
#         may be five eigenfunctions associated with the SIGMA.g
#         irrep, but all of the SIGMA.g functions are associated with
#         exactly the same orbitals (e.g. H P:x, etc.).  The
#         association is stored with each irrep (SYM), rather than
#         with each MO. 
# Additionally, if the underlying orbital is symmetric in some way
#         (like the oxygen orbitals in linear molecule PuO2), the KS
#         MO stores only one coefficient, and the irrep (SYM) stores
#         the overall coefficient and the mapping to the underlying
#         orbitals.  For example, the KS MO might be expanded with
#         coefficients a and b. b might be associated with 2 oxygen
#         orbitals P:x and this assoc is stored with the irrep along
#         with the overall 1/sqrt(2) extra coeff (SFO_proto). 
# SUMMARY: SYM (irrep) (irrep name, underlying orbitals and their
#         overall coefficients)
#
# Each SFO is associated with a particular fragment file. In our
#         studies, each fragment is associated with one atom, and
#         therefore, the SFOs are just the atomic orbitals. This would
#         not be true if the fragment file was associated with a
#         molecule. Further, these SFOs are orthonormalized through a
#         Lowdin orthogonalization procedure.  We want the Lowdin
#         coefficients and ignore the fact that there is admixing
#         between other atomic orbitals. 
# <======================================================================>

import kfwrapper
import kf
import pdb
import numpy as np
import ADF
import cPickle as pickle

import Plist

# print np.__version__

import glob

# <======================================================================>
# SFOPropertyList: a mapping between ADF and our internal format
# Fields:
# - idList - list of valid SFOs 
# - idLastIndex - place holder into list
# Methods:
# - validSFO: compares ADF SFO to one of the atomic orbitals in the atomic
# .. orb property list
# <======================================================================>
class SFOPropertyList:
    """
    Class SFOPropertyList - A class to map ADF SFO's to our internal 
    atomic orbitals based on the property list data.
    """
    # ------------------------------------------------------------
    def __init__(self):
        self.idList = []
        self.idLastIndex = -1

    # <======================================================================>
    # validSFO - check to see if this SFO is in the allowed list
    # .... return a unique ID associated with this SFO
    # .... returns 0 if SFO is not valid
    # <======================================================================>
    def validSFO(self, aoList, f, fragID, sfo, norb, spin):

        #print 'Validating', f.atomtype, fragID, sfo.sym, norb, spin

        retList = [-1,-1]
        
        # search through the atomic orbitals listed in the property list
        # if matched, save the index, 0  means not in property lsit
        if spin == 'A':
            s = 0.5
        else:
            s = -0.5
            
        ao = aoList.validSFO(f.atomtype, fragID, sfo.sym, norb, spin)
        if (ao is not None):
            retList[0] = ao.aoIndex

        # search through the global id list and assign a unique global ID
        for aoIndex, aoTuple in enumerate(self.idList):

            if ((aoTuple[Plist.AoEnum.orbital] == sfo.sym) and
                (aoTuple[Plist.AoEnum.atom] == f.atomtype) and 
                (aoTuple[Plist.AoEnum.fragNo] == fragID) and
                (aoTuple[Plist.AoEnum.orbitalID] == norb) and
                (aoTuple[Plist.AoEnum.spin] == spin)):
                
                retList[1] = aoIndex
                return retList
                
        # add this orbital
        self.idList.append((f.atomtype, fragID, sfo.sym, norb, spin))
        self.idLastIndex = self.idLastIndex + 1
        retList[1] =  self.idLastIndex
        return retList

    
# ------------------------------------------------------------
# routine to extract key file data from *.t15 and *.t21 (TAPE15 
# and TAPE21) files ...
# ...should be one of each in the directory!
# For example, rename the TAPE21 file to projName.t21, etc.
# Extracted data is stored in pickle files:
# - Uo, Uks, Occ, Hdft, Nks
# Uses ADF supplied kf.py and requires valid ADF license (check!) 
# ... the error may be hard to see if invalid license.
# ------------------------------------------------------------
def KFextract():

    kf15Path = glob.glob('*.t15')
    print kf15Path
    kf21Path = glob.glob('*.t21')
    print kf21Path
    if len(kf15Path) != 1 or len(kf21Path) != 1:
        print 'One and only one of t15 and t21 files allowed in directory!'
        return
    
    # <======================================================================>
    # DebugFlag Hook:
    # Use the kfwrapper which dumps everything or you can use kffile
    # silently. Don't change kf cause that is provided (and supported) by
    # ADF directly:
    debugFlag = True
    if debugFlag :
        kf15 = kfwrapper.kffileWrapper(kf15Path[0])
        kf21 = kfwrapper.kffileWrapper(kf21Path[0])
    else:
        kf15 = kf.kffile(kf15Path[0])
        kf21 = kf.kffile(kf21Path[0])
        
    # <======================================================================>
        
    Symmetry_symlab = kf15.read ('Symmetry', 'bb')  # Read the names of the KS orbitals
    print Symmetry_symlab

    moList = []  # Will contain an array of MOs
    irrepList = [] # Will contain an array of info associated with the symmetry
        
    # <======================================================================>
    # Extract the data from the key files associated with each KS orbital
    #
    # Eash KS Symmetry (irrep includes:
    # - symmetry name (irrepName)
    # - the expansion coefficients for duplicate atoms (e.g. 1/sqrt(2))
    # .... (SFO_proto)
    # - list of bas functions that support the symmetry (npart)
    #
    # Each KS Orbital (mo) includes:
    # - unique symmetry index associated with the sym array
    # - symmetry name (irrepName)
    # - instance ID within symmetry name and spin
    # - spin A or B 
    # - energy in eV (escale .. converted from AU in key file) 
    # - coefficients in terms of Lowdins (Eigen-Low)
    # - occupancy (froc)
    #
    # <======================================================================>
    
    irrepID = 0
        
    for irrepName in Symmetry_symlab:
            
        irrepID = irrepID + 1
        
        # --------------------------------------------------------------------
        # First extract the bulk data
        # --------------------------------------------------------------------
        
        # --------------------------------------------------------------------
        # Array of SFOs associated with this irrep
        SFO_proto = kf15.read(irrepName, 'SFO_proto') # SFO expansion coeff
        SFO_nofr = kf15.read(irrepName, 'nofr') # fragment indices
        SFO_nofrs = kf15.read(irrepName, 'nofrs') # SFO indices within fragment
        SFO_nrcoef = kf15.read(irrepName, 'nrcoef') # number of coeff for symmetry
        irrepList.append(ADF.Irrep(irrepName, SFO_proto, SFO_nofr, SFO_nofrs, SFO_nrcoef))  # <==== irrep 
        # --------------------------------------------------------------------
            
        # The expansion of this KS MO in Lowdin orbitals
        EigenLow_A = kf15.read(irrepName, 'Eigen-Low_A') 
        EigenLow_B = kf15.read(irrepName, 'Eigen-Low_B')  
            
        # The number of MOs associated with this symmetry
        nmo_A = kf15.read(irrepName, 'nmo_A')
        nmo_B = kf15.read(irrepName, 'nmo_B')
            
        # The energy eigenvalue associated with this MO
        escale_A = kf21.read (irrepName, 'escale_A')
        escale_B = kf21.read (irrepName, 'escale_B')
            
        # The electronic occupancy of this MO
        froc_A = kf21.read(irrepName,'froc_A')
        froc_B = kf21.read(irrepName,'froc_B')
            
        # --------------------------------------------------------------------
        # Now, for each KS orbital, build an mo instance by spin type A/B
        # --------------------------------------------------------------------
        moInst = 1
            
        # --------------------------------------------------------------------
        for moIndex in range(nmo_A):
            strt = moIndex*nmo_A
            end = strt + nmo_A

            moList.append(ADF.KSOrbital(irrepID, moInst,  # <==== KS MO
                                    'A', 
                                    escale_A[moIndex], 
                                    EigenLow_A[strt:end],
                                    froc_A[moIndex]))  
            moInst = moInst + 1
        # --------------------------------------------------------------------
                
        # --------------------------------------------------------------------
        moInst = 1
        for moIndex in range(nmo_B):
            strt = moIndex*nmo_B
            end = strt + nmo_B
            moList.append(ADF.KSOrbital(irrepID, moInst, # <=============
                                    'B', 
                                    escale_B[moIndex], 
                                    EigenLow_B[strt:end],
                                    froc_B[moIndex]))  
            moInst = moInst + 1
        # --------------------------------------------------------------------
                    
                    
    # <======================================================================>
    # Extract the information about the SFOs from the fragment files
    # 
    # Each sfo includes:
    # - the sfo symmetry name e.g. P:x (ftyp_symlab)
    # - the sfo atom type (e.g. H, O, Pu) (ftyp_atomtyp)
    # - list of bas functions that support the symmetry (npart) 
    # ... linked to the same in the sym array above
    # <======================================================================>
            
    # --------------------------------------------------------------------
    # read the number of fragment types .. really the number of different
    # atoms since each atom has a unique fragment file in our implementation
    # --------------------------------------------------------------------
    ntyp = kf21.read('Geometry', 'ntyp') 
                    
    # --------------------------------------------------------------------
    fragList = []
    fragID = 0
                    
    # Loop through each fragment type
    # Each fragment type is associated with one atom type
    # If that changes, then this extraction code probably won't work
                    
    for fragTypeIndex in range(ntyp):
        fragID = fragID + 1
                        
        # Probably works as long as the number of frags is less than 10
        ftyp_str = 'Ftyp '+str(fragTypeIndex + 1)  
        
        # Atomic Label Pu,H,O,etc.
        ftyp_atomtype =  kf21.read(ftyp_str, 'atomtype') 
                        
        # Symbolic label - the SFO label like P:x, D:xy etc.
        ftyp_symlab = kf21.read(ftyp_str, 'bb')
                        
        # The number of orbitals for each SFO
        # ADF says size of the fock matrix but doesn't make sense
        ftyp_norb = kf21.read(ftyp_str, 'norb')
                        
        ftyp_nfrag = kf21.read(ftyp_str, 'nfrag')
        # Force a fragment instance for each atom
        # eventhough they are the same type
        for n in range(ftyp_nfrag):
            fragList.append(ADF.Frag(ftyp_atomtype[0], ftyp_symlab, ftyp_norb))
                    
    # --------------------------------------------------------------------
    sfoGlobalID = 1
    for ks in moList:
        sfoGlobalID = ks.link(irrepList, fragList, sfoGlobalID)
    # --------------------------------------------------------------------

    for f in fragList:
        f.dump()
    
    # <======================================================================>
    # Read in the property list for the appropriate SFO's to keep
    # <======================================================================>
    pl = Plist.AtomicOrbitalPlist()
    pl.readPlist()
    aoList = pl.aoList
    nStates = aoList.nStates
    
    sfoPlist = SFOPropertyList()

    # <======================================================================>
    nKSOrbitals = len(moList)
    print 'Number of KS Orbitals = ', nKSOrbitals
    Hdft = np.mat(np.zeros((nKSOrbitals, nKSOrbitals), float))
    moIndex = 0
    Uo =  np.mat(np.zeros((nStates, nKSOrbitals), float))
    Uks =  np.mat(np.zeros((nKSOrbitals, nKSOrbitals), float))
    Occ = np.zeros((nKSOrbitals), float)

    for mo in moList:

        # <====================== Energy Matrix =================================>
        irrep = irrepList[mo.irrepID - 1]
        Hdft[moIndex,moIndex] = mo.escale
        print 'E(', irrep.irrepName, '.', mo.instance,'.', mo.spin, ') = ', mo.escale, '(', mo.occ, ')'

        # <====================== Occupation List =================================>
        Occ[moIndex] = mo.occ

        # <======================================================================>
        moIndex = moIndex + 1
        # End for mo in moList



    moIndex = 0
    row = 0
    for moIndex, mo in enumerate(moList):
        coeffIndex = 0
        nrcoefIndex = 0
        count = 0
        irrep = irrepList[mo.irrepID - 1]

        # <======================================================================>
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

            saveCoeffIndex = coeffIndex

            # <======================================================================>
            # Loop thru the sfo's
            # norb - multiple SFOs per SFO, as it were
            # <======================================================================>
            for norb in range(sfo.norb):
                # Generate a unique ID associated with the atomic orbital
                [rowO, rowKS]  = sfoPlist.validSFO(aoList, f, fragID, sfo, norb+1, mo.spin)
                #print rowO, rowKS, f.atomtype, fragID, sfo.sym, norb+1, mo.spin
                if (rowO >= 0):
                    Uo[rowO,moIndex] = SFO_proto*mo.eigenLow[coeffIndex]
                    #print 'Uo[', rowO,'(', f.atomtype, fragID, sfo.sym, norb+1,'),',moIndex,'] = ', Uo[rowO,moIndex]

                if (rowKS == -1):
                    raise(KohnShamError('rowKS is not valid!'))

                Uks[rowKS,moIndex] =  SFO_proto*mo.eigenLow[coeffIndex]
                #print 'Uks[', rowKS,'(', f.atomtype, fragID, sfo.sym, norb+1,'),',moIndex,'] = ', Uks[rowKS,moIndex]

                coeffIndex = coeffIndex + 1

            # <======================================================================>
            # end of sfo loop - looping thru the atomic orbitals
            # <======================================================================>

        # <======================================================================>
        # end of frag loop 
        # <======================================================================>

    # <======================================================================>
    # end of mo loop ... looping through the KS orbitals
    # <======================================================================>

    Nks = np.multiply(np.resize(Occ, Uo.shape), np.multiply(Uo,Uo.conjugate())).sum()
    print 'Nks = ', Nks

    f = open('Uo.pkl', 'w')
    pickle.dump(Uo, f)
    f.close()
    f = open('Uks.pkl', 'w')
    pickle.dump(Uks, f)
    f.close()
    f = open('Occ.pkl', 'w')
    pickle.dump(Occ, f)
    f.close()
    f = open('Hdft.pkl', 'w')
    pickle.dump(Hdft, f)
    f.close()
    f = open('Nks.pkl', 'w')
    pickle.dump(Nks, f)
    f.close()
    
    #f = open('irrep.pkl', 'w')
    #pickle.dump(irrep, f)
    #f.close()
                                
    #f = open('mo.pkl', 'w')
    #pickle.dump(mo, f)
    #f.close()
                                
    #f = open('frag.pkl', 'w')
    #pickle.dump(frag, f)
    #f.close()

    return


                                
