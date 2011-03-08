#!/usr/bin/env python

# Coulomb
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
# Coulomb processing
# ... contains the following classes...
# => GauntCoefficient, and associated orbtial specific classes
# - GauntPP, GauntFF
# - Coulomb Interaction - the I interaction tensor indexed by lz
# - CoulombTypeEnum - Coulomb or Mean Field
# - GenCoulombOps - Build the two-particle operator(s)
# => CoulombOperator (TwoParticleOperator) and assoc. orbital specific classes
# => MeanFieldOperator (OneParticleOperator) and assoc. orbital specific classes
# => CoulombAndMFOperator build both of them
# 
# <======================================================================>

import os.path
import numpy as np
from math import sqrt
import cPickle as pickle

import Plist
import Operators as op

# <======================================================================>
# GauntCoefficients are indexed by k, ml,mr where 
# k = 0,2,4,6 for f electrons
# m is the angular momentum
# Look at Gaunt Coefficient Table to understand code below (Slater Condon Text).
# Fields: 
# - k: array of denominators
# - c: array of coefficients, 
# - offset: offset for indexing, e.g. 1 for p, and 3 for f electrons
# Methods: set and get
# example: set(ml, mr, kIndex, x, sgn) => (ml and mr are indexed by offset)
# => c(ml,mr,kIndex) = sgn*sqrt(x/k[kIndex])
# => c(-ml,-mr,kIndex) = sgn*sqrt(x/k[kIndex])
# => c(mr,ml,kIndex) = -1^(|ml-mr|) sgn*sqrt(x/k[kIndex])
# => c(-mr,-ml,kIndex) = -1^(|ml-mr|) sgn*sqrt(x/k[kIndex])
# <======================================================================>
class GauntCoefficient:
    """
    Gaunt Coefficient Base Class
    """

    # ------------------------------------------------------------
    def __init__(self):
        self.debugFlag = False
        return


    # ------------------------------------------------------------
    def set(self, ml, mr, kIndex, x, sgn):

        v = sgn * sqrt(x/self.k[kIndex])

        self.c [ml+self.offset, mr+self.offset, kIndex] = v
        self.c [-ml+self.offset, -mr+self.offset,kIndex] = v
        if self.debugFlag:
            print 'Setting Gaunt ', ml, mr, kIndex, '=', v
            print 'Setting Gaunt ', -ml, -mr, kIndex, '=', v
        if (abs(ml) != abs(mr)):
            self.c [mr+self.offset, ml+self.offset, kIndex] = pow(-1, abs(ml-mr)) * v
            self.c [-mr+self.offset, -ml+self.offset, kIndex] = pow(-1, abs(ml-mr)) * v
            if self.debugFlag:
                print 'Setting Gaunt ', mr, ml, kIndex, '=', self.c [mr+self.offset, ml+self.offset, kIndex] 
                print 'Setting Gaunt ', -mr, -ml, kIndex, '=', self.c [-mr+self.offset, -ml+self.offset, kIndex] 

    # ------------------------------------------------------------
    def get(self, ml, mr, kIndex):
        return self.c[ml+self.offset,mr+self.offset,kIndex]

# <======================================================================>
# Gaunt coefficients for pp interaction
# <======================================================================>
class GauntPP(GauntCoefficient):
    """
    Gaunt coefficient for 2 electrons in p orbital
    """


    # ------------------------------------------------------------
    def __init__(self):
        self.offset = 1 # for p electrons
        GauntCoefficient.__init__(self)
        
        self.k = np.zeros( (3), float)
        self.k[0] = 1.0
        self.k[2] = 25.0


        self.c = np.zeros((3,3,3), float)
        self.set(1,1,0,1,1)
        self.set(1,1,2,1,-1)
        self.set(1,0,2,3,1)
        self.set(0,0,0,1,1)
        self.set(0,0,2,4,1)
        self.set(1,-1,2,6,-1)

# <======================================================================>
# Gaunt coefficients for ff interaction
# <======================================================================>
class GauntFF(GauntCoefficient):
    """
    Gaunt coefficient for 2 electrons in f orbital
    """

    # ------------------------------------------------------------
    def __init__(self):
        self.offset = 3 # for f electrons l = 3
        GauntCoefficient.__init__(self)
        
        self.k = np.zeros( (7), float) # indexed from l = -3 to 3

        self.k[0] = 1.0
        self.k[2] = 225.0
        self.k[4] = 1089.0
        self.k[6] = 184041

        self.c = np.zeros((7,7,7), float)

        self.set(3,3,0,1,1)
        self.set(3,3,2,25.0,-1)
        self.set(3,3,4,9.0,1)
        self.set(3,3,6,25.0,-1)

        self.set(3,2,2,25.0,1)
        self.set(3,2,4,30.0,-1)
        self.set(3,2,6,175.0,1)

        self.set(3,1,2,10.0,-1)
        self.set(3,1,4,54.0,1)
        self.set(3,1,6,700.0,-1)

        self.set(3,0,4,63.0,-1)
        self.set(3,0,6,2100.0,1)
        
        self.set(2,2,0,1.0,1)
        self.set(2,2,4,49.0,-1)
        self.set(2,2,6,900,1)
        
        self.set(2,1,2,15.0,1)
        self.set(2,1,4,32.0,1)
        self.set(2,1,6,2625.0,-1)

        self.set(2,0,2,20.0,-1)
        self.set(2,0,4,3.0,-1)
        self.set(2,0,6,5600,1)

        self.set(1,1,0,1.0,1)
        self.set(1,1,2,9.0,1)
        self.set(1,1,4,1.0,1)
        self.set(1,1,6,5625.0,-1)

        self.set(1,0,2,2.0,1)
        self.set(1,0,4,15.0,1)
        self.set(1,0,6,8750.0,1)
        
        self.set(0,0,0,1.0,1)
        self.set(0,0,2,16.0,1)
        self.set(0,0,4,36.0,1)
        self.set(0,0,6,10000.0,1)

        self.set(3,-3, 6,23100.0,-1)

        self.set(3,-2,6,11550.0,1)

        self.set(3,-1,4,42.0,1)
        self.set(3,-1,6,5250.0,-1)

        self.set(2,-2,4,70.0,1)
        self.set(2,-2,6,12600.0,1)

        self.set(2,-1,4,14.0,-1)
        self.set(2,-1,6,9450.0,-1)

        self.set(1,-1,2,24.0,-1)
        self.set(1,-1,4,40.0,-1)
        self.set(1,-1,6,10500,-1)

# <======================================================================>
# Class CoulombInteraction
# Creates the Coulomb Strength Integral I(orb)[m1,m2,m3,m4] as documented
# by Hotta Orbital ordering phenomena in d and f-electron systems p71, equation 106.
# Fields:
# - CIdebugFlag: turns on debug printing
# - onlyF0: performance optimization used if F^0 is only non-zero Slater Integral
# - L: angular momentum (e.g. f electron L=3)
# - I: 4 dimensional interaction matrix indexed by m.
# Methods:
# - get: retrieve I[m1,m2,m3.m4] where -l <= m <= l
# <======================================================================>
class CoulombInteraction():
    """
    Coulomb Interaction
    """
    # ------------------------------------------------------------
    def __init__(self, L, scp):

        if len(scp) != 2*L+1:
            raise Exception('Mismatch between L and SCP in CoulombInteraction!')

        self.CIdebugFlag = False

        self.onlyF0 = True  # helps optimize building of ham

        self.L = L
        
        oDim = L*2 + 1
        oDimList = range(oDim)
        oStart = -L
        oEnd = L + 1
        offset = L

        if L == 1:
            gaunt = GauntPP()
        elif L == 3:
            gaunt = GauntFF()
        else:
            raise Exception('Unimplemented L value in CoulombInteraction!')

        # ------------------------------------------------------------
        # Allocate the interaction strength matrix I
        # ------------------------------------------------------------
        self.I = np.zeros((oDim,oDim,oDim,oDim), float)

        # ------------------------------------------------------------
        # I = delta(m1+m2 = m3+m4) sum k=0..6 F(k) c(m1,m4,k) c(m2,m3,k)
        # ... for m=(-3...3) ...
        # ------------------------------------------------------------        
        for m1Index in oDimList:
            m1 = m1Index-offset
            for m2Index in oDimList:
                m2 = m2Index - offset
                for m3Index in oDimList:
                    m3 = m3Index - offset
                    for m4Index in oDimList:
                        m4 = m4Index - offset

                        # enforce delta function m1+m2 = m3+m4
                        # conserves angular momentum
                        if m1+m2 == m3+m4 :

                            value = 0
                            for k in range(0,len(scp),2):

                                if k > 0 and scp[k] != 0:  # helps optimize building of ham
                                    self.onlyF0 = False

                                value = value + scp[k]*gaunt.get(m1,m3,k)*gaunt.get(m4,m2,k)
                                if False:
                                    print  'k',k, 'm', m1, m2, m3, m4
                                    print scp[k], gaunt.get(m1,m3,k), gaunt.get(m4,m2,k)
                                    print 'value', value
                            self.I[m1Index,m2Index,m3Index,m4Index] = value
                            
        # ------------------------------------------------------------        
        # Debug dump of I and verify values
        # ------------------------------------------------------------        
        if self.CIdebugFlag:
            for m1 in oDimList:
                for m2 in oDimList:
                    for m3 in oDimList:
                        for m4 in oDimList:
                            value = self.I[m1,m2,m3,m4]
                            if self.CIdebugFlag and value:
                                print 'I[', m1-offset,',',m2-offset,',',m3-offset,',',m4-offset,'] = ', value
                            if abs(value -   self.I[m2,m1,m4,m3]) > 1e-10:
                                print '*** ERROR ***: I[', m2,',',m1,',',m4,',',m3,'] = ',  self.I[m2,m1,m4,m3]
                                print self.I[m2,m1,m4,m3]
                                raise Exception('Invalid Coulomb Interaction Matrix')

    # ------------------------------------------------------------        
    # get: retrieve I[m1,m2,m3.m4] where -l <= m <= l
    # ------------------------------------------------------------        
    def get(self, m1, m2, m3, m4):
        return self.I[m1+self.L, m2+self.L, m3+self.L, m4+self.L]


# <======================================================================>
# Coulomb Type Enumeration
# <======================================================================>
class CoulombTypeEnum:
    coulomb, meanfield, both = range(3)

# <======================================================================>
# Generate Coulomb Operator(s)
# - Coulomb and/or Mean-Field
# - Indexed by m and spin, uses Coulomb interaction which is indexed by m.
# Fields:
# - coulombOpsDebugFlag: turns on debug printing
# - vsList: valid states list; a list of states that participate in this Coulomb interaction
#  .... a performance optimization which saves searching states that don't matter
# - Hcoul and Hmf ... the respective operators
# Methods: 
# Pickle Extracts:
# - Uo: unitary matrix rows are the atomic states, cols are the KS orbitals,
#  ... entries are the coefficients
# - Occ: occupation of the KS orbitals
# <======================================================================>
class GenCoulombOps():
    """
    Generate the Coulomb Operator(s)
    """

    # ------------------------------------------------------------        
    def __init__(self, type, scpTuplesList):
        self.coulombOpsDebugFlag = False

        self.diagonal = True  # Optimization for using F0 only

        # Build Coulomb Operator?
        if type == CoulombTypeEnum.coulomb or type == CoulombTypeEnum.both:
            cFlag = True
        else:
            cFlag = False

        # Build Mean Field Operator?
        if type == CoulombTypeEnum.meanfield or type == CoulombTypeEnum.both:
            mfFlag = True
        else:
            mfFlag = False


        # ------------------------------------------------------------
        # Instantiate the property list and extract the atomic orbitals
        # ------------------------------------------------------------
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()        # read information from Ham.plist
        aoList = pl.aoList    # local copies of the atomic orbital list etc.
        nStates = aoList.nStates 

        # An optimization to only scan states that are involved in the Coulomb interaction
        self.vsList = np.zeros((nStates), int)
        #self.vsList = [] # valid state list 

        # If not specified by the user, read the property list values
        # This allows us to iterate using different scp values
        if scpTuplesList == []:
            scpTuplesList = pl.slaterCondonParameters()
        print 'Slater Condon Parameters: ', scpTuplesList

        # ------------------------------------------------------------
        # Convert the generic interaction matrix to one that uses our
        # basis set and supports all of the atomic states used  for the
        # 2 particle operator
        # ------------------------------------------------------------
        if os.path.isfile('Uo.pkl') and os.path.isfile('Occ.pkl'):
            genExpValue = True
        else:
            if mfFlag:
                raise Exception('MeanField operator without Uo.pkl or Occ.pkl')
            genExpValue = False

        if genExpValue:
            f = open('Uo.pkl', 'r')
            Uo = pickle.load(f)
            f.close()
            f = open('Occ.pkl', 'r')
            Occ = pickle.load(f)
            f.close()
            #
            # Uocc is a nStates x nKSStates matrix ...
            # ... where the unoccupied KS states (columns) are zero
            #
            Uocc = np.mat(np.multiply(Uo, Occ))
            #
            # A matrix of ground state expectation values i.e. <f^\dagger_i f_j>_0
            #
            gsExpMx = Uocc*Uocc.transpose().conjugate()
            if False:
                gsExpMx = np.zeros((nStates, nStates))
                for row in range(gsExpMx.shape[0]):
                    for col in range(gsExpMx.shape[1]):
                        if row and col < 7:
                            gsExpMx[row, col] = .1
                for row in range(gsExpMx.shape[0]):
                    for col in range(gsExpMx.shape[1]):
                        if abs(gsExpMx[row, col]) > 1e-10: 
                            print 'gsExpMx[', row, col, '] =', gsExpMx[row, col]
            
        if cFlag:
            self.Hcoul=np.zeros((nStates, nStates, nStates, nStates), complex)

        if mfFlag:
            self.Hmf=np.zeros((nStates, nStates), complex)
            J = np.zeros((nStates, nStates), complex)
            K = np.zeros((nStates, nStates), complex)

        # ------------------------------------------------------------
        # Loop through all the Slater Condon List, building the operator 
        # all Coulombic interactions. scp = Slater Condon Parameter
        # scpTuplesList contains 1 or more scpTuples
        # ... each scpTuple contains the ang mo L and a list of associated scp's 
        # ------------------------------------------------------------
        for scpTuple in scpTuplesList:
            L = scpTuple[0]
            scp = scpTuple[1]
                
            I = CoulombInteraction(L, scp)
            if not I.onlyF0:
                print 'Non-diagonal Operator'
                self.diagonal = False
            else:
                print 'Diagonal Operator'

            # ------------------------------------------------------------
            # Create Coulomb Interaction Matrix I for the basis
            # ... basically a copy of the I into Ifop converting the values as 
            # ... necessary.  Lots of looping for just a copy.
            # only loop through ao's that are f-orbitals
            # ------------------------------------------------------------
            for ao1 in aoList.aoByLGen(L):
                self.vsList[ao1.aoIndex] = 1
                for sao1 in ao1.sphericalExpGen():
                    lz1 = sao1.ao.lz
                    c1 = sao1.c.conjugate()# !! complex conjugate of the coeff
                    
                    for ao2 in aoList.aoByLGen(L):
                        for sao2 in ao2.sphericalExpGen():
                            lz2 = sao2.ao.lz
                            c2 = sao2.c.conjugate()# !! complex conjugate of the coeff

                            if ao1 is not ao2:
                                for ao3 in aoList.aoByLGen(L):
                                    for sao3 in ao3.sphericalExpGen():
                                        lz3 = sao3.ao.lz
                                        c3 = sao3.c 
                                        for ao4 in aoList.aoByLGen(L):
                                            for sao4 in ao4.sphericalExpGen():
                                                if ao1.sz == ao3.sz and ao2.sz == ao4.sz and \
                                                        ao3 is not ao4:
                                                    lz4 = sao4.ao.lz
                                                    c4 = sao4.c                                                         
                                                        #if value and self.coulombOpsDebugFlag:
                                                            #print 'I[', ao1.aoIndex,ao2.aoIndex,ao3.aoIndex,ao4.aoIndex, ']', \
                                                                #value, lz1, lz2, lz3, lz4

                                                    value = I.get(lz1, lz2, lz3, lz4)
                                                    coeff = c1*c2*c3*c4

                                                # convert the generic value to our basis set
                                                    if cFlag and value:
                                                        # Coulomb operators conserves spin - spin delta function
                                                        # Don't need the 1/2 factor because the operator does not 
                                                        # double count
                                                        self.Hcoul[ao1.aoIndex,ao2.aoIndex,ao3.aoIndex,ao4.aoIndex] = \
                                                            (self.Hcoul[ao1.aoIndex,ao2.aoIndex,ao3.aoIndex,ao4.aoIndex] \
                                                                 + value * coeff)

                                                    if mfFlag:
                                                    #print ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, coeff, value
                                                    
                                                        # exchange K and direct J integrals for mean field
                                                        # K forces the spins to be the same
                                                        # The 1/2 factor is added in later
                                                    
                                                        J[ao1.aoIndex, ao3.aoIndex] = \
                                                            J[ao1.aoIndex, ao3.aoIndex] \
                                                            + value \
                                                            * gsExpMx[ao2.aoIndex, ao4.aoIndex] \
                                                            * coeff

                                                        K[ao1.aoIndex, ao4.aoIndex] = \
                                                            K[ao1.aoIndex, ao4.aoIndex] \
                                                            + value \
                                                            * gsExpMx[ao2.aoIndex, ao3.aoIndex] \
                                                            * coeff
        # ------------------------------------------------------------
        # robustness check
        # ------------------------------------------------------------
        if cFlag:
            if genExpValue:
                expValue = 0
            statesList = range(nStates)
            for mu1Index, mu1 in enumerate(statesList):
                for mu2 in statesList:
                    for mu3Index, mu3 in enumerate(statesList):
                        for mu4 in statesList:
                            
                            if self.Hcoul[mu1,mu2,mu3,mu4]:
                                if genExpValue:
                                    expValue = expValue + 0.5*self.Hcoul[mu1,mu2,mu3,mu4]*(gsExpMx[mu1,mu3]*gsExpMx[mu2,mu4] - gsExpMx[mu1,mu4]*gsExpMx[mu2,mu3])
                                delta =  self.Hcoul[mu1,mu2,mu3,mu4] - self.Hcoul[mu2,mu1,mu4,mu3]
                                if abs(delta) > 1e-15:
                                    print ' *** ERROR ****'
                                    print 'Hcoul[', mu1,mu2,mu3,mu4,'] = ', self.Hcoul[mu1,mu2,mu3,mu4]
                                    print 'Hcoul[', mu2,mu1,mu4,mu3,'] = ', self.Hcoul[mu2,mu1,mu4,mu3]
                                    print 'Delta I[', mu1,mu2,mu3,mu4,'] = ',delta
                                    raise Exception()
            if genExpValue:
                print '<Hcoul> = ', expValue

        # ------------------------------------------------------------
        # Create One Particle Operator for Mean-field
        # ------------------------------------------------------------
        if mfFlag:
            expValue = 0
            for row in xrange(nStates):
                for col in xrange(nStates):
                    if abs(J[row, col]) < 1e-10:
                        J[row,col] = 0
                    if abs(K[row, col]) < 1e-10:
                        K[row,col] = 0
                    if abs(J[row,col].imag) < 1e-10:
                        J[row,col] = J[row,col].real
                    if abs(K[row,col].imag) < 1e-10:
                        K[row,col] = K[row,col].real

                    self.Hmf[row,col] = self.Hmf[row,col] + 0.5*(J[row,col] - K[row,col])

                    if False:
                        if J[row, col] != 0:
                            print 'J[', row, col, '] = ', J[row, col]
                        if K[row, col] != 0:
                            print 'K[', row, col, '] = ', K[row, col]

                    expValue = expValue + self.Hmf[row,col]*gsExpMx[row,col]

            # compute gs expectation value
            print '<Hmf> = ', expValue
            
# <======================================================================>
# Coulomb Operator - Two Particle Operator
# - Optional Slater Condon Parameters provided in tuples of 
# ... (L, (F0, F2, ...))
# Fields: 
# Methods: 
# - dbgProcessMatrixElements: calls operators' processors
# <======================================================================>
class CoulombOperator(op.TwoParticleOperator):
    # ------------------------------------------------------------        
    def __init__(self, scpTuplesList = []):
        co = GenCoulombOps(CoulombTypeEnum.coulomb, scpTuplesList)
        op.TwoParticleOperator.__init__(self, co.Hcoul, 'Hcoul')
        self.setValidStates(co.vsList)
        self.setOptimized()
        if co.diagonal:
            self.setDiagonalOp()

    # ------------------------------------------------------------        
    def dbgProcessMatrixElements(self, mbHam):
        mbHam = self.processMatrixElements(mbHam)
        #mbHam.dump()
        return mbHam

# <======================================================================>
# Mean Field Operator - One Particle Operator
# - Optional Slater Condon Parameters provided in tuples of 
# ... (L, (F0, F2, ...))
# Fields: 
# Methods: 
# - dbgProcessMatrixElements: calls operators' processors
# <======================================================================>
class MeanFieldOperator(op.OneParticleOperator):

    # ------------------------------------------------------------        
    def __init__(self, scpTuplesList = []):
        co = GenCoulombOps(CoulombTypeEnum.meanfield , scpTuplesList)
        op.OneParticleOperator.__init__(self, co.Hmf, 'Hmf')

    # ------------------------------------------------------------        
    def dbgProcessMatrixElements(self, mbHam):
        mbHam = self.processMatrixElements(mbHam)
        return mbHam

# <======================================================================>
# Coulomb and Mean Field Operators
# - Optional Slater Condon Parameters provided in tuples of 
# ... (L, (F0, F2, ...))
# Fields: 
# - Hcoul - Coulomb Operator
# - Hmf - Mean-Field Operator
# Methods: 
# - dump: calls operator dump
# - processMatrixElements: calls operators' processors
# <======================================================================>
class CoulombAndMFOperators():
    """
    Coulomb and Mean Field Operators
    """
    # ------------------------------------------------------------        
    def __init__(self, scpTuplesList = []):
        co = GenCoulombOps(CoulombTypeEnum.both, scpTuplesList)
        self.Hcoul = op.TwoParticleOperator(co.Hcoul, 'Hcoul')
        self.Hcoul.setValidStates(co.vsList)
        if co.diagonal:
            self.Hcoul.setDiagonalOp()
        self.Hcoul.setOptimized()
        self.Hmf = -op.OneParticleOperator(co.Hmf, 'Hmf')

    # ------------------------------------------------------------        
    def dump(self):
        self.Hcoul.dump()
        self.Hmf.dump()

    # ------------------------------------------------------------        
    def processMatrixElements(self, mbHam):
        self.Hcoul.processMatrixElements(mbHam)
        self.Hmf.processMatrixElements(mbHam)
        return mbHam

    


    
