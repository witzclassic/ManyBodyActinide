#!/usr/bin/env python

# HamClasses.py 
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
# Contains classes to support building a many body Hamiltonian:
# 
# - ParticleOperator base class of:
# -- OneParticleOperator
# -- TwoParticleOperator
# <======================================================================>

import struct
import numpy as np
import math
import Bit
import cPickle as pickle
import Plist
import ManyBodyState

# <======================================================================>
# InteractionTensor: 4 Dimensional 
# Fields: 
# - I: the tensor
# Methods:
# - add: adds value to entry
# <======================================================================>
class InteractionTensor:
    # ------------------------------------------------------------
    def __init__(self, nStates, type=complex):
        self.I = np.zeros((nStates,nStates,nStates,nStates), type)

    # ------------------------------------------------------------
    def add(self, i1,i2,i3,i4,value):
        self.I[i1,i2,i3,i4] = self.I[i1,i2,i3,i4] + value

# <======================================================================>
# --------------------ParticleOperator-----------------------------------
# Base class for all operators
# Fields:
# - debug flags
# - op: 2D or 4D tensor
# - nStates: number of sps (single particle states)
# - validStatesList: list of states to consider for this operator
# - diagonalOp: optimization for performance
# - name: text name
# - optimized:  optimization for performance
# Methods:
# - setOptimized
# - setValidStates
# - processMatrixElements: loop thru mb states, generate <l | Op | r>
# - (un)pickle: (restore) save op
# <======================================================================>
class ParticleOperator:
    """
    Virtual base class for the real particle operators
    - processMatrixElements drives the non-zero matrix element subroutines
    """

    # -----------------------------------------------------------------------
    def __init__(self, op, name):

        # True dumps out a lot of info to terminal
        #self.particleOperatorDebugFlag = True
        self.particleOperatorDebugFlag = False
        self.particleOperatorDeepDebugFlag = False

        self.op = op           # 2D or 4D nState Tensor
        self.nStates = self.op.shape[0] # number of one body states

        # Assume all states are valid for this operator
        self.validStatesList = np.ones(self.nStates) 
        
        # Assume that the operator is not diagonal in this basis set
        self.diagonalOp = False

        # Useful for debugging
        self.name = name
        
        # Another performance hack
        self.optimized = False

    # -----------------------------------------------------------------------
    def setOptimized(self):
        print self.name, 'Optimized!'
        self.optimized = True

    # ------------------------------------------------------------
    # Helps to process only the states specific to this operator
    # Performance optimization only
    # ------------------------------------------------------------
    def setValidStates(self, vsList):
        self.validStatesList = vsList

    # -----------------------------------------------------------------------
    # processMatrixElements - loops through mb states calling the sub-class
    #  .. generator of non-zero matrix elements <lState| Op | rState>
    # -----------------------------------------------------------------------
    def processMatrixElements(self, mbHam):

        # Read in the state table pickles
        f = open('stateTable.pkl', 'r')
        stateTable = pickle.load(f)
        f.close()
        f = open('stateTableDict.pkl', 'r')
        stateTableDict = pickle.load(f)
        f.close()


        statusThreshold = 10000
        nextStatusIndex = 0

        # Initialize the structures of the many body state
        rmbs = ManyBodyState.ManyBodyState(self.nStates, self.validStatesList, self.optimized)

        print
        print '<-------- Processing Matrix Elements (', self.name, ')------->'
        print
        
        stateRange = 0

        # loop through each mb state processing the non-zero matrix elements
        for rStateIndex, rState in enumerate(stateTable):

            if nextStatusIndex == rStateIndex: # print some status
                state = Bit.BitStr(self.nStates, rState)
                print  'StateIndex %d value %s' % (rStateIndex, state.display())
                nextStatusIndex = nextStatusIndex + statusThreshold

            if self.particleOperatorDebugFlag: # use to debug
                state = Bit.BitStr(self.nStates, rState)
                print  'StateIndex %d value %s' % (rStateIndex, state.display())
            
            if stateRange:
                stateRange = stateRange - 1
            else:
                # generate any non-zero matrix elements associated with this state and ham
                stateRange = rmbs.setMBS(rState)
                #rmbs.dump()
                #print 'Debug Operators stateRange = ', stateRange
                self.genNonZeroMatrixElements(rmbs, rStateIndex, mbHam, stateTableDict, stateRange)
                stateRange = stateRange - 1

    # ------------------------------------------------------------
    def pickle(self):
        f = open(self.name+'.pkl', 'w')
        pickle.dump(self, f)
        f.close()

    # ------------------------------------------------------------
    def unpickle(self):
        f = open(self.name+'.pkl', 'r')
        self = pickle.load(f)
        f.close()

# <======================================================================>
# ---------------------OneParticleOperator---------------------------
# Inherits ParticleOperator
# Fields:
# Methods:
# Builtins __add__, __iadd__, __sub__, __isub__, __neg__
# - Eigenvalues: diagonalizes single particle operator (sps)
# - dump
# - modify: set entry to value
# - isHermitian: robustness check to check Hermiticity
# - genNonZeroMatrixElements: called by processMatrixElements
# - ksExpValue: 
# ... uses pickle files Uo and Occ, generates exp value for operator using 
# ... KS ground state    
# <======================================================================>
class OneParticleOperator(ParticleOperator):
    """
    A one particle operator (Fermionic): sum_{i,j} O_{ij} c^\dagger_i c_j
    """

    # ------------------------------------------------------------
    def __init__(self, op, name='OneParticleOperator'):
        ParticleOperator.__init__(self, op, name)
        if  not self.isHermitian():
            raise Exception('Non Hermitian Operator')

    # ------------------------------------------------------------
    def __add__(self, other): # op + op 
        if self.nStates != other.nStates:
            raise Exception('mismatched operators')
        print 'Adding ', self.name, other.name
        self.name = self.name + other.name
        self.op = self.op + other.op
        return self

    # ------------------------------------------------------------
    def __iadd__(self, other):  # += 
        if self.nStates != other.nStates:
            raise Exception('mismatched operators')
        self.name = self.name + other.name
        self.op = self.op + other.op
        return self

    # ------------------------------------------------------------
    def __sub__(self, other): # op - op
        if self.nStates != other.nStates:
            raise Exception('mismatched operators')
        self.name = self.name + other.name
        self.op = self.op - other.op
        return self

    # ------------------------------------------------------------
    def __isub__(self, other): # -= operation
        if self.nStates != other.nStates:
            raise Exception('mismatched operators')
        self.name = self.name + other.name
        self.op = self.op - other.op
        return self

    # ------------------------------------------------------------
    def __neg__(self):  # change the sign of the operator
        self.op = -self.op
        return self

    # ------------------------------------------------------------
    # generate and print the single state eigenvalues
    # ------------------------------------------------------------
    def eigenvalues(self, nFermions=0):
        
        print
        print self.name, '<=== One Particle Eigenvalues: ===> '

        [self.e,self.v]=np.linalg.eigh(self.op)
        if nFermions == 0:
            print self.e
        else:
            n = math.floor(nFermions)
            print 
            print self.name, 'Ground State Energy(N = ', nFermions, ') = ', sum(self.e[0:n]) + (nFermions-n)*self.e[n+1]

        for e in self.e:
            print e
        
    # ------------------------------------------------------------
    def dump(self):
        print '<<<<< -------- ', self.name, '------------- >>>>'
        for row in range(self.nStates):
            for col in range(self.nStates):
                if abs(self.op[row,col]) > 1e-10:
                    print '[', row, ',', col,']', '=',self.op[row,col]
        
    # ------------------------------------------------------------
    # May be used to modify operator matrix entries
    # ------------------------------------------------------------
    def modify(self, r, c, v):
        self.op[r,c] = v

        
    # ---------------- hermitian ? -----------------------------------
    # Good robustness check
    # ------------------------------------------------------------

    def isHermitian(self):
        stateList = range(self.nStates)
        for index, r in enumerate(stateList):
            for c in stateList[index:]:
                if abs(self.op[r,c] - self.op[c,r].conjugate()) > 1e-10:
                    print r,c, self.op[r,c], self.op[c,r]
                    return False
        return True
                
                
    # ----------------genNonZeroMatrixElements--------------------
    # Two cases: 
    #     1. Diagonal elements which sum over all of the occupied states
    #     2. Off-diagonal elements which can differ only by one occ
    # state
    #     Note, only generate half the entries since MBHam is
    # Hermitian.  Therefore, lState > rState, always
    # 
    # PERFORMANCE is king!  Be careful about this routine, it needs to scale
    # to very large Many Body Hamiltonians.
    # ------------------------------------------------------------
    def genNonZeroMatrixElements(self, rmbs, rStateIndex, mbHam, stateTableDict, stateRange):

        # Be careful ... swap changes rmbs value, use rmbs.getValue() !!not!! rmbs.value

        lowestZero = 0
        diagValue = 0

        # --------------------------------------------------
        # Loop through the occupied states 
        # --------------------------------------------------
        for occState,nextEmptyStateIndex in zip(rmbs.occList, rmbs.nextEmptyList):

            # accrue diagonal matrix element while we go along
            diagValue = diagValue + self.op[occState,occState]

            # Move the occupied state to next highest unocc state
            for emptyState in rmbs.emptyList[nextEmptyStateIndex:]:

                # process if h is non-zero
                m1Value =   self.op[emptyState,occState]
                if m1Value:
                    if (rmbs.parity[occState] - rmbs.parity[emptyState] == 0):
                        m1Value = -m1Value

                    # lStateIndex is retrieved by generating the mbs by swapping the 
                    # empty and occupied states and then looking up the quantum ID (state index)
                    # in the state table dictionary
                    lStateIndex = stateTableDict[rmbs.swapTwoStates(emptyState, occState)] 
                        
                    mbHam.setMatrixElement(lStateIndex,rStateIndex,m1Value)
                    if self.particleOperatorDebugFlag:
                        print 'm1Value:', m1Value, stateTableDict[lmbs.value], rStateIndex

        # Finally, add diagonal elements to the mbHam
        if diagValue:
            if self.particleOperatorDebugFlag:
                print 'm0Value:', diagValue, rStateIndex
            mbHam.setDiagMatrixElement(rStateIndex, diagValue)

    # ----------------------------------------------------------------------
    # Debug routine to generate expectation values, not sure it still works
    # ----------------------------------------------------------------------
    def ksExpValue(self):
        self.dump()
        f = open('Uo.pkl', 'r')
        Uo =pickle.load(f)
        f.close()
        f = open('Occ.pkl', 'r')
        occList =pickle.load(f)
        f.close()
        
        expValue = 0.0
        for ks in occList:
            if ks:
                for st1 in range(self.nStates):
                    if Uo[st1, ks] != 0.0:
                        for st2 in range(self.nStates):
                            print self.op[st1, st2], Uo[st1,ks], Uo[st2,ks]
                            expValue = expValue + self.op[st1, st2]*Uo[st1,ks]*Uo[st2,ks].conjugate()

        print 'expValue = ', expValue

# <======================================================================>
# ---------------------TwoParticleOperator ------------------------------
# 
# Creates an array of non-zero matrix elements (nzMtx)
# associated with this single-body op and this state:
#
#                   <lState | Op | rState>
# ------------------------------------------------------------
# Inherits ParticleOperator
# Fields:
# Methods:
# - setDiagonalOp: set flag to optimize processing
# - dump
# - modify - set entry to value
# - genNonZeroMatrixElements - called by base class processMatrixEntries
# - ksExpValue: find KS expectation value, uses pickle Uo and Occ files
# <======================================================================>
class TwoParticleOperator(ParticleOperator):
    """
    TwoParticleOperator: (Fermionic)
    sum_{i,j,k,l} Op_{i j k l} c^\dagger_i c^\dagger_j c_l c_k
    """

    # ------------------------------------------------------------
    def __init__(self, op, name='TwoParticleOperator'):
        ParticleOperator.__init__(self, op, name)

    # ------------------------------------------------------------
    # Useful for things like F0 only in Coulomb operators ...
    # much faster processing
    # ------------------------------------------------------------
    def setDiagonalOp(self):
        self.diagonalOp = True

    # ------------------------------------------------------------
    def dump(self):
        print '<<<<< -----------', self.name, ' ------------ >>>>'
        print 'nStates = ', self.nStates
        for m1 in range(self.nStates):
            for m2 in range(self.nStates):
                for m3 in range(self.nStates):
                    for m4 in range(self.nStates):
                        if self.op[m1,m2,m3,m4]:
                            print 'op[', m1, m2, m3, m4,'] = ', self.op[m1,m2,m3,m4]
        print 

    # ------------------------------------------------------------
    # Used for creating the Op
    # ------------------------------------------------------------
    def modify(self, m1,m2,m3,m4, v):
#        print 'Setting ',r,c,'to',v
        self.op[m1,m2,m3,m4] = v
        
    # ------------------------------------------------------------
    # Meat of the operators, generate the non-zero matrix elements
    # 
    # We are looking thru the occupied states of |rState>
    # always generating the !! higher !! numbered state <lState|.
    # There are three cases (1) diagonal elements, (2) one state
    # changes and one stays the same and (3) both states change.
    # 
    # Processing only higher states cuts the routine in half as the Ham is Hermitian 
    # and we could (if needed) copy one half to the other half.  Generally, 
    # this isn't needed as the diag routines understand the sparse mtx format.
    #
    # Two Body Op generator of non-zero matrix elements:
    #
    # ------------------------------------------------------------
    def genNonZeroMatrixElements(self, rmbs, rStateIndex, mbHam, stateTableDict, stateRange):

        # Be careful ... swap changes value, use rmbs.getValue() !!not!! rmbs.value

        mx0Value = 0  # matrix value zero states change ... i.e. the diagonal entry

        # --------------------------------------------------
        # Loop through all of the occupied states 
        # --------------------------------------------------
        for sIndex, s in enumerate(rmbs.occList):
            if self.particleOperatorDeepDebugFlag:
                print 's=', s
            for r in rmbs.occList[:sIndex]:
                if self.particleOperatorDeepDebugFlag:
                    print 'outer r=', r
                # ==============> Diagonal Entry <===============
                # sum_{s>r} (Self.Op[r,s,r,s] - Self.Op[r,s,s,r])
                #
                # Nothing else for it, it is a N*N-1 loop for occ states
                # --------------------------------------------------
                mx0Value = mx0Value + (self.op[r,s,r,s]-self.op[r,s,s,r])
                if mx0Value and self.particleOperatorDebugFlag:
                    print 'mx0Value = ', mx0Value, r,s
                # ==============> Diagonal Entry <===============
            
            # ---------------------------------------------------------
            # If we know that it is a diagonal operator, just move on
            # much, much faster
            # ---------------------------------------------------------
            if self.diagonalOp:
                continue
            
            # ---------------------------------------------------------
            # Loop through all of the unoccupied states greater than s
            # !!! t > s !!!
            # ---------------------------------------------------------
            nextEmptyStateIndex = rmbs.nextEmptyList[sIndex]
            for tIndex, t in enumerate(rmbs.emptyList[nextEmptyStateIndex:]):
                if self.particleOperatorDeepDebugFlag:
                    print 't=', t
                # --------------------------------------------------
                mx1Value = 0 # matrix value one state changes s->t, r remains the same
                for r in rmbs.occList:
                    if s == r:  # double destructor is always zero
                        continue
                    if self.particleOperatorDeepDebugFlag:
                        print 'inner r=', r

                    # ==============> Swap one Particle <===============
                    # swap s & t
                    # -(-1^q sum_{s\ne r,t} (O_{r,t,r,s} - O_{r,t,s,r}))
                    # Overall minus because we need -parity[t]
                    # --------------------------------------------------
                    mx1Value = (self.op[r,t,s,r]-self.op[r,t,r,s])
                    # --------------------------------------------------
                    # If we processed mx1Value, add it to the mbHam
                    # --------------------------------------------------
                    if mx1Value:
                        if rmbs.parity[t]+rmbs.parity[s] == 0:
                            mx1Value = -mx1Value

                        if self.particleOperatorDebugFlag:
                            print 'mx1Value ', rStateIndex,  '=', mx1Value, s,t,r

                        lStateIndex = stateTableDict[rmbs.swapTwoStates(t, s)]
                        for offset in xrange(stateRange):
                            mbHam.setMatrixElement(lStateIndex+offset,
                                                   rStateIndex+offset, 
                                                   mx1Value)
                    # ==============> Swap one Particle <===============

                    if s > r:
                        if self.particleOperatorDeepDebugFlag:
                            print 's > r'
                        # --------------------------------------------------
                        # Loop thru unoccupied states looking to swap two particles
                        # t(unocc) > s(occ) > r(occ)
                        #
                        # ====> since we loop u only up to t, u < t <=====
                        #
                        # --------------------------------------------------
                        for u in rmbs.emptyList[:nextEmptyStateIndex+tIndex]:

                            if self.particleOperatorDeepDebugFlag:
                                print u, t, r, s
                            # since t > s > r, we always get a higher state
                            # t > u
                            #
                            # substitute in QM Book v->r, u->s, s->t, r->u
                            # <> = \pm (-1)^(t+u) (-1)^(r+s) ([u,t,s,r] - [u,t,r,s])
                            # plus  if (t < u && s < r) || (u < t && r < s)
                            # minus if (t < u && s > r) || (u < t && r > s)
                            #
                            # In our case, u<t && r<s ==> plus
                            #
                            # Parity is computed on rState, t and u need to be adusted for lState
                            # ==> r && s parity remains the same
                            # 3 cases
                            # if u < r, u stays
                            # if r < u < s, u changes
                            # if s < u, u stays
                            # t parity changes if s < u
                            #
                            mx2Value = (self.op[t,u,s,r] - self.op[t,u,r,s])
                            if mx2Value:
                                lStateIndex = stateTableDict[rmbs.swapFourStates(t, u, r, s)]
                                
                                # compute the parity
                                pty =  rmbs.parity[r] + rmbs.parity[s] - rmbs.parity[t] + rmbs.parity[u]
                                if pty == 2 or pty == -2:
                                    mx2Value = -mx2Value

                                # if r < u < s:
                                if r < u and u < s:
                                    mx2Value = -mx2Value

                                
                                if self.particleOperatorDebugFlag:
                                    print 'mx2Value ', lStateIndex, rStateIndex, '=', mx2Value, r,s,t,u
                                for offset in xrange(stateRange):
                                    mbHam.setMatrixElement(lStateIndex+offset,
                                                           rStateIndex+offset, 
                                                           mx2Value)
                                # there can only be one set of 4 different indices
                                # no need to search for more
                                # break 
                    # --------------------------------------------------
                    #  End for r in rmbs.ones 
                    # --------------------------------------------------

        # --------------------------------------------------
        # Flush out diag 
        # --------------------------------------------------
        if (mx0Value):
            for offset in xrange(stateRange):
                mbHam.setDiagMatrixElement(rStateIndex+offset, mx0Value)
                                   
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    def ksExpValue(self):
        f = open('Uo.pkl', 'r')
        Uo =pickle.load(f)
        f.close()
        f = open('Occ.pkl', 'r')
        occList =pickle.load(f)
        f.close()
        
        expValue = 0.0
        for ks in occList:
            if ks:
                for st1 in range(self.nStates):
                    if Uo[st1, ks] != 0.0:
                        for st2 in range(self.nStates):
                            if Uo[st2, ks] != 0.0:
                                for st3 in range(self.nStates):
                                    if Uo[st3, ks] != 0.0:
                                        for st4 in range(self.nStates):
                                            if Uo[st4, ks] != 0.0:
                                                expValue = expValue + self.op[st1, st2, st3, st4]*Uo[st1,ks]*Uo[st2,ks]*Uo[st3,ks].conjugate()*Uo[st4,ks].conjugate()


        print 'expValue = ', expValue

