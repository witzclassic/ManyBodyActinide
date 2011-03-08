#!/usr/bin/env python

# stateTable.py - Implementation of a Python interface to build a
#         many-body state table.  Each state consists of a bit
#         string where the occupied states are represented by a '1'. 
#
# Copyright (C) 2009 Brown University Physics, JB Marston
# Author: Steve Horowitz
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
# Loop through all of the possible many-body states adding a key,value
#         pair to the python dictionary "stateTable". Each key is the
#         bit string representing the state occupancy and the value is
#         the cardinal ordering of the state (1,2,3, ...)
# Classes:
# - ManyBodyState
# <======================================================================>

import CG
import Bit
import sys
import cPickle as pickle
import Plist
import numpy as np

# <======================================================================>
# Class ManyBodyState
# 
# Used to inherity from Bit string, but the peformance hit is too high
# Best to call init once to preallocate the structures and then call setValue
#
# Fields: 
# - nStates: total number of single particle states (sps)
# - parity: used for swapping sign
# - vsList: list of valid single particle states
# ... optimization for operators to skip some states
# - optStates: another optimization to move the single particle
# ... states that do not participate in 2-body interactions to the low
# ... order .. see code below for more description.
# - occList, nextEmptyList, emptyList ... lists that aid in operator processing
# ... occ is which valid states are occupied
# ... emptyList is list of sps that are empty
# ... nextEmptyList is the list empty positions after this occ. pos.
# - mbs - the many body state ... essentially a bit string representing 
# ... the many body state in the occupation basis.
# Pickle:
# .. nParticles: number of particles in many body state
# Methods:
# - dump
# - setMBS ... init structures associated with this state
# - getMBS - return mbs
# - setWithParity (pos): mark the sps occupied in the mbs, return fermion parity
# - swapTwoStates(c, d) : mark sps c occupied and clear sps d.
# - swapFourStates(c1, c2, d1, d2) c becomes occupied and d becomes clear
# ... c for constructor operator (e.g. f^\dagger) and d for destructor (e.g. f)
# - support routines: set, clear ... slow be wary of calling
# - get: returns the occupation of the sps
# - display: print out the mbs
# - __getitem__ : supports builtin a = mbs[index]
# <======================================================================>
class ManyBodyState():
    """
    ManyBodyState 
    - input vsList - valid state list
    - occList - a list of occupied states merged with vsList
    - emptyList - a list of empty states merged with vsList
    - nextEmptyList - for this occupied state, an index into emptyList that is past this state
    """

    # ------------------------------------------------------------
    def __init__(self, nStates, vsList, optimized=False):
        
        f = open('nParticles.pkl', 'r')
        self.nParticles = pickle.load(f)
        f.close()

        self.nStates = nStates
        self.parity = np.zeros((nStates), int)
        self.vsList = vsList


        # If the stateTable is optimized for this operator (a real
        # hack, but it works), then we can calculate one non-zero
        # element nze = <LHS | Op | RHS> and just write nze to all of
        # the other states which just vary in the non-operated states.
        
        # More specifically, the f-coulomb operator only operates on f
        # electrons and their interactions.  We can organize the p
        # electrons to be the lower order single particle states (sps)
        # such state MBS = |f electrons, p electrons>.  Then, the
        # following values are equal:
        # <LHS | Op | RHS>, <LHS+1|Op|RHS+1>, ..., <LHS+n|Op|RHS+n>
        # where n is the number of states that can be formed by
        # varying the p-electrons among the p states.

        # So, we need to find n = p electrons choose p electron states
        # p electrons is nParticles - f electrons
        # p states is fixed (well, based on the vsList)

        if optimized: 
            print 'ManyBodyState is optimized.'
            self.optStates = nStates - sum(vsList)
        else:
            self.optStates = 0

    # ------------------------------------------------------------
    def dump(self):
        print self.nStates, self.optStates
        print self.occList
        print self.nextEmptyList
        print self.emptyList
        print self.parity
        print self.vsList
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    def setMBS(self, mbs):
            
        self.mbs = mbs

        parity = 1
        occPos = 0
        empPos = 0
        self.occList = []
        self.nextEmptyList = []
        self.emptyList = []

        # Loop through each of the sps's:
        for s in xrange(self.nStates):
            self.parity[s] = parity # set the parity

            if self.mbs & (0x1 << s): # if s is an occupied state

                if self.vsList[s]:  # if it is a valid state for the operator
                    self.occList.append(s)
                    self.nextEmptyList.append(empPos)
                    occPos = occPos + 1

                parity = -parity # if occupied, change the next state parity

            elif self.vsList[s]:
                self.emptyList.append(s) # mark unoccupied
                empPos = empPos + 1

        # So, we need to find n = p electrons choose p electron states
        # p electrons is nParticles - f electrons
        # p states is fixed (well, based on the vsList)
        # see discussion in __init__

        if self.optStates > 0:
            #print self.optStates, occPos
            return int(CG.binomialCoeff(self.optStates, self.nParticles - occPos))
        else:
            return 1
            

    # ------------------------------------------------------------
    def getMBS(self):
        return self.mbs
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    def setWithParity(self, pos):
        parity = 1
        for st in range(self.nStates-1, pos, -1):
            #if self.get(st):
            if self.mbs & (0x1 << st):
                parity = -parity
        #self.set(pos)
        self.mbs = self.mbs | (0x1 << pos)
        return parity
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    def swapTwoStates(self, c, d):
        st = self.mbs
        #self.set(c)
        st = st | (0x1 << c)
        #self.clear(d)
        st = st & ~(0x1 << d)
        return st
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    def swapFourStates(self, c1, c2, d1, d2):
        st = self.mbs
        #self.set(c1)
        #self.set(c2)
        st = st | ((0x1 << c1) + (0x1 << c2))
        #self.clear(d1)
        #self.clear(d2)
        st = st & ~((0x1 << d1) + (0x1 << d2))
        return st
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # DONT CALL THE FOLLOWING if performance is key .. do it inline
    # ------------------------------------------------------------
    def clear(self, pos):
        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr clear) Pos'+str(pos)+'out of range '+str(self.len)

        # --------------------------------------------------------------------

        self.mbs = self.mbs & ~(0x1 << pos)

    # ------------------------------------------------------------
    def set(self, pos):

        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr set) Pos'+str(pos)+'out of range '+str(self.len)
        # --------------------------------------------------------------------

        self.mbs = self.mbs | (0x1 << pos)

    # ------------------------------------------------------------
    def get(self, pos):
        # --------------------------------------------------------------------
        #if pos  >= self.len:
        #raise Exception('(BitStr get) Pos'+str(pos)+'out of range '+str(self.len)
        if self.mbs & (0x1 << pos):
            return True
        else:
            return False

    # ------------------------------------------------------------
    def display(self):
        # Slow, so rewrite if we need performance
        bstr = '0b'
        bOccList = []
        for b in range(self.nStates,0,-1):
            if self.get(b-1):
                bstr += '1'
                bOccList.append(b-1)
            else:
                bstr += '0'
        return bstr+str(bOccList)

    # ------------------------------------------------------------
    def __getitem__(self, index):
        return self.get(index)


# <======================================================================>
# ----------------------------------------------------------------------
# Analyze particle occupations in the ground state eigenvector
# ----------------------------------------------------------------------
# <======================================================================>
def analyzeOccupations(vec):

    # Read in the state table pickles
    f = open('stateTable.pkl', 'r')
    stateTable = pickle.load(f)
    f.close()

    f = open('nParticles.pkl', 'r')
    nParticles = pickle.load(f)
    f.close()

    pl = Plist.AtomicOrbitalPlist()
    pl.readPlist()        # read information from Ham.plist
    aoList = pl.aoList    # local copies of the atomic orbital list etc.
    nStates = aoList.nStates 

    occ = np.zeros(nStates, dtype=float)
    polHist = np.zeros((2*7)+1, dtype=float)
    fOcc = np.zeros((2*7)+1, dtype=float)

    vsList = np.ones(nStates)
    # Initialize the structures of the many body state
    mbs = ManyBodyState(nStates, vsList)
    weight = np.multiply(vec,vec.conjugate())
    for stateIndex, state in enumerate(stateTable):
        stateRange = mbs.setMBS(state)
        nf = 0
        pol = 7
        for sbs in mbs.occList:
            occ[sbs] = occ[sbs] + weight[stateIndex]
            if aoList[sbs].l == 3:
                nf = nf + 1
                if aoList[sbs].sz > 0:
                    pol = pol + 1
                else:
                    pol = pol - 1
        fOcc[nf] = fOcc[nf] + weight[stateIndex]
        polHist[pol] = polHist[pol] + weight[stateIndex]

    print '========================================'
    print 'Occupation'
    for o in occ:
        print o
    print nParticles,'?=?', sum(occ)
    print '========================================'
    print 'F Occ Probability'
    for fo in fOcc:
        print fo
    print 'Sum probability:', sum(fOcc) 
    print '========================================'
    print 'Polarization Prob'
    for pl in polHist:
        print pl
    print 'Sum probability:', sum(polHist)
    print '========================================'

    polar = 0
    for ao in aoList:
        if abs(occ[ao.aoIndex]) > 1e-10:
            print ao.orb, 'spin', ao.sz, 'occ=', occ[ao.aoIndex] 
            if ao.sz < 0:
                polar = polar - occ[ao.aoIndex]
            else:
                polar = polar + occ[ao.aoIndex]
    print 'Polarization = ', polar
    print '========================================'

# <======================================================================>
# ----------------------------------------------------------------------
# 
# recursiveCalcExpList - recursive routine to expand a many-body state (mbs)
# ... into another basis set by processing each single-particle state (sps)
#
# result is a list of mbs and coefficients
# Useful for debugging 
# ----------------------------------------------------------------------
# <======================================================================>
def recursiveCalcExpList(aoList, state, nStates, nParticlesLeft, lastState):

    # ----------------------------------------------------------------------
    # Check if we are at the end of the list ... 
    # ... initialize the expansion list
    # ----------------------------------------------------------------------
    if nParticlesLeft == 0:
        newExpList = []
        newExpList.append((0, 1.0)) # tuple is the new mbs and 
        return newExpList

    # ----------------------------------------------------------------------
    # loop from the current sps lower
    # ----------------------------------------------------------------------
    for st in range(lastState,-1,-1):

        if state.get(st): # if occupied

            expList = recursiveCalcExpList(aoList, state, nStates, nParticlesLeft-1, st-1)

            ao = aoList[st]
            newExpList = []

            
            # ----------------------------------------------------------------------
            # expand this sps into the other basis (switchExpGen takes care of that)
            # ----------------------------------------------------------------------
            for eao in ao.switchExpGen():

                # Resuse the aoList to find the appropriate sps in the other basis
                # for example, p:x will return l=1, lz=1 and vice versa.
                newAO = aoList.get(ao.atom, ao.fragNo, ao.orb, ao.orbID, eao.ao.l, eao.ao.lz, ao.sz)
                if newAO is None:
                    print ao.atom, ao.fragNo, ao.orb, ao.orbID, eao.ao.l, eao.ao.lz, ao.sz
                    raise Exception('!!! ERROR !!! Null atomic orbital', )

                # for each sps in the expansion, 
                for exp in expList:
                    newState = ManyBodyState(nStates, exp[0], np.ones(nStates))
                    if not newState.get(newAO.aoIndex):
                        parity = newState.setWithParity(newAO.aoIndex)
                        newC = exp[1]*eao.c*parity
                        newExpList.append((newState.value, newC))
            return newExpList

# ----------------------------------------------------------------------
# Generate a transformation matrix between one many body state basis to another
# ... very useful for debugging if the basis sets don't get the same results
# ----------------------------------------------------------------------
def genTransformationMatrix():
    # Read in the state table pickles
    f = open('stateTable.pkl', 'r')
    stateTable = pickle.load(f)
    f.close()
    f = open('stateTableDict.pkl', 'r')
    stateTableDict = pickle.load(f)
    f.close()
    f = open('mbDimension.pkl', 'r')
    mbDim = pickle.load(f)
    f.close()

    plist = Plist.AtomicOrbitalPlist()
    plist.readPlist()
    aoList = plist.aoList
    nStates = plist.aoList.nStates
    nParticles = plist.nParticles()

    # loop through each state
    
    U = np.matrix(np.zeros((mbDim, mbDim), complex))
    for oldStateIndex, oldState in enumerate(stateTable):
        state = Bit.BitStr(nStates, oldState)
        print  'StateIndex %d value %s (%d)' % (oldStateIndex, state.display(), state.value)
        expList = recursiveCalcExpList(aoList, state, nStates, nParticles, nStates-1)
        for newState, coeff in expList:
            newStateBitStr = Bit.BitStr(nStates, newState)
            print newStateBitStr.display(), coeff
            newStateIndex = stateTableDict[newState]
            print 'New State', newStateIndex
            U[oldStateIndex, newStateIndex] = U[oldStateIndex, newStateIndex] + coeff
            
    if False:
        for r in range(mbDim):
            for c in range(mbDim):
                if U[r,c] != 0:
                    print 'U[', r,c, '] =', U[r,c]
        Ut = U.transpose().conjugate()
        for r in range(mbDim):
            for c in range(mbDim):
                if Ut[r,c] != 0:
                    print 'Ut[', r,c, '] =', Ut[r,c]
    return U

# ------------------------------------------------------------
# genFermionStateTable - generates the many body state table from the
# number of fermions and the total number of sb states.  Generates the
# following pickle files:
# - stateTable ... the mb state table
# - stateTableDict ... a reverse look up python dictionary mbs -> qID
# - mbDimension  ... the total number of many body states
# - nParticles ... the total number of particles in this model
# ------------------------------------------------------------
def genFermionStateTable(nFermions):

    plist = Plist.AtomicOrbitalPlist()
    plist.readPlist()
    nStates = plist.aoList.nStates

    debugFlag = True

    # In real situations, more than 32 states will overflow memory
    # But in small situations, this is allowed
    if nStates > 32:
        print '*** WARNING: *** nStates is larger than 32!! Continuing ... '
        

    # The number of fermions cannot be more than the number of states
    if nFermions > nStates:
        raise Exception('Number of particles exceeds the number of states! '+str(nFermions)+' '+str(nStates))
        return

    # The dimension is the total number of many-body states
    dim = CG.binomialCoeff(nStates, nFermions)
    print '<------------------------ State Table -------------------->'
    print 'Fermions = ', nFermions, ' ... States = ', nStates, '... Dimension = ', dim

    # Create an empty python dictionary
    stateTableDict = dict()
    stateTable = []

    # initialize the first state 
    # the low nFermions bits are set to 1
    state = Bit.BitStr(nStates)
    for ferm in range(nFermions):
        state.set(ferm)
        
    # cursor is the position of the highest bit in this group of ones
    cursor = nFermions - 1
    # ones is the total number of ones in this group of ones
    ones = nFermions

    # --------------------------------------------------------------------
    # Loop through each many-body state
    # --------------------------------------------------------------------
    for mbState in range(dim):

        stateTableDict[state.value] = mbState # set the key & value
        stateTable.append(state.value)
        if dim == 1:
            break

        # move the high bit up
        state.clear(cursor)
        cursor = cursor
        state.set(cursor+1)
        ones = ones - 1 # decrement the number of ones in the group
        
        if ones > 0:
            # --------------------------------------------------------------------
            # move all of the lower ones down
            # --------------------------------------------------------------------
            for sbState in range(cursor):
                if sbState < ones:
                    state.set(sbState)
                else:
                    state.clear(sbState)
            # --------------------------------------------------------------------
            cursor = ones - 1  # move the cursor to end of group

        else: # no more ones in the group
            # --------------------------------------------------------------------
            # find the next one to move
            for sbState in range(cursor+1,nStates):
                if not state.get(sbState):
                    cursor = sbState - 1
                    break
                else:
                    ones = ones + 1
            # --------------------------------------------------------------------
            # end for sbState in range(cursor+1,nStates+1):
            # --------------------------------------------------------------------

    f = open('./stateTable.pkl', 'w')
    pickle.dump(stateTable, f)
    f.close()

    f = open('./stateTableDict.pkl','w')
    pickle.dump(stateTableDict, f)
    f.close()

    del stateTable
    del stateTableDict

    # --------------------------------------------------------------------
    # Dump out the table for debugging
    # --------------------------------------------------------------------
    if debugFlag:
        # Read in the state table pickles
        f = open('stateTable.pkl', 'r')
        dbgStateTable = pickle.load(f)
        f.close()

        f = open('stateTableDict.pkl', 'r')
        dbgStateTableDict = pickle.load(f)
        f.close()

        for rStateIndex, rState in enumerate(dbgStateTable):
            state = Bit.BitStr(nStates, rState)
            print  'State %d = %s' % (rStateIndex, state.display())
            if rStateIndex != dbgStateTableDict[rState]:
                raise Exception('Mismatch between state index and dictionary!'+str(rStateIndex))
                
            

    f = open('./mbDimension.pkl','w')
    pickle.dump(dim, f)
    f.close()

    f = open('./nParticles.pkl','w')
    pickle.dump(nFermions, f)
    f.close()

    print '<------------------------ End State Table -------------------->'
    print
# <======================================================================>
# Allows calling from the command line:
if __name__ == '__main__':
    import sys
    genFermionStateTable()
# <======================================================================>
