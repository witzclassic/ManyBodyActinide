#!/usr/bin/env python

# ManyBodyHam
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
# - ManyBodyHam base class of:
# -- ManyBodyHamMTX generates an MTX file
# -- ManyBodyHamMatrix generates a matrix (not good for large matrices due to memory reqs.)
# The following were attempts at optimization, but are not fully tested.
# -- ManyBodyHamSparseMatrix
# -- ManyBodyHamMatrixBinaryFile
# Also:
# - LowLyingStates: reads in eigenvector from LowLyingStates.mtx file
# - EigenStates: formats spectrum of diagonalization calculation
#
# <======================================================================>

import struct
from scipy import sparse
import numpy as np
import Bit
import cPickle as pickle
import Plist
import copy
import types

# <======================================================================>
# Low Lying States: reads in eigenvector from LowLyingStates.mtx
# Fields: 
# - v: the eigenvector
# Methods:
# - dump: print debug stuff
# Pickle loads: 
# - mbDimension: length of the eigenvector
# - LowLyingStates.mtx: eigenvector output of sparseMatrix diagonlization utility.
# .. format is a rough market matrix format, entries are complex format coeffs.
# <======================================================================>
class LowLyingStates:
    """
    Low Lying States
    """

    # ------------------------------------------------------------
    def __init__(self):
        f = open('mbDimension.pkl', 'r')
        mbDim = pickle.load(f)
        f.close()

        f = open('LowLyingStates.mtx', 'r')
        header = f.readline()
        entries, gsEnergy = header.split()
        vList = np.zeros(mbDim, complex)
        
        for index, line in enumerate(f):
            real, imag = line.split()
            freal = float(real)
            fimag = float(imag)
            if freal < 1e-10 and freal > -1e-10:
                freal = 0
            if fimag < 1e-10 and fimag > -1e-10:
                fimag = 0
            vList[index] = complex(float(real), float(imag))
            
        f.close()
        self.v = np.matrix(map(lambda x: [x], vList))

    # ------------------------------------------------------------
    def dump(self):
        for vec in self.v:
            print vec
        
# <======================================================================>
# EigenStates: Used for printing spectrum from diagonlization
# Fields: 
# - e: the energy eigenvalues list
# - v: the coeffs of the many-body state, the eigenvector list
# Methods:
# - spectrum: generate the spectrum printing
# ... optional list of manyBodyHams for computing expectation values.
# <======================================================================>
class EigenStates:
    """
    Eigenstates
    """

    # ------------------------------------------------------------
    def __init__(self, e,v):
        self.e = e
        self.v = np.matrix(v)

    # --------------------------------------------------
    # spectrum ... compute the energy spectrum
    # ... optional list of manyBodyHams for computing expectation values.
    # --------------------------------------------------
    def spectrum(self, mbHamList=[]):

        print '----------------------------------------'
        print 'Spectrum processing'
        print '----------------------------------------'

        numExpValues = len(mbHamList)

        base = self.e[0]
        last = self.e[0]
        mult = 0
        uniqueIndex = 1

        expValue = np.zeros((numExpValues, len(self.e)), float)
        for i in range(numExpValues):
            expValue[i] = np.diagonal(self.v.transpose() * mbHamList[i].mbh * self.v.conjugate())

        # --------------------------------------------------
        # loop thru each eigenvalue
        # --------------------------------------------------
        for eIndex, energy in enumerate(self.e):

            print
            if (energy - last) > 1e-10:
                print '----------------------------------------'
                print 'Diff' , last-base, '(', mult, ')'
                print uniqueIndex, round(last, 3)
                print '----------------------------------------'
                print
                last = energy
                mult = 1
                uniqueIndex = uniqueIndex + 1
            else:
                mult = mult + 1
            print 'Eigenvalue', eIndex, round(energy, 3)
            for i in range(numExpValues):
                print mbHamList[i].name, 'expValue', i, ' = ', round(expValue[i, eIndex], 2)
            #for vecIndex, vec in enumerate(self.v[:,eIndex]):
                #if abs(vec) > 1e-10:
                    #print 'State: ', vecIndex, '=', vec
        print
        print '----------------------------------------'
        print 'Diff' , last-base, '(', mult, ')'
        print uniqueIndex, last
        print '----------------------------------------'
        print

# <======================================================================>
# ManyBodyHam: The base class for all many body ham classes
# Pickle Load: 
# - mbDimension: the many body dimensions
# <======================================================================>
class ManyBodyHam:
    """
    Many Body Ham Base Class
    """
    # ------------------------------------------------------------
    def __init__(self):
        f = open('mbDimension.pkl', 'r')
        self.mbDim = pickle.load(f)
        f.close()

# <======================================================================>
# ManyBodyHamMTX - generates an MTX file
# Inherits ManyBodyHam
# Use this format for very large many body hams
# Header includes mbDim, mbDim, realCount, complexCount
# Fields: 
# - formats for writing entries to mtx file
# - fname: file name without .mtx suffix
# - entryCount: total number of entries
# - diagCount: total number of diagonal entries (row == col)
# - realCount: total number of real entries
# - complexCount: total number of complex entries
# - open, read flags: supports case when file is already opened, and read/write
# Methods:
# - diag: Does not diagonalize .. only writes out mtx file.  Use sparseMatrix to diag.
# - open: Opens file
# - close 
# - setDiagMatrixElement: writes out diagonal entry
# - setMatrixElement: writes out entry for upper half triangle (lState > rState)
# - getFile: returns self.file
# - readline: calls file.readline()
# - readHeader: skips first 3 lines
# - readEntry: reads next entry in open file, returns [row, col, real, cmplx]
# - diagInPlace: Debug routine to read in mtx and diagonalize using ManyBodyHamMatrix 
# ... do not use for large matrices
# - expValue: compute the expectation value of supplied state.  Obsolete Routine.
# - dump .. not implemented
# <======================================================================>
class ManyBodyHamMTX(ManyBodyHam):
    """
    Many Body Ham in modified Market Matrix (mtx) format.
    """

    # ------------------------------------------------------------
    def __init__(self, fname, readFlag=False):
        ManyBodyHam.__init__(self)

        self.hdrFormat = '%15d %15d %15d %d\n'
        self.realEntryFormat = '%15d %15d %.15f %d\n'
        self.imagEntryFormat = '%15d %15d %.15f %.15f\n'

        self.fname = fname
        self.entryCount = 0
        self.openFlag = False
        self.readFlag = False

        self.open(readFlag)

    # ------------------------------------------------------------
    def open(self, readFlag):

        if self.readFlag != readFlag and self.openFlag == True:
            self.close()  # open in a different mode

        if self.openFlag:
            return  # already open

        self.readFlag=readFlag

        if readFlag:
            self.file = open(self.fname+'.mtx', 'r')
            self.file.readline()
            self.file.readline()
            self.file.readline()
        else:
            self.file = open(self.fname+'.mtx', 'w')
            self.file.write('%%MatrixMarket matrix coordinate complex hermitian\n')
            self.file.write('%rows, cols, nz entries real and imag\n')
            self.savePos = self.file.tell()
            self.file.write(self.hdrFormat % (0,0,0,0))
            self.complexCount = 0
            self.realCount = 0
            self.diagCount = 0

        self.openFlag = True

            
    # ------------------------------------------------------------
    def setDiagMatrixElement(self, stateIndex, value):
        self.file.write(self.realEntryFormat % (stateIndex+1, stateIndex+1, value, 0))
        self.diagCount = self.diagCount + 1
        
    # ------------------------------------------------------------
    # Store only the upper half triangle entries 
    # ------------------------------------------------------------
    def setMatrixElement(self, lStateIndex, rStateIndex, value):
        
        if value.imag > 1e-10 or value.imag < -1e-10:
            self.complexCount = self.complexCount + 1
            self.file.write(self.imagEntryFormat % (rStateIndex+1, lStateIndex+1, value.real, -value.imag))
        elif value.real > 1e-10 or value.real < -1e-10:
            self.realCount = self.realCount + 1
            self.file.write(self.realEntryFormat % (rStateIndex+1, lStateIndex+1, value, 0))

    # ------------------------------------------------------------
    def getFile(self):
        return self.file
    
    # ------------------------------------------------------------
    def readLine(self):
        return self.file.readline()
    
    # ------------------------------------------------------------
    def readHeader(self):
        self.readLine()
        self.readLine()
        self.readLine()

    # ------------------------------------------------------------
    # readEntry: reads next entry in open file, returns [row, col, real, cmplx]
    # ------------------------------------------------------------
    def readEntry(self):
        line = self.readLine()
        if line == "":
            return [-1,-1,0]

        [rowTxt, colTxt, realTxt, imagTxt] = line.split()
        row = int(rowTxt)
        col = int(colTxt)
        entry = complex(float(realTxt), float(imagTxt))
        return int(rowTxt)-1, int(colTxt)-1, complex(float(realTxt), float(imagTxt))

    # ------------------------------------------------------------
    # - diagInPlace: Debug routine to read in mtx and diagonalize using ManyBodyHamMatrix 
    # ... do not use for large matrices
    # ------------------------------------------------------------
    def diagInPlace(self):
        readFlag = True
        self.open(readFlag)
        m = ManyBodyHamMatrix('self.name')
        while True:
            row, col, entry = self.readEntry()
            if row < 0: 
                e,v = m.diag()
                return e,v
            if row == col:
                m.setDiagMatrixElement(row,entry)
            else:
                m.setMatrixElement(row, col, entry)
            
        
    # ------------------------------------------------------------
    # Expectation Value of this mbHam and the supplied state vector
    # ... Obsololete routine ...
    # ------------------------------------------------------------
    def expValue(self, vec):
        readFlag = True
        self.open(readFlag)
        exp = 0
        while True:
            row, col, entry = self.readEntry()
            if row < 0: 
                return round(exp,2)
            exp = exp + vec[row].conjugate()*entry*vec[col]
            
            if row != col:
                exp = exp + vec[col].conjugate()*entry.conjugate()*vec[row]

    # ------------------------------------------------------------
    def dump(self):
        print 'dump not implemented yet.'
        return

    # ------------------------------------------------------------
    def close(self):
        if self.openFlag == False:
            return

        self.file.close()
        self.openFlag = False

    # ------------------------------------------------------------
    # MTX doesn't diag ... just writes the file
    # ------------------------------------------------------------
    def diag(self):
        print 'Wrote out mtx file: ', self.mbDim, self.realCount, self.complexCount
        self.file.seek(self.savePos)
        self.file.write(self.hdrFormat % (self.mbDim, self.mbDim,self.realCount, self.complexCount))
        self.close()

# <======================================================================>
# Class ManyBodyHamSparseMatrix
# Integration with sparseMatrix diag routine is incomplete, use with caution.
# Inherits ManyBodyHam
# Fields: 
# mbhd, mbhc, mbhr: lists of entries assoc with diag, cmplx and real entries.
# name: text name
# Methods:
# - setDiagMatrixElement: write to diag array
# - setMatrixElement: appends entry to approp. list, row, col, value
# - dump: prints total number of entries
# - genBinary: writes out in binary format
# - diag: does not diag ... calls genBinary.  Use sparseMatrix to diag.
# <======================================================================>
class ManyBodyHamSparseMatrix(ManyBodyHam):
    """
    Many Body Ham Sparse Matrix format
    Uses lists.  Supports a dump directly to the binary file format.
    """
    # ------------------------------------------------------------
    def __init__(self, name):
        ManyBodyHam.__init__(self)
        self.mbhd = np.zeros((self.mbDim))
        self.mbhc = []
        self.mbhr = []

        self.name = name
        
    # ------------------------------------------------------------
    def setDiagMatrixElement(self, stateIndex, value):
        self.mbhd[stateIndex] = self.mbhd[stateIndex] + value.real

    # ------------------------------------------------------------
    def setMatrixElement(self, lStateIndex, rStateIndex, value):
        if value.imag:
            self.mbhc.append((lStateIndex, rStateIndex, value.real, value.imag))
        else:
            self.mbhr.append((lStateIndex, rStateIndex, value))


    # ------------------------------------------------------------
    def dump(self):
        print 'mbHam: ', self.name
        print 'Real:', self.mbhr.nnz
        print 'Complex:', self.mbhc.nnz

    # ------------------------------------------------------------
    def genBinary(self):
        numReal = len(self.mbhr)
        numComplex = len(self.mbhc)

        print 'Generating ', self.name+'.bin', self.mbDim, numReal, numComplex
        f = file(self.name+'.bin', 'w')
        f.write(struct.pack('@I', self.mbDim))
        f.write(struct.pack('@I', self.mbDim))
        f.write(struct.pack('@I', numReal))
        f.write(struct.pack('@I', numComplex))

        for d in range(self.mbDim):
            f.write(struct.pack('@f', self.mbhd[d]))

        for entry in self.mbhr:
            f.write(struct.pack('LLf', entry[0], entry[1], entry[2]))
            
        for entry in self.mbhc:
            f.write(struct.pack('LLff', entry[0], entry[1], entry[2], entry[3]))

            

    # ------------------------------------------------------------
    def diag(self):
        print '<------------------------ Diag (', self.name, ') --------------------------->'
        self.genBinary()


# <======================================================================>
# Class ManyBodyHamMatrix - A full n x n matrix
# Inherits ManyBodyHam
# Fields:
# - mbh : many body ham in numpy matrix format
# - name: text name
# - real : real entry count
# - cpx : complex entry count
# - diags : diagonal entry count
# Methods:
# - setDiagMatrixElement: add value to diagonal entry
# - setMatrixElement: mbh[r,c] = value, mbh[c,r] = value.conjugate()
# - dump
# - diag: calls numpy linalg.eigh and uses EigenStates to generate spectrum
# - expValue: generates expectation value from supplied vector state
# <======================================================================>
class ManyBodyHamMatrix(ManyBodyHam):

    # ------------------------------------------------------------
    def __init__(self, name):
        ManyBodyHam.__init__(self)
        self.mbh = np.matrix(np.zeros((self.mbDim, self.mbDim), complex))
        self.name = name
        self.real = 0
        self.cpx = 0
        self.diags = 0

    # ------------------------------------------------------------
    def setDiagMatrixElement(self, stateIndex, value):
        self.mbh[stateIndex,stateIndex] = self.mbh[stateIndex,stateIndex] + value
        self.diags = self.diags + 1

    # ------------------------------------------------------------
    def setMatrixElement(self, lStateIndex, rStateIndex, value):
        if value.imag > 1e-10 or value.imag < -1e-10:
            self.mbh[lStateIndex,rStateIndex] = self.mbh[lStateIndex,rStateIndex] + value
            self.mbh[rStateIndex,lStateIndex] = self.mbh[lStateIndex,rStateIndex].conjugate()
            self.cpx = self.cpx + 1
        else:
            self.mbh[lStateIndex,rStateIndex] = self.mbh[lStateIndex,rStateIndex] + value.real
            self.mbh[rStateIndex,lStateIndex] = self.mbh[lStateIndex,rStateIndex].conjugate()
            self.real = self.real + 1

            

    # ------------------------------------------------------------
    def dump(self):
        print 'mbHam: ', self.name
        for row in range(self.mbDim):
            for col in range(self.mbDim):
                if abs(self.mbh[row,col]) > 1e-10:
                    print '[', row, ',', col,']', '=',self.mbh[row,col]
#        print self.mbh

    # ------------------------------------------------------------
    def diag(self, mbHamList=[]):
        print '<------------------------ Diag (', self.name, ') --------------------------->'
        print self.diags, self.real, self.cpx
        print
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()
        aoList = pl.aoList
        nStates = aoList.nStates
        e,v = np.linalg.eigh(self.mbh)
        print self.name, 'Eigenvalues', np.round(e, 10)
        #print v
        print
        eStates = EigenStates(e,v)
        eStates.spectrum(mbHamList)

        return (e,v)

    # ------------------------------------------------------------
    def expValue(self, vec):
        v = np.matrix(vec)
        print v.shape
        return v* self.mbh*v.H

# <======================================================================>
# Class ManyBodyHamMatrixBinaryFile - Generate optimized binary file
# Attempt to use sparse.lil_matrix but probably does not work.
# <======================================================================>
class ManyBodyHamMatrixBinaryFile(ManyBodyHam):
    def __init__(self, name):
        ManyBodyHam.__init__(self)
        self.mbhd = np.zeros((self.mbDim))
        self.mbhc = sparse.lil_matrix((self.mbDim, self.mbDim), dtype=complex)
        self.mbhr = sparse.lil_matrix((self.mbDim, self.mbDim), dtype=float)

        self.name = name
        
    def setMatrixElement(self, lStateIndex, rStateIndex, value):
        if abs(value) > 1e-10:
            if lStateIndex == rStateIndex:
                self.mbhd[lStateIndex] = self.mbhd[lStateIndex] + value.real
            else :
                if abs(value.imag) > 1e-10:
                    self.mbhc[lStateIndex,rStateIndex] = self.mbhc[lStateIndex,rStateIndex] + value
                else:
                    self.mbhr[lStateIndex,rStateIndex] = self.mbhr[lStateIndex,rStateIndex] + value

    # ------------------------------------------------------------
    def dump(self):
        print 'mbHam: ', self.name
        print 'Real:', self.mbhr.nnz
        print 'Complex:', self.mbhc.nnz

    # ------------------------------------------------------------
    def genBinary(self):
        rrList, rcList = self.mbhr.nonzero()
        crList, ccList = self.mbhc.nonzero()
        numReal = len(rrList)
        numComplex = len(crList)

        print 'Generating ', self.name+'.bin', self.mbDim, numReal, numComplex
        f = file(self.name+'.bin', 'w')
        f.write(struct.pack('@I', self.mbDim))
        f.write(struct.pack('@I', self.mbDim))
        f.write(struct.pack('@I', numReal))
        f.write(struct.pack('@I', numComplex))

        for d in range(self.mbDim):
            f.write(struct.pack('@f', self.mbhd[d]))

        for row, col in zip(rrList, rcList):
            if row != col:
                f.write(struct.pack('LLf', row.__int__(), col.__int__(), self.mbhr[row, col]))

        for row, col in zip(crList, ccList):
            if row != col:
                f.write(struct.pack('LLff', row.__int__(), col.__int__(), self.mbhc[row, col].real, self.mbhc[row,col].imag))

            

    # ------------------------------------------------------------
    def diag(self):
        print '<------------------------ Diag (', self.name, ') --------------------------->'
        self.genBinary()




