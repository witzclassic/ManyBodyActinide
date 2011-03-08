#!/usr/bin/env python

# AngMo - Angular Momentum
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
# Angular Momentum operators
# S^2
# L^2
# J^2
# 
# A^2 = sum_i A_i^2 + sum_{i \ne j} (A_i^+ A_j^- + A_i^- A_j^+) + 2A_i^zA_j^z
# 
# <======================================================================>
# Interaction Tensor: 4 dimensional structure
# AngMo2Op: 
# LzOp, JzOp: Lz, Jz operator
# J2Op, L2Op, S2Op: Ang Mom. operators
# Routine buildOperators which builds A2Op (A=J,L,S) operators.
# <======================================================================>
import numpy as np
from math import sqrt

import Plist
import ManyBodyHam as mbh
import Operators

class InteractionTensor:
    def __init__(self, nStates, type=complex):
        self.I = np.zeros((nStates,nStates,nStates,nStates), type)
        self.nStates = nStates

    def add(self, i1,i2,i3,i4,value):
        self.I[i1,i2,i3,i4] = self.I[i1,i2,i3,i4] + value

    def cleanup(self):
        for i1 in xrange(self.nStates):
            for i2 in xrange(self.nStates):
                for i3 in xrange(self.nStates):
                    for i4 in xrange(self.nStates):
                        if (self.I[i1, i2, i3, i4] != 0) and (abs(self.I[i1, i2, i3, i4].real) < 1e-14 ):
                            self.I[i1, i2, i3, i4] = 0


class AngMo2Op():
    def __init__(self, L):
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()
        self.aoList = pl.aoList
        self.nStates = self.aoList.nStates

        self.vsList = np.zeros((self.nStates), int)
        for ao in self.aoList.aoByLGen(L):
            self.vsList[ao.aoIndex] = 1
            
        
        # class AngMoExp ... l, lz, c representing ang mo, ang mo z component and coeff

        self.I = InteractionTensor(self.nStates, complex)
        self.T = np.zeros((self.nStates, self.nStates), complex)


class JzOp(Operators.OneParticleOperator):
    def __init__(self, L):
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()
        aoList = pl.aoList
        nStates = aoList.nStates

        T = np.zeros((nStates, nStates), complex)
        for ao1 in aoList.aoByLGen(L):
            for jao1 in ao1.jExpGen():
                j1 = jao1.ao.j  
                jz1 = jao1.ao.jz
                for ao2 in aoList.aoByLGen(L):
                    for jao2 in ao2.jExpGen():
                        j2 = jao2.ao.j
                        jz2 = jao2.ao.jz
                        
                        if (jz1 == jz2 and j1 == j2):

                            coeff = jao1.c.conjugate()*jao2.c
                            T[ao1.aoIndex, ao2.aoIndex] = T[ao1.aoIndex, ao2.aoIndex] + (jz1 * coeff)
                            #print ao1.aoIndex, ao2.aoIndex, jz1 * coeff

        Operators.OneParticleOperator.__init__(self, T, 'Jz')

class LzOp(Operators.OneParticleOperator):
    def __init__(self, L):
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()
        aoList = pl.aoList
        nStates = aoList.nStates

        T = np.zeros((nStates, nStates))
        for ao1 in aoList.aoByLGen(L):
            for sao1 in ao1.sphericalExpGen():
                l1 = sao1.ao.l
                lz1 = sao1.ao.lz

                for ao2 in aoList.aoByLGen(L):
                    for sao2 in ao2.sphericalExpGen():
                        l2 = sao2.ao.l
                        lz2 = sao2.ao.lz
                        
                        coeff = sao1.c.conjugate()*sao2.c
                        
                        if (lz1 == lz2):
                            T[ao1.aoIndex, ao2.aoIndex] = T[ao1.aoIndex, ao2.aoIndex] + (lz1 * coeff)
                            
        Operators.OneParticleOperator.__init__(self, T, 'Jz')

        
        
# <======================================================================>
# J2Op --- J^2
# <======================================================================>
class J2Op(AngMo2Op):
    def __init__(self, L):
        AngMo2Op.__init__(self, L)

        for ao1 in self.aoList.aoByLGen(L):
            for jao1 in ao1.jExpGen():
                j1 = jao1.ao.j  
                jz1 = jao1.ao.jz
                for ao2 in self.aoList.aoByLGen(L):
                    for jao2 in ao2.jExpGen():
                        j2 = jao2.ao.j
                        jz2 = jao2.ao.jz

                        if (ao1 is not ao2) and (jao1 is not jao2): 

                            for ao3 in self.aoList.aoByLGen(L):
                                for jao3 in ao3.jExpGen():
                                    j3 = jao3.ao.j
                                    jz3 = jao3.ao.jz

                                    if (j1 == j3):

                                        for ao4 in self.aoList.aoByLGen(L):
                                            for jao4 in ao4.jExpGen():
                                                j4 = jao4.ao.j
                                                jz4 = jao4.ao.jz
                                                if (ao3 is not ao4) and (jao3 is not jao4) and (j2 == j4):
                                                
                                                    coeff = jao1.c.conjugate()*jao2.c.conjugate()*jao3.c*jao4.c

                                                    # ------------------------------------------------------------
                                                    # j(1)+ j(2)- = 
                                                    #     sqrt[(j-jz(1))(j+jz+1)] * sqrt[(j+jz(2))(j-jz(2)+1)
                                                    #         jz(1) -> jz(1)+1 , jz(2) -> jz(2) - 1
                                                    # ------------------------------------------------------------
                                                    if (jz1 == jz3 + 1) and (jz2 == jz4 - 1):
                                                    
                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                       sqrt( (j3*(j3+1)-jz3*(jz3+1)) * (j4*(j4+1) - jz4*(jz4-1))) * coeff)

                                                    # ------------------------------------------------------------
                                                    #  j(1)- j(2)+ = 
                                                    #     sqrt[(j+jz(1))(j-jz(1)+1)] * sqrt[(j-jz(2))(j+jz(2)+1)
                                                    #         jz(1) -> jz(1) - 1 , jz(2) -> jz(2) + 1
                                                    # ------------------------------------------------------------
                                                    if (jz1 == jz3 - 1) and (jz2 == jz4 + 1):

                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                       sqrt( (j3*(j3+1)-jz3*(jz3-1)) * (j4*(j4+1) - jz4*(jz4+1))) * coeff)

                                                    # ------------------------------------------------------------
                                                    #  2*j(1)z j(2)z 
                                                    # ------------------------------------------------------------
                                                    if (jz1 == jz3 and jz2 == jz4):

                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                  2*jz1*jz2 * coeff)



        self.exchangeOp = Operators.TwoParticleOperator(self.I.I, 'J2ExchangeOp')
        self.exchangeOp.setValidStates(self.vsList)
        #self.exchangeOp.setOptimized()
        #self.exchangeOp.dump()

        for ao1 in self.aoList.aoByLGen(L):
            for jao1 in ao1.jExpGen():
                j1 = jao1.ao.j
                jz1 = jao1.ao.jz
                for ao2 in self.aoList.aoByLGen(L):
                    for jao2 in ao2.jExpGen():

                        j2 = jao2.ao.j
                        jz2 = jao2.ao.jz
                        # ------------------------------------------------------------
                        #  j(n) dot j(n) = j(j+l)
                        # ------------------------------------------------------------
                        if (jao1.ao.j == jao2.ao.j) and (jz1 == jz2):
                            self.T[ao1.aoIndex, ao2.aoIndex] = self.T[ao1.aoIndex, ao2.aoIndex] + \
                                (jao1.ao.j + 1)*jao1.ao.j * jao1.c.conjugate()*jao2.c

        self.directOp = Operators.OneParticleOperator(self.T, 'J2DirectOp')
        self.directOp.setValidStates(self.vsList)
        #self.directOp.dump()
        
    def processMatrixElements(self, mbHam):
        self.exchangeOp.processMatrixElements(mbHam)
        self.directOp.processMatrixElements(mbHam)
        return mbHam


    def dump(self):
        self.directOp.dump()
        self.exchangeOp.dump()

# <======================================================================>
# L2Op --- L^2
# <======================================================================>
class L2Op(AngMo2Op):
    def __init__(self, L):
        AngMo2Op.__init__(self, L)

        # sao is class AngMoExp ... l, lz, c representing ang mo, ang mo z component and coeff

        for ao1 in self.aoList.aoByLGen(L):
            for sao1 in ao1.sphericalExpGen():
                l1 = sao1.ao.l
                lz1 = sao1.ao.lz
                c1 = sao1.c.conjugate()

                for ao2 in self.aoList.aoByLGen(L):
                    for sao2 in ao2.sphericalExpGen():
                        l2 = sao2.ao.l
                        lz2 = sao2.ao.lz
                        c2 = sao2.c.conjugate()

                        if (ao1 is not ao2) and (sao1 is not sao2) :

                            for ao3 in self.aoList.aoByLGen(L):
                                for sao3 in ao3.sphericalExpGen():
                                    l3 = sao3.ao.l
                                    lz3 = sao3.ao.lz
                                    c3 = sao3.c

                                    # l and sz is conserved .. already looping by L
                                    if (sao1.ao.sz == sao3.ao.sz):

                                        for ao4 in self.aoList.aoByLGen(L):
                                            for sao4 in ao4.sphericalExpGen():
                                                l4 = sao4.ao.l
                                                lz4 = sao4.ao.lz
                                                c4 = sao4.c

                                                coeff = c1*c2*c3*c4

                                                # l and sz is conserved ... note we're looping by L
                                                # double destructor
                                                if (ao3 is not ao4) and (sao3 is not sao4) and (sao2.ao.sz == sao4.ao.sz): 

                                                    # ------------------------------------------------------------
                                                    # l(1)+ l(2)- = 
                                                    #     sqrt[l(l+1) - lz(1)(lz(1)+1)] *sqrt[l(l+1) - lz(1)(lz(1)-1)]
                                                    #         lz(1) -> lz(1)+1 , lz(2) -> lz(2) - 1
                                                    # ------------------------------------------------------------
                                                    if (lz1 == lz3 + 1) and (lz2 == lz4 - 1):
                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                       sqrt( (l3*(l3+1) - lz3*(lz3 + 1)) * (l4*(l4+1) - lz4*(lz4-1))) * coeff)

                                                    # ------------------------------------------------------------
                                                    #  l(1)- l(2)+ = 
                                                    #     sqrt[l(l+1) - lz(1)(lz(1)-1)] *sqrt[l(l+1) - lz(1)(lz(1)+1)]
                                                    #         lz(1) -> lz(1) - 1 , lz(2) -> lz(2) + 1
                                                    # ------------------------------------------------------------
                                                    if (lz1 == lz3 - 1) and (lz2 == lz4 + 1):
                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                       sqrt( (l3*(l3+1) - lz3*(lz3 - 1)) * (l4*(l4+1) - lz4*(lz4 + 1))) * coeff)
                                            
                                                    # ------------------------------------------------------------
                                                    #  2*l(1)z l(2)z 
                                                    # ------------------------------------------------------------
                                                    if (lz1 == lz3) and (lz2 == lz4):
                                                        self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                                                       2 * lz1 * lz2 * coeff)


        self.exchangeOp = Operators.TwoParticleOperator(self.I.I, 'L2ExchangeOp')
        self.exchangeOp.setValidStates(self.vsList)
        self.exchangeOp.setOptimized()

        for ao1 in self.aoList.aoByLGen(L):
            for sao1 in ao1.sphericalExpGen():

                for ao2 in self.aoList.aoByLGen(L):
                    for sao2 in ao2.sphericalExpGen():

                        # ------------------------------------------------------------
                        #  l(n) dot l(n) = l(l+l)
                        # ------------------------------------------------------------
                        if (sao1.ao.lz == sao2.ao.lz and ao1.sz == ao2.sz):
                            self.T[ao1.aoIndex, ao2.aoIndex] = self.T[ao1.aoIndex, ao2.aoIndex] + \
                                (sao1.ao.l + 1)*sao1.ao.l*sao1.c.conjugate()*sao2.c

        self.directOp = Operators.OneParticleOperator(self.T, 'L2DirectOp')
        self.directOp.setValidStates(self.vsList)
        #self.directOp.dump()

        
    def processMatrixElements(self, mbHam):
        self.exchangeOp.processMatrixElements(mbHam)
        self.directOp.processMatrixElements(mbHam)
        #mbHam.dump()
        return mbHam

    def dump(self):
        self.directOp.dump()
        self.exchangeOp.dump()
        
# <======================================================================>
# S2Op --- S^2
# <======================================================================>
class S2Op(AngMo2Op):
    def __init__(self, L):
        AngMo2Op.__init__(self, L)

        s = 0.5

        for ao1 in self.aoList.aoByLGen(L):
            for ao2 in self.aoList.aoByLGen(L):
                for ao3 in self.aoList.aoByLGen(L):
                    for ao4 in self.aoList.aoByLGen(L):
                        # l and lz are not affected
                        if (ao1.lz == ao3.lz and ao2.lz == ao4.lz):

                            # ------------------------------------------------------------
                            # s(1)+ s(2)- = 
                            #     sqrt[(s-sz(1))(s+sz(1)+1)] * sqrt[(s+sz(2))(s-sz(2)+1)
                            #         sz(1) -> sz(1)+1 , sz(2) -> sz(2) - 1
                            # ------------------------------------------------------------
                            if (ao1.sz == ao3.sz + 1) and (ao2.sz == ao4.sz - 1):
                                                        
                                self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                          sqrt((s-ao3.sz)*(s+ao3.sz+1)) \
                                          * sqrt((s+ao4.sz)*(s-ao4.sz+1)) )

                            # ------------------------------------------------------------
                            #  s(1)- s(2)+ = 
                            #     sqrt[(1+sz(1))(1-sz(1)+1)] * sqrt[(1-sz(2))(1+sz(2)+1)
                            #         sz(1) -> sz(1) - 1 , sz(2) -> sz(2) + 1
                            # ------------------------------------------------------------
                            if (ao1.sz == ao3.sz - 1) and (ao2.sz == ao4.sz + 1):

                                self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                          sqrt((s+ao3.sz)*(s-ao3.sz+1))\
                                          * sqrt((s-ao4.sz)*(s+ao4.sz+1)))
                                            
                            # ------------------------------------------------------------
                            #  2*l(1)z l(2)z 
                            # ------------------------------------------------------------
                            if (ao1.sz == ao3.sz and ao2.sz == ao4.sz):
                                self.I.add(ao1.aoIndex, ao2.aoIndex, ao3.aoIndex, ao4.aoIndex, \
                                          2*ao1.sz*ao2.sz )


        self.exchangeOp = Operators.TwoParticleOperator(self.I.I, 'S2ExchangeOp')
        self.exchangeOp.setValidStates(self.vsList)
        self.exchangeOp.setOptimized()
        

        for ao1 in self.aoList.aoByLGen(L):
            for ao2 in self.aoList.aoByLGen(L) :
            
                # ------------------------------------------------------------
                #  s(n) dot s(n) = s(s+l)
                # ------------------------------------------------------------
                if (ao1.lz == ao2.lz and ao1.sz == ao2.sz):
                    self.T[ao1.aoIndex, ao2.aoIndex] = self.T[ao1.aoIndex, ao2.aoIndex] + (s + 1)*s

        self.directOp = Operators.OneParticleOperator(self.T, 'S2DirectOp')
        self.directOp.setValidStates(self.vsList)
        #self.directOp.dump()
        
    def processMatrixElements(self, mbHam):
        self.exchangeOp.processMatrixElements(mbHam)
        self.directOp.processMatrixElements(mbHam)
        #mbHam.dump()
        return mbHam

def buildOperators():
    L = 3
    s2Ham = mbh.ManyBodyHamMTX('SSquaredOp')
    s2op = S2Op(L)
    s2op.processMatrixElements(s2Ham)
    s2Ham.diag()
    l2Ham = mbh.ManyBodyHamMTX('LSquaredOp')
    l2op = L2Op(L)
    l2op.processMatrixElements(l2Ham)
    l2Ham.diag()
    j2Ham = mbh.ManyBodyHamMTX('JSquaredOp')
    j2op = J2Op(L)
    j2op.processMatrixElements(j2Ham)
    j2Ham.diag()

    
