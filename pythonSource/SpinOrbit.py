#!/usr/bin/env python

# SO.py - Spin Orbit
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
# Spin-Orbit Code
# ... contains the SpinOrbitOperator class, a one particle operator
# 
# 
# 
# <======================================================================>

import Operators
import Plist
import AO
import numpy
import cPickle as pickle

# <======================================================================>
# SpinOrbitOperator Class
# Inherits One particle operator
# Reads in spin orbit coupling constants (soCouplingConstants) from Property List 
# (Plist). Generates one particle operator, Hso = soXi L \cdot S + soAlpha
# soAlpha is a shift in the energy spectrum, soXi is the splitting 
# <======================================================================>

class SpinOrbitOperator(Operators.OneParticleOperator):
    """
    SpinOrbitOperator
    """
    # ------------------------------------------------------------
    def __init__(self):
        self.spinOrbitDebugFlag = False
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()
        aoList = pl.aoList        # atomic orbital list
        nStates = aoList.nStates  # number of states

        soCouplingConstants = pl.soCouplingConstants() # soZeta and soAlpha for each L

        soI = numpy.zeros((nStates, nStates), complex) # the spin-orbit interaction matrix

        # -------------------------------------------------------
        # for each spin-orbit coupling constant tuple, caclulate the interactions
        # between states in the j, jz basis:
        # -------------------------------------------------------
        for soCCTuple in soCouplingConstants:
            
            # extract l from the tuple
            soCCl = soCCTuple[Plist.SoCouplingEnum.l] 

            jp = soCCl + 0.5 # e.g. for f jp = 7/2
            jm = soCCl - 0.5 # and jm = 5/2

            # calculate the so coefficient j^2 - l^2 - s^2 
            # for example, for l=3, jpp = 0.5*((3.5*4.5) - 12 - 0.75) = 6.5
            # for l=1, jpp = 0.5*(15/4 - 8/4 - 3/4) = 1/2
            # jmm = 0.5*(3/4 - 8/4 - 3/4) = -1
            jpp = 0.5*(jp*(jp+1) - soCCl*(soCCl+1) - 0.75) 
            jmm = 0.5*(jm*(jm+1) - soCCl*(soCCl+1) - 0.75) 
            
            # calculate the so interaction operator for jp and jm
            jpOp = (soCCTuple[Plist.SoCouplingEnum.soXi] * jpp
                    + soCCTuple[Plist.SoCouplingEnum.soAlpha])
            jmOp = (soCCTuple[Plist.SoCouplingEnum.soXi] * jmm
                    + soCCTuple[Plist.SoCouplingEnum.soAlpha])

            print '<----------- Building Spin Orbit Operator...', 'l=', soCCl, \
                'xi =' ,soCCTuple[Plist.SoCouplingEnum.soXi],  \
                'alpha =' ,soCCTuple[Plist.SoCouplingEnum.soAlpha],  \
                'j+ = ', jpOp, 'j-=', jmOp, ') ------------>'

            # -------------------------------------------------------
            # loop through each ao ... find its j basis buddy and compute interaction
            # <j, jz | jOp | j, jz>
            # -------------------------------------------------------
            for ao1 in aoList.aoByLGen(soCCl):
                for jao1 in ao1.jExpGen():
                    j1 = jao1.ao.j
                    jz1 = jao1.ao.jz
                    for ao2 in aoList.aoByLGen(soCCl):
                        for jao2 in ao2.jExpGen():

                            j2 = jao2.ao.j
                            jz2 = jao2.ao.jz
                            # ------------------------------------------------------------
                            #  j(n) dot j(n) = j(j+l)
                            # ------------------------------------------------------------
                            if (jao1.ao.j == jao2.ao.j) and (jz1 == jz2):

                                if jao1.ao.j == jp:
                                    op = jpOp
                                else:
                                    op = jmOp

                                # I = I + op * expansion coefficients
                                soI[ao1.aoIndex, ao2.aoIndex] = \
                                    soI[ao1.aoIndex, ao2.aoIndex] + (op * jao1.c.conjugate() * jao2.c)


        Operators.OneParticleOperator.__init__(self, soI, 'Hso')
        if self.spinOrbitDebugFlag:
            self.dump()


