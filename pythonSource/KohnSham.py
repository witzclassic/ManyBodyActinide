#!/usr/bin/env python

# Hks.py - Implementation of a Python interface to generate the
#         Kohn-Sham effective Hamiltonian.
#
# Copyright (C) 2009 Brown University Dept. of Physics
# Brad Marston, Steve Horowitz author
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

import cPickle as pickle
import Operators as op
import Plist 

# <======================================================================>
# Main routine - build the Kohn-Sham operator, inherits OneParticleOperator
# Fields: 
# - KohnShamOperatorDebugFlag: turns on debug printing
# Methods:
# Pickle Extracts:
# - Uo: unitary matrix rows are the atomic states, cols are the KS orbitals,
#  ... entries are the coefficients
# - Hdft: the energy eigenvalues for the KS orbitals
# - Nks: Number of electrons in the KS orbitals
# <======================================================================>
class KohnShamOperator(op.OneParticleOperator):
    """
    Class KohnShamOperator - A one particle operator extracted from the results of the DFT calculation
    and projected onto the atomic orbitals in the property list.
    """

    # ------------------------------------------------------------
    def __init__(self):
        self.KohnShamOperatorDebugFlag = False
        f = open('Uo.pkl', 'r')
        Uo = pickle.load(f)
        f.close()
        f = open('Hdft.pkl', 'r')
        Hdft = pickle.load(f)
        f.close()
        f = open('Nks.pkl', 'r')
        Nks = pickle.load(f)
        f.close()

        op.OneParticleOperator.__init__(self, Uo*Hdft*Uo.transpose().conjugate(), 'Hks')
        if self.KohnShamOperatorDebugFlag:
            self.dump()
        self.eigenvalues(Nks)
        pl = Plist.AtomicOrbitalPlist()
        pl.readPlist()        # read information from Ham.plist
        aoList = pl.aoList    # local copies of the atomic orbital list etc.
        for ao in aoList:
            print ao.fragNo, ao.orb, ao.orbID, self.op[ao.aoIndex, ao.aoIndex]
