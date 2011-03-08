import cPickle as pickle
import Coulomb
import ManyBodyHam as mbh

F0 = 1.0
F2 = 0.0
F4 = 0.0
F6 = 0.0

print '<------------------------------> F0 = ', F0

mbHam = mbh.ManyBodyHamMTX('Coulomb')

scpTuplesList = [(3, [F0, 0.0, F2, 0.0, F4, 0.0, F6])]
cmf = Coulomb.CoulombAndMFOperators(scpTuplesList)
cmf.processMatrixElements(mbHam)
mbHam.diag()
