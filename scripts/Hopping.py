import cPickle as pickle
import ManyBodyHam as mbh
import Coulomb


HoppingHam = mbh.ManyBodyHamMTX('Hopping')

f = open('Hopping.pkl', 'r')
H = pickle.load(f)
f.close()

F0 = 0.0
F2 = 0.0
F4 = 0.0
F6 = 0.0

if F2 != 0 or F4 != 0 or F6 != 0:
    scpTuplesList = [(3, [F0, 0.0, F2, 0.0, F4, 0.0, F6])]
    cmf = Coulomb.CoulombAndMFOperators(scpTuplesList)
    Hopping = H + cmf.Hmf
    cmf.Hcoul.processMatrixElements(HoppingHam)
else:
    Hopping = H


Hopping.processMatrixElements(HoppingHam)


HoppingHam.diag()
