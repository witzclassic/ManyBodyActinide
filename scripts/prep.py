import KohnSham
import SpinOrbit
import extractkf
import cPickle as pickle

f = open('Nks.pkl', 'r')
Nks = pickle.load(f)
f.close()

Hso = SpinOrbit.SpinOrbitOperator()
Hso.name = 'Hso'
Hso.pickle()
Hks =  KohnSham.KohnShamOperator() 
Hks.pickle()
Hopping = Hks + Hso
Hopping.name = 'Hopping'
Hopping.eigenvalues(Nks)
Hopping.pickle()
