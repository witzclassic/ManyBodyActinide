import cPickle as pickle

f = open('Hks.pkl', 'r')
Hks = pickle.load(f)
f.close()
f = open('Hopping.pkl', 'r')
Hop = pickle.load(f)
f.close()
f = open('nParticles.pkl', 'r')
nParticles = pickle.load(f)
f.close()
Hks.eigenvalues(nParticles)
Hop.eigenvalues(nParticles)

