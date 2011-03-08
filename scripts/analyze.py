import numpy as np
import AngMo as am
import ManyBodyState as mbs
import ManyBodyHam as mbh

for iter in range(8,11):
    
    fname = 'LLS.Veal.F0.'+str(iter*0.5)+'.mtx'
    print
    print '<====== ' + fname + '========>'
    print
    vecFile = open(fname, 'r')

    [vectorSizeTxt, gsEnergyTxt] = vecFile.readline().split()  
    print vectorSizeTxt, gsEnergyTxt
    vectorSize = int(vectorSizeTxt)
    vector = np.zeros(vectorSize, dtype=complex)
    index = 0
    
    for line in vecFile:
        [realTxt, imagTxt] = line.split()
        #print realTxt, imagTxt
        vector[index] = complex(float(realTxt), float(imagTxt))
        index = index + 1

    vecFile.close()

    mbs.analyzeOccupations(vector)
    readFlag = True

    s2Ham = mbh.ManyBodyHamMTX('SSquaredOp', readFlag)
    exp = s2Ham.expValue(vector)
    print 'S2 exp:', exp
    l2Ham = mbh.ManyBodyHamMTX('LSquaredOp', readFlag)
    exp = l2Ham.expValue(vector)
    print 'L2 exp:', exp
    j2Ham = mbh.ManyBodyHamMTX('JSquaredOp', readFlag)
    exp = j2Ham.expValue(vector)
    print 'J2 exp:', exp
