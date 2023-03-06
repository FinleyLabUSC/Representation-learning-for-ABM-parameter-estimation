import numpy as np
import sys
import os

fld = sys.argv[1]
nSims = int(sys.argv[2])

# force params
m = 50
k = 12 
ol = 0.2
d = 10

params = []
for _ in range(nSims):
    kill = np.random.uniform(1e-3,1e-1)
    infil = np.random.uniform(0.0,1.0)
    pdl1m = np.random.uniform(1e-3,1e-1)
    pdl1g = np.random.uniform(1e-7,1e-4)
    influ = np.random.uniform(30,150)
    recRate = np.random.uniform(1e-2,1e-1)
    recRatio = np.random.uniform(0.1,0.5)

    params.append([kill,
                   infil,
                   pdl1m,
                   pdl1g,
                   influ,
                   recRate,
                   recRatio])

params = np.vstack(params)



os.system('mkdir -p ' + fld)

np.savetxt(fld + '/mcParams.csv', params, delimiter=',')
