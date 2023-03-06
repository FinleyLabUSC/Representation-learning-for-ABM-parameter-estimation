import numpy as np
import sys
import os

fld = sys.argv[1]
nSims = int(sys.argv[2])

# force params
m = 50  # 70.0#np.random.choice(mu)
k = 12  # 6.0
ol = 0.2  # np.random.choice(maxOverlap)
d = 10  # 9.0

params = []
for _ in range(nSims):
    kill = np.random.uniform(1e-3,1e-1)
    infil = np.random.uniform(0.0,1.0)
    death = 1/np.random.uniform(2*24, 10*24)
    hypoxL = np.random.uniform(1e2, 1e3)
    hypoxP = np.random.uniform(1e-4,1e-1)

    params.append([kill,
                   infil,
                   death,
                   hypoxL,
                   hypoxP])

params = np.vstack(params)

os.system('mkdir -p ' + fld)

np.savetxt(fld + '/mcParams.csv', params, delimiter=',')
