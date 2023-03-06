import numpy as np
import sys
import os

paramSet = int(sys.argv[1])
saveFld = sys.argv[2] + '/params'
fld = sys.argv[3]

# force params
m = 50 
k = 12
ol = 0.2
d = 10

mcParams = np.loadtxt(fld+'/mcParams.csv', delimiter=',')

kill = mcParams[paramSet, 0]
infil = mcParams[paramSet, 1]
death = mcParams[paramSet, 2]
hypoxL = mcParams[paramSet, 3]
hypoxP = mcParams[paramSet, 4]

cellParams = np.zeros((10, 2))

# cancer params
cellParams[0, 0] = m  # mu
cellParams[1, 0] = k  # kc
cellParams[2, 0] = d  # damping
cellParams[3, 0] = ol  # overlap
cellParams[4, 0] = 1 / 35  # div probability (hours)
cellParams[5, 0] = death # death probability (hours)
cellParams[6, 0] = 20  # diameter (um)
cellParams[7, 0] = hypoxL # hypoxicL -> lambda for exponential decay for hypoxia
cellParams[8, 0] = hypoxP # max probability for death from hypoxia

# cd8 params
cellParams[0, 1] = m  # mu
cellParams[1, 1] = k  # kc
cellParams[2, 1] = d  # damping
cellParams[3, 1] = ol  # overlap
cellParams[4, 1] = 0 # death probability
cellParams[5, 1] = 240  # migration speed um/hr
cellParams[6, 1] = kill # killProb
cellParams[7, 1] = infil  # max infiltration distance
cellParams[8, 1] = 0.3  # migration bias 
cellParams[9, 1] = 10  # diameter (um)

recParams = np.zeros((3, 1))
recParams[0] = 0.0 # cd8RecRate
recParams[1] = 20  # cd8Ratio
recParams[2] = 500  # recDist (recruit a uniform distribution recDist away from the tumor edge)

envParams = np.zeros((3, 1))
envParams[0] = 3  # simulation duration (days)
envParams[1] = 0 # 3d? 0 - no, 1 - yes
envParams[2] = 75 # initial radii number of cells

os.system('mkdir -p ' + saveFld)

np.savetxt(saveFld + '/cellParams.csv', cellParams, delimiter=',')
np.savetxt(saveFld + '/recParams.csv', recParams, delimiter=',')
np.savetxt(saveFld + '/envParams.csv', envParams, delimiter=',')

print('done generating params')