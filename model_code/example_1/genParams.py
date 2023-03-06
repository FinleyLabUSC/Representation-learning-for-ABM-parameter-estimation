import numpy as np
import sys
import os

print()
print('launched gen params')
print()

paramSet = sys.argv[1]
saveFld = sys.argv[2]

# force params
m = 50  # 70.0#np.random.choice(mu)
k = 12  # 6.0
ol = 0.2  # np.random.choice(maxOverlap)
d = 10  # 9.0


modelParams = np.loadtxt('fittingInfo/paramsModel.csv', delimiter=',')
kill = modelParams[paramSet, 0]
infil = modelParams[paramSet, 1]
pdl1m = modelParams[paramSet, 2]
pdl1g = modelParams[paramSet, 3]
influ = modelParams[paramSet, 4]

cellParams = np.zeros((12, 2))

# cancer params
cellParams[0, 0] = m  # mu
cellParams[1, 0] = k  # kc
cellParams[2, 0] = d  # damping
cellParams[3, 0] = ol  # overlap
cellParams[4, 0] = 1 / 35  # div probability (hours) 
cellParams[5, 0] = 1/(24*10) # death probability (hours) 
cellParams[6, 0] = pdl1m  # pdl1 when expressed
cellParams[7, 0] = pdl1g  # prob of gaining pdl1
cellParams[8, 0] = 20  # diameter (um)

# cd8 params
cellParams[0, 1] = m  # mu
cellParams[1, 1] = k  # kc
cellParams[2, 1] = d  # damping
cellParams[3, 1] = ol  # overlap
cellParams[4, 1] = 1/(24*3) # death probability 
cellParams[5, 1] = 240  # migration speed um/hr
cellParams[6, 1] = kill # killProb 
cellParams[7, 1] = influ  # influence distance
cellParams[8, 1] = infil  # max infiltration distance
cellParams[9, 1] = 0.3  # migration bias 
cellParams[10, 1] = 0.1  # decrease of influence when suppressed
cellParams[11, 1] = 10  # diameter (um) 

recParams = np.zeros((3, 1))
recParams[0] = 0.01 # cd8RecRate
recParams[1] = 0.3  # cd8Ratio
recParams[2] = 200  # recDist (recruit a uniform distribution recDist away from the tumor edge)

envParams = np.zeros((2, 1))
envParams[0] = 40  # simulation duration (days)
envParams[1] = 0 # 3d? 0 - no, 1 - yes

os.system('mkdir -p ' + saveFld)

np.savetxt(saveFld + '/cellParams.csv', cellParams, delimiter=',')
np.savetxt(saveFld + '/recParams.csv', recParams, delimiter=',')
np.savetxt(saveFld + '/envParams.csv', envParams, delimiter=',')

print('done generating params')