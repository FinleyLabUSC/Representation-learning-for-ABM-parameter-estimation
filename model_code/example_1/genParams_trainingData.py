import numpy as np
import sys
import os

fld = sys.argv[1]
paramSet = int(sys.argv[2])
saveFld = sys.argv[3]+'/params'

# force params
m = 50
k = 12
ol = 0.2
d = 10

mcParams = np.loadtxt(fld+'/mcParams.csv', delimiter=',')

kill = mcParams[paramSet,0]
infil = mcParams[paramSet,1]
pdl1m = mcParams[paramSet,2]
pdl1g = mcParams[paramSet,3]
influ = mcParams[paramSet,4]
recRate = mcParams[paramSet,5]
recRatio = mcParams[paramSet,6]

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
recParams[0] = recRate #0.01 # cd8RecRate
recParams[1] = recRatio #0.3  # cd8Ratio
recParams[2] = 200  # recDist (recruit a uniform distribution recDist away from the tumor edge)

envParams = np.zeros((2, 1))
envParams[0] = 40  # simulation duration (days)
envParams[1] = 0 # 3d? 0 - no, 1 - yes

os.system('mkdir -p ' + saveFld)

np.savetxt(saveFld + '/cellParams.csv', cellParams, delimiter=',')
np.savetxt(saveFld + '/recParams.csv', recParams, delimiter=',')
np.savetxt(saveFld + '/envParams.csv', envParams, delimiter=',')