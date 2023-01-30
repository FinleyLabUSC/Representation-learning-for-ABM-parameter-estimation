'''
example file for formatting model simulations
in this example, simulations were saved to ~/training_sims/simulation_IDXOFSIMULATION/set_IDXOFREPLICATE

params
------
- nSims: number of parameters in Monte Carlo set
- nReps: number of replicates for each parameter set
- gridSize: size of one grid space (the cell diameter)
- imSize: number of pixels in the final formatted image that go into the neural network
- bindingLayer: manually set to 0, which, in this example, was cancer cells

formats all of the simulations, then saves them in one file
'''

import numpy as np
from parseData import *
from data_discretizeFunctions import *
import os

nSims = 1000
nReps = 10
gridSize = 20
imSize = (50,50)
bindingLayer = 0

sims = []
for i in range(nSims):
    for j in range(nReps):
        layers = loadSingle('training_sims/simulation_'+str(i)+'/set_'+str(j))
        ds = DiscreteImg(gridSize, layers, bindingLayer)
        sg = ds.smallGrids(imSize)

        sims.append(sg)

print(len(sims))
os.system('mkdir -p trainingData')
np.save('trainingData/simulations.npy', sims)