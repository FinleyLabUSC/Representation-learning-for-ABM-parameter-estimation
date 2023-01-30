'''
same as formatData.py, but formats the simulations used in the parameter estimation steps
'''


import numpy as np
from parseData import *
from data_discretizeFunctions import *
import os
import sys

nSims = int(sys.argv[1])
gridSize = int(sys.argv[2])
imSize = (int(sys.argv[3],int(sys.argv[3]))
bindingLayer = 0

sims = []
for i in range(nSims):
    layers = loadSingle('modelPredictions/simulation_'+str(i)+'/set_0')
    ds = DiscreteImg(gridSize, layers, bindingLayer)
    sg = ds.smallGrids(imSize)

    sims.append(sg)

print(len(sims))
np.save('modelPredictions/simulations.npy', sims)