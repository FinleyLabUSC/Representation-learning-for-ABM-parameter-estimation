'''
example file for formatting a tumor image

this is run after the extraction of cell coordinates from the image
saves to 'modelPredictions/base.npy' to be used later by calculate_scores.py
'''

import numpy as np
from parseData import *
from data_discretizeFunctions import *
import os

gridSize = 20
imSize = (50,50)

layers = loadSingle('PATH/TO/TUMOR/CELL/COORDINATES')
ds = DiscreteImg(gridSize, layers, 0)
sg = ds.smallGrids(imSize)

os.system('mkdir -p modelPredictions')
np.save('modelPredictions/base.npy', sg)