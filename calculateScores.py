'''
calculates the distance between model predictions and the image the simulations are being fit to

input
-----
- nModels: number of neural networks in the ensemble


loads simulations from format_simulations.py
loads the formated image (here, modelPredictions/base.npy)
projects the image and the simulations into low-dimensional space using the neural networks
takes the euclidean distance between each simulation and the image
    averages the distance calculated by each neural network in the ensemble
saves these distances to fittingInfo/scores.csv
adds these scores to fittingInfo/allScores.csv (which tracks the scores for all fitting steps)
'''

import numpy as np
from sn_snModel import loadModel
import sys
import os

nModels = int(sys.argv[1])

modelSimulations = np.load('modelPredictions/simulations.npy', allow_pickle=True)
shape = [len(modelSimulations)]
for i in modelSimulations.shape:
    shape.append(i)
xTest = np.zeros(shape)
for i in range(len(modelSimulations)):
    xTest[i] = modelSimulations[i]

base = np.load('modelPredictions/base.npy', allow_pickle=True)
shape = [1]
for i in base.shape:
    shape.append(i)
base = np.reshape(base, shape)


results = []
for m in range(nModels):
    print(m)
    model = loadModel('sn_models/projector/model_'+str(m))
    pointsTest = model.predict(xTest)
    pointBase = model.predict(base)
    if len(pointsBase.shape) == 1:
        pointBase = np.reshape(pointBase, [1, pointBase.shape[0]])

    distances = []
    for i in range(pointsTest.shape[0]):
        dx2 = []
        for j in range(pointsTest.shape[1]):
            dx2.append((pointsTest[i,j] - pointBase[0,j])**2)
        distances.append(np.sqrt(np.sum(dx2)))
    results.append(distances)

results = np.vstack(results)
results = np.mean(results, axis=0)
np.savetxt('fittingInfo/scores.csv', results, delimiter=',')

try:
    allScores = np.loadtxt('fittingInfo/allScores.csv', delimiter=',')
    if len(allScores.shape) == 1:
        allScores = np.reshape(allScores, [-1,1])
    allScores = np.concatenate((allScores, scores), axis=1)
    np.savetxt('fittingInfo/allScores.csv', allScores, delimiter=',')
except:
    np.savetxt('fittingInfo/allScores.csv', scores, delimiter=',')