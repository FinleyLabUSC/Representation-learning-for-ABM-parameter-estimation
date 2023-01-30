'''
loads and formats model simulations
this example is for a T cell - tumor model with PD-L1

outputs a list of arrays, where the list index corresponds with the color channel in the formated image
each array is two dimensional: (number of cells of that type, 3),
                               where the second dimension is [x coordinate, y coordinate, value]
                               for cell types, value is 1 to indicate that a cell is present at that location
                               for continuous properties, the property is normalized to 0-1
'''

import numpy as np

def loadSingle(fld):
    '''
    loads cells and organizes them into a list of arrays, one array per cell type/property
    cell array - [x, y, value]
    value - 1 for cell is present, continuous value for pdl1
    normalizes pdl1 to 0-1

    arrays:
        cancer
        cd8
        suppressed cd8
        pdl1
    '''

    cancer = np.loadtxt(fld+'/cancerCells.csv', delimiter=',')
    if cancer.size == 0:
        return []
    cd8 = np.loadtxt(fld+'/cd8Cells.csv', delimiter=',')

    cancerArray = []
    cd8Array = []
    cd8sArray = []
    pdl1Array = []

    for i in range(cancer.shape[0]):
        cancerArray.append([cancer[i,0],
                      cancer[i,1],
                      1])
        pdl1Array.append([cancer[i,0],
                      cancer[i,1],
                      cancer[i,4]/np.max(cancer[:,4])])

    for i in range(cd8.shape[0]):
        if cd8[i,4] == 0:
            cd8Array.append([cd8[i,0],
                          cd8[i,1],
                          1])
        elif cd8[i,4] == 1:
            cd8sArray.append([cd8[i,0],
                          cd8[i,1],
                          1])
        else:
            exit("incorrect cd8 state")

    cancerArray = np.vstack(cancerArray)
    cd8Array = np.vstack(cd8Array)
    cd8sArray = np.vstack(cd8sArray)
    pdl1Array = np.vstack(pdl1Array)

    return [cancerArray, cd8Array, cd8sArray, pdl1Array]


def loadSimulations(fld, nSims):
    simulations = []

    for i in range(nSims):
        sim = loadSingle(fld+'/set_'+str(i))
        if len(sim) == 0:
            return []
        simulations.append(sim)

    return simulations
