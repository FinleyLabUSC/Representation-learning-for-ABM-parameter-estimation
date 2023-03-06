'''
This file converts lists of different cells into an image
This is done by discretizing the continuous cell coordinates to a grid, then shrinking the grid to a set size
This transforms the explicit cell locations into cell densities and facilitates the comparison between model simulations
and actual tumor images. While some spatial detail is lost in the process, it is unlikely that the model captures
enough biological detail to fully replicate in vivo intricacies.

If the model is simulating on a 1-to-1 spatial scale with the image data, than resizing is not necessary
and createGrids() can be run in place of smallGrids(). However, care must be taken to ensure that
these grids end up the same size as each other.

inputs to class
---------------
- gridSize: the physical length of one side of a grid space (usually the size of one cell diameter)
- layers: list containing arrays of formated cell properties
          each element in the list is an array, and the element index corresponds to the color channel
            in the outputted image
          this list is produced by parseData.py
- bindingLayer: which layer to use to set the edges of the image
                for example, if the bindingLayer corresponds to cancer cells, the image will be bound
                to the edge of the tumor. If a different bindingLayer is chosen, it may bind the image
                to cells outside of the tumor. If bindingLayer == -1, image is bound to the outermost
                cells, regardless of type

smallGrids(imSize)
------------------
- this function discretizes all cells to grids to create an image, then shrinks it since all of the
  images need to be the same size to be inputs to the neural network
- input imSize: the number of x,y pixels of the shrunken image
- output: an (imSize[0], imSize[1], len(layers)) array that can be treated as an image
'''

import numpy as np
import cv2

class DiscreteImg:
    def __init__(self, gridSize, layers, bindingLayer):
        '''
        layers: list of cell types/properties. each list element is an array w/ columns [x, y, val]
        bindingLayer: layer that provides the edges of the environment. if -1, no binding layer
        '''
        self._layers = layers
        self._bl = bindingLayer
        self._scaledLayers = self.scaleLayers(gridSize)

    def scaleLayers(self, size):
        '''
        scales cell position to a discrete grid-size
        binds the grid to self._bl
        '''
        minX = 1e8
        minY = 1e8
        maxX = -1e8
        maxY = -1e8

        if self._bl == -1:
            for i in range(len(self._layers)):
                minX = np.min((minX, np.min(self._layers[i][:,0])))
                minY = np.min((minY, np.min(self._layers[i][:,1])))

                maxX = np.max((maxX, np.max(self._layers[i][:,0])))
                maxY = np.max((maxY, np.max(self._layers[i][:,1])))
        else:
            minX = np.min((minX, np.min(self._layers[self._bl][:,0])))
            minY = np.min((minY, np.min(self._layers[self._bl][:,1])))

            maxX = np.max((maxX, np.max(self._layers[self._bl][:,0])))
            maxY = np.max((maxY, np.max(self._layers[self._bl][:,1])))

        scaledLayers = []
        for l in self._layers:
            if l.shape[0] == 0:
                scaledLayers.append(l)
                continue

            scaled = []
            for i in range(l.shape[0]):
                if l[i,0] < minX or l[i,0] > maxX or l[i,1] < minY or l[i,1] > maxY:
                    continue
                x = (l[i,0] - minX)/size
                y = (l[i,1] - minY)/size
                v = l[i,2]
                if np.isnan(x):
                    exit('x')
                
                scaled.append([x,y,v])
            scaledLayers.append(np.vstack(scaled))

        for l in scaledLayers:
            if l.shape[0] == 0:
                continue
            if np.isnan(np.max(l)):
                exit('asdfadsf')

        return scaledLayers

    def createGrids(self):
        '''
        creates a grid from scaledLayers
        treats each grid-site as a pixel in an image
        each layer is a channel, akin to RGB
        '''
        sizeX = 0
        sizeY = 0
        for l in self._scaledLayers:
            if l.shape[0] == 0:
                continue
            sizeX = int(np.max((sizeX, np.max(l[:,0]))))
            sizeY = int(np.max((sizeY, np.max(l[:,1]))))

        sizeX += 1
        sizeY += 1

        grid = np.zeros((sizeX, sizeY, len(self._scaledLayers)))
        for k in range(len(self._scaledLayers)):
            for idx in range(self._scaledLayers[k].shape[0]):
                i = int(self._scaledLayers[k][idx, 0])
                j = int(self._scaledLayers[k][idx, 1])
                v = self._scaledLayers[k][idx, 2]
                grid[i,j,k] = v

        return grid

    def smallGrids(self, imSize):
        '''
        creates cell grid then shrinks it to imSize
        '''
        grid = self.createGrids()
        if imSize[0] == 0:
            return grid

        gridSmall = cv2.resize(grid, imSize, interpolation=cv2.INTER_AREA)
        if len(self._scaledLayers) == 1:
            gridSmall = gridSmall.reshape(gridSmall.shape[0], gridSmall.shape[1], 1)

        for i in range(gridSmall.shape[2]):
            if np.max(gridSmall[:,:,i]) == 0:
                continue
            gridSmall[:,:,i] = gridSmall[:,:,i]/np.max(gridSmall[:,:,i])

        return gridSmall
