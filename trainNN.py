'''
File to be run to train a neural network, assuming that model simulations are formatted like an image
	x- and y-dimensions represent the spatial size of the simulation
	color channels are cell types/properties

the argument taken in, n_nn, is which neural network in the ensemble is being trained

data is loaded from a numpy array of shape:
	[number_of_simulations, x_dimension, y_dimension, number_of_color_channels]


NEURAL NETWORK PARAMETERS
- cnnFilters: list of the number of convolutional filters at each layer. each element in the list is a 
              layer, with the integer value being the number of filters
- kernelSize: size of the convolutional kernel
- fullyConnected: list of the number of neurons in the fully-connected neural network that follows the
				  convolutional layers. each element is a hidden layer, with the integer value being the
				  number of neurons in that layer. the last layer is the projection layer
- activation: the activation function for the convolutional and hidden layers. note: the activation for
			  the last fullyConnected layer is linear (defined in genNN.py)
- dropout: neuron dropout during training
- learning_rate: training learning rate
- epochs: maximum number of training epochs
- patience: number of epochs without improvement before early stopping
- batchSize: number of simulations per training batch
- input_shape: size of a simulation (x-dimension, y-dimension, number of color channels)
- temperature: temperature parameter for NT-Xent loss
- batchNorm: True/False for turning batch normalization on
'''

from genNN import createTrainedModel
import os
import gc
import numpy as np
import sys

n_nn = int(sys.argv[1])

distortion = 0.01

data = np.load('trainingData/simulations.npy', allow_pickle=True)
input_shape = data[0].shape 
x_shape = [len(data)]
for i in input_shape:
	x_shape.append(i)
x = np.zeros(x_shape)
for i in range(len(data)):
	x[i] = data[i]

params = {'cnnFilters': [16,32],
		  'kernelSize': 3,
		  'fullyConnected': [128, 128, 128, 2],
		  'activation': 'selu',
		  'dropout': 0.0,
		  'learning_rate': 1e-4,
		  # 'rampup_epochs': 10,
		  # 'sustain_epochs': 0,
		  # 'lr_decay_exp': 0.95,
		  'epochs': 1000,
		  'patience': 50,
		  'batchSize': 5000,
		  'input_shape': input_shape,
		  'temperature': 1e-1,
		  'batchNorm': True}


os.system('mkdir -p sn_models/projector')
os.system('mkdir -p sn_models/loss')
projector, loss = createTrainedModel(params, x)
projector.save('sn_models/projector/model_'+str(n_nn))
np.savetxt('sn_models/loss/loss_'+str(n_nn)+'.csv', loss, delimiter=',')
np.savetxt('sn_models/acc/acc_'+str(n_nn)+'.csv', acc, delimiter=',')


