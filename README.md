# Representation-learning-for-ABM-parameter-estimation

Code from https://www.biorxiv.org/content/10.1101/2023.01.12.523847v1

Files that are run:

   - formatTumorImage.py - formats the tumor image that the model is being fit to. Should be run before fitting
   
   - formatData.py - formats model simulations into simplified images for training the neural network
   
   - format_simulations.py - formats model simulations to be inputted into the trained neural network
   
   - trainNN.py - trains a neural network to project 
   
   - calculateScores.py - uses the simulations from format_simulations.py to calculate the distance of each simulation from the image
   
 Files that are called by other files:
 
   - parseData.py - turns lists of cells into a list of arrays, with each array containing cell locations and the value for a cell property
   
   - data_discretizeFunctions.py - turns the lists generated by parseData.py into a simplified image
   
   - genNN.py - generates and trains a neural network for representation learning
   
Note: code for the genetic algorithm is not provided as it was written specifically to run on our university's computing cluster and interface between the neural network code (written in Python) and the test models (written in C++). It does not run on a local desktop without modification.
