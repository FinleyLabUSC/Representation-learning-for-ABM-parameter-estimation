'''
code for generating and training a neural network for projecting model simulations formatted as images
to low-dimensional space

the function that is run is createTrainedModel(params, data)
is outputs a trained encoder and training loss

training method is taken from SimCLR: https://arxiv.org/abs/2002.05709
code based on https://keras.io/examples/vision/semisupervised_simclr/
'''

import tensorflow as tf 
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np


def prepare_dataset(data, batch_size):
    dataset = tf.data.Dataset.from_tensor_slices(data).shuffle(10*batch_size).batch(batch_size)
    return dataset


def getC2D(input_shape, conv, kern, fc, activation, dropout, batchNorm):
    model = [keras.Input(shape=input_shape)]
    for i in conv:
        model.append(layers.Conv2D(filters=i, kernel_size=kern, activation=activation))
        model.append(layers.MaxPooling2D(pool_size=2))
    model.append(layers.Flatten())
    for i in range(len(fc)-1):
        model.append(layers.Dense(fc[i], activation=activation))
        if batchNorm == True:
            model.append(layers.BatchNormalization())
        model.append(layers.Dropout(dropout))
    model.append(layers.Dense(fc[-1]))
    return keras.Sequential(model)


def get_encoder(params):
    return getC2D(params['input_shape'],
                  params['cnnFilters'],
                  params['kernelSize'],
                  params['fullyConnected'],
                  params['activation'],
                  params['dropout'],
                  params['batchNorm'])

class RandomRot90(layers.Layer):
    def __init__(self, **kwargs):
        super(RandomRot90, self).__init__(**kwargs)

    def call(self, inputs):
        rot = tf.random.uniform((),0,4,dtype=tf.int32)
        return tf.image.rot90(inputs, k=rot)


def get_augmenter(input_shape):    
    rotLayer = RandomRot90()
    return keras.Sequential([
            keras.Input(shape=input_shape),
            layers.RandomFlip('horizontal'),
            rotLayer,
            ])

class ContrastiveModel(keras.Model):
    def __init__(self, params):
        super().__init__()

        self.temperature = params['temperature']
        self.contrastive_augmenter = get_augmenter(params['input_shape'])
        self.encoder = get_encoder(params)

    def compile(self, contrastive_optimizer, **kwargs):
        super().compile(**kwargs)

        self.contrastive_optimizer = contrastive_optimizer
        self.contrastive_loss_tracker = keras.metrics.Mean(name="c_loss")
        self.contrastive_accuracy = keras.metrics.SparseCategoricalAccuracy(
            name="c_acc"
        )

    @property
    def metrics(self):
        return [
            self.contrastive_loss_tracker,
            self.contrastive_accuracy,
        ]

    def euclidean_similarity(self, p1, p2):
        # Get the dot product between all embeddings
        # shape (batch_size, batch_size)
        dot_product = tf.matmul(p1, tf.transpose(p2))

        # Get squared L2 norm for each embedding. We can just take the diagonal of `dot_product`.
        # This also provides more numerical stability (the diagonal of the result will be exactly 0).
        # shape (batch_size,)
        square_norm = tf.linalg.diag_part(dot_product)

        # Compute the pairwise distance matrix as we have:
        # ||a - b||^2 = ||a||^2  - 2 <a, b> + ||b||^2
        # shape (batch_size, batch_size)
        distances = tf.expand_dims(square_norm, 1) - 2.0 * dot_product + tf.expand_dims(square_norm, 0)

        # Because of computation errors, some distances might be negative so we put everything >= 0.0
        distances = tf.maximum(distances, 0.0)

        # Because the gradient of sqrt is infinite when distances == 0.0 (ex: on the diagonal)
        # we need to add a small epsilon where distances == 0.0
        mask = tf.cast(tf.equal(distances, 0.0), tf.float32)
        distances = distances + mask * 1e-16
        #distances = distances + 1e-16

        distances = tf.sqrt(distances)

        # Correct the epsilon added: set the distances on the mask to be exactly 0.0
        distances = distances * (1.0 - mask)

        # https://stats.stackexchange.com/questions/53068/euclidean-distance-score-and-similarity
        similarity = 1/(1 + distances)

        return similarity

    def contrastive_loss(self, projections_1, projections_2):
        # InfoNCE loss (information noise-contrastive estimation)
        # NT-Xent loss (normalized temperature-scaled cross entropy)

        # Cosine similarity: the dot product of the l2-normalized feature vectors
        #projections_1 = tf.math.l2_normalize(projections_1, axis=1)
        #projections_2 = tf.math.l2_normalize(projections_2, axis=1)
        #similarities = (
        #    tf.matmul(projections_1, projections_2, transpose_b=True) / self.temperature
        #)
        similarities = self.euclidean_similarity(projections_1, projections_2)

        # The similarity between the representations of two augmented views of the
        # same image should be higher than their similarity with other views
        batch_size = tf.shape(projections_1)[0]
        contrastive_labels = tf.range(batch_size)
        self.contrastive_accuracy.update_state(contrastive_labels, similarities)
        self.contrastive_accuracy.update_state(
            contrastive_labels, tf.transpose(similarities)
        )

        # The temperature-scaled similarities are used as logits for cross-entropy
        # a symmetrized version of the loss is used here
        loss_1_2 = keras.losses.sparse_categorical_crossentropy(
            contrastive_labels, similarities, from_logits=True
        )
        loss_2_1 = keras.losses.sparse_categorical_crossentropy(
            contrastive_labels, tf.transpose(similarities), from_logits=True
        )
        return (loss_1_2 + loss_2_1) / 2

    def train_step(self, data):
        x_train = data

        # Each image is augmented twice, differently
        augmented_train_1 = self.contrastive_augmenter(x_train, training=True)
        augmented_train_2 = self.contrastive_augmenter(x_train, training=True)
        with tf.GradientTape() as tape:
            projections_1 = self.encoder(augmented_train_1, training=True)
            projections_2 = self.encoder(augmented_train_2, training=True)
            contrastive_loss = self.contrastive_loss(projections_1, projections_2)
        gradients = tape.gradient(
            contrastive_loss,
            self.encoder.trainable_weights,
        )
        self.contrastive_optimizer.apply_gradients(
            zip(
                gradients,
                self.encoder.trainable_weights,
            )
        )
        self.contrastive_loss_tracker.update_state(contrastive_loss)

        return {m.name: m.result() for m in self.metrics}


def lrfn(epoch, max_lr, min_lr=0.0, rampup_epochs=10, start_lr=0.0, sustain_epochs=0, exp_decay=0.95):
    if epoch < rampup_epochs:
        return (max_lr - start_lr)/rampup_epochs * epoch + start_lr
    elif epoch < rampup_epochs + sustain_epochs:
        return max_lr
    else:
        return (max_lr - min_lr)*exp_decay**(epoch - rampup_epochs - sustain_epochs) + min_lr


def get_model(params):
    model = ContrastiveModel(params)
    return model

def createTrainedModel(params, data):
    model = get_model(params)
    stopEarly = tf.keras.callbacks.EarlyStopping(monitor='c_loss', patience=params['patience'])
    '''
    lr_decay = tf.keras.callbacks.LearningRateScheduler(lambda epoch: lrfn(epoch,
                                                                           params['learning_rate'],
                                                                           rampup_epochs=params['rampup_epochs'],
                                                                           sustain_epochs=params['sustain_epochs'],
                                                                           exp_decay=params['lr_decay_exp']), verbose=True)
    '''
    
    model.compile(contrastive_optimizer=keras.optimizers.Adam(params['learning_rate']))
    dataset = prepare_dataset(data, params['batchSize'])
    history = model.fit(dataset, 
                        epochs=params['epochs'],
                        callbacks=[stopEarly], 
                        verbose=True)
    return model.encoder, history.history['c_loss']

def loadModel(fld):
    return keras.models.load_model(fld)
