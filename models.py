import numpy as np
import tensorflow as tf
from utils import get_dataset

def conv2Dmodel():
    model = tf.keras.Sequential()
    for i in range(2):
        model.add(tf.keras.layers.Conv2D(filters=32*(i+1),
                                         kernel_size=(6, 6),
                                         strides=(1, 1)))
        # model.add(tf.keras.layers.Dropout(rate=0.2))
        # model.add(tf.keras.layers.BatchNormalization(renorm=True))
        model.add(tf.keras.layers.LeakyReLU())
        model.add(tf.keras.layers.MaxPool2D(pool_size=(2, 2)))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(3, activation='softmax'))
    
    model.compile(loss='sparse_categorical_crossentropy', 
                  optimizer=tf.optimizers.Adam(lr=5e-4, amsgrad=True), 
                  metrics=['SparseCategoricalAccuracy'])
    return model
