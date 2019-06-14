import numpy as np
import tensorflow as tf
from utils import get_dataset

def conv2Dmodel(input_shape=None):

    model = tf.keras.Sequential([
        tf.keras.layers.Conv2D(#input_shape=input_shape,
                               filters=64,
                               kernel_size=(3, 3),
                               activation='relu'),
        tf.keras.layers.MaxPool2D(pool_size=(2, 2)),
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(2, activation='sigmoid')
    ])


    model.compile(loss='binary_crossentropy', 
                  optimizer='Adam', 
                  metrics=['Accuracy'])

    return model

if __name__ == '__main__':
    BATCH_SIZE = 32
    ds, n_images = get_dataset(batch_size=BATCH_SIZE)
    model = conv2Dmodel()
    model.fit(ds,
              steps_per_epoch=np.ceil(n_images/BATCH_SIZE),
              epochs=10)
