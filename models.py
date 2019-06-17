import numpy as np
import tensorflow as tf
from utils import get_dataset

tf.compat.v1.disable_eager_execution()

def conv2Dmodel():
    model = tf.keras.Sequential()
    for i in range(3):
        model.add(tf.keras.layers.Conv2D(filters=32*(i+1),
                                         kernel_size=(6, 6),
                                         activation=tf.keras.layers.LeakyReLU()))
        # model.add(tf.keras.layers.Dropout(rate=0.2))
        # model.add(tf.keras.layers.BatchNormalization())
        model.add(tf.keras.layers.MaxPool2D(pool_size=(2, 2)))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(3, activation='softmax'))
    
    model.compile(loss='sparse_categorical_crossentropy', 
                  optimizer=tf.optimizers.Adam(lr=2e-3), 
                  metrics=['SparseCategoricalAccuracy'])
    return model

if __name__ == '__main__':
    BATCH_SIZE = 64
    training_set, n_train = get_dataset(batch_size=BATCH_SIZE, training='train')
    test_set, n_test = get_dataset(batch_size=BATCH_SIZE, training='test')


    callbacks = [tf.keras.callbacks.ReduceLROnPlateau(patience=3),
                 tf.keras.callbacks.EarlyStopping(patience=5)]

    model = conv2Dmodel()
    model.fit(training_set,
              steps_per_epoch=np.ceil(n_train/BATCH_SIZE),
              validation_data=test_set,
              validation_steps=np.ceil(n_test/BATCH_SIZE),
              epochs=10,
              callbacks=callbacks)
