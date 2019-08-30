import tensorflow as tf
import numpy as np
import utils

class StochasticWeightAveraging(tf.keras.callbacks.Callback):
    # def __init__(self):
    #     super().__init__()

    def on_train_begin(self, logs=None):
        self.w_swa = self.model.get_weights()

    def on_epoch_end(self, epoch, logs=None):
        print(f"Epoch {epoch} Updating average weights...")
        for i, _ in enumerate(self.model.layers):
            self.w_swa[i] = ((self.w_swa[i]*(epoch + 1) + self.model.get_weights()[i]) / (epoch + 2))

    def on_train_end(self, logs=None):
        self.model.set_weights(self.w_swa)
        if self.params["verbose"] == 1:
            print("Setting SWA weights.")
            print("WARNING: Make sure to run forward pass in train mode to recompute batch norm statistics")


def fit_with_swa(model, dataset, steps_per_epoch, save_path=None, lr=0.001, 
                 epochs=5, label_smoothing=0.05,):

    """
    Given a pretrained model.  Perform stochastic weight averaging with 
    training tf.data.Dataset.  Upon completion, replace model weights, then
    performs one last forward pass on the dataset in train mode to recompute
    the batch normalization statistics with the new weights.
    """

    model.compile(
        loss=tf.keras.losses.CategoricalCrossentropy(
            label_smoothing=label_smoothing), 
        optimizer=tf.keras.optimizers.SGD(learning_rate=lr),
        metrics=['CategoricalAccuracy']
    )

    # replace weights with the running average weights over course of training
    model.fit(
        dataset,
        steps_per_epoch=steps_per_epoch,
        epochs=epochs,
        callbacks=[StochasticWeightAveraging()],
        verbose=1
    )

    # forward pass to recompute batch norm statistics
    for x, _ in dataset.take(steps_per_epoch):
        model(x, training=True)

    if save_path:
        model.save(save_path)
    return model


if __name__ == '__main__':
    # TODO set up arg parsing
    BATCH_SIZE = 32
    
    dataset, n = utils.get_dataset(
        batch_size=32,
        data_dir='data/high_cov',
        training='train',
        augmentation=True,
    )

    steps_per_epoch = np.ceil(n/BATCH_SIZE)

    model = tf.keras.models.load_model('./saved_models/CNN_7_16.h5')

    fit_with_swa(
        model, dataset, steps_per_epoch,
        epochs=10,
        save_path='saved_models/CNN_7_16_SWA.h5',
        lr=0.1)





