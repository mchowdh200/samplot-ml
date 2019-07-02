import tensorflow as tf


class Conv2DBlock:
    """
    Block containing a conv2d layer followed by dropout, batchnorm (optional), and pooling
    """
    def __init__(self, filters=128, kernel_size=(6, 6),
                 strides=(1, 1), dilation_rate=(1, 1),
                 dropout_rate=0.0, pool_size=(2, 2),
                 kernel_regularizer=None,
                 data_format='channels_last',
                 padding='valid',
                 batch_norm=False):
        self.filters = filters
        self.kernel_size = kernel_size
        self.strides = strides
        self.dilation_rate = dilation_rate
        self.dropout_rate = dropout_rate
        self.pool_size = pool_size
        self.kernel_regularizer = kernel_regularizer
        self.data_format = data_format
        self.padding = padding
        self.batch_norm = batch_norm

    def __call__(self, x):
        x = tf.keras.layers.Conv2D(
            filters=self.filters, kernel_size=self.kernel_size,
            strides=self.strides, dilation_rate=self.dilation_rate,
            data_format=self.data_format, padding=self.padding,
            kernel_regularizer=self.kernel_regularizer,
            kernel_initializer='glorot_uniform')(x)
        x = tf.keras.layers.LeakyReLU()(x)
        x = tf.keras.layers.Dropout(rate=self.dropout_rate)(x)
        if self.batch_norm:
            x = tf.keras.layers.BatchNormalization()(x)
        return tf.keras.layers.MaxPool2D(pool_size=self.pool_size)(x)


def CNN():
    """
    Construct and return an (uncompiled) conv2d model out of Conv2DBlocks.
    """
    inp = tf.keras.Input(shape=(None, None, 3))
    x = inp

    for i in range(5):
        x = Conv2DBlock(filters=32*(i+1), kernel_size=(3, 3),
                        pool_size=(2, 2), batch_norm=True,
                        padding='same', dropout_rate=0.2)(x)
        # TODO try 1x1 convolutions
    x = tf.keras.layers.GlobalAveragePooling2D()(x)
    for i in range(2):
        x = tf.keras.layers.Dense(1024)(x)
        x = tf.keras.layers.LeakyReLU()(x)
        x = tf.keras.layers.Dropout(0.5)(x)
    out = tf.keras.layers.Dense(3, activation='softmax')(x)
    return tf.keras.Model(inputs=inp, outputs=out)


