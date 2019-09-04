import functools
import tensorflow as tf


class Conv2DBlock:
    """
    Composition of 2D convolutions
    * Why isn't this a keras Layer?  I found that it just makes 
      whole model saving + loading a real pain so I just wrote this
      as if it was a keras Layer but did not inherit from Layer
    """
    def __init__(self, n_channels=32, n_layers=1, kernel_regularizer=None,
                 kernel_size=(3, 3), dilation_rate=(1, 1)):
        # layer properties
        self.n_channels = n_channels
        self.n_layers = n_layers
        self.kernel_regularizer = kernel_regularizer
        self.kernel_size = kernel_size
        self.dilation_rate = dilation_rate

        # sub layers
        self.conv_layers = [
            tf.keras.layers.Conv2D(
                filters=self.n_channels, kernel_size=self.kernel_size,
                dilation_rate=self.dilation_rate,
                kernel_regularizer=self.kernel_regularizer, padding='same')
            for i in range(self.n_layers)]
        self.bnorm_layers = [
            tf.keras.layers.BatchNormalization(renorm=True)
            for i in range(self.n_layers)]
        self.leaky_relu_layers = [
            tf.keras.layers.LeakyReLU()
            for i in range(self.n_layers)]

    def __call__(self, x):
        for conv, bnorm, leaky_relu in zip(
            self.conv_layers, self.bnorm_layers, self.leaky_relu_layers):
            x = conv(x)
            x = bnorm(x)
            x = leaky_relu(x)
        return x


class ResidualBlock(Conv2DBlock):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add = tf.keras.layers.Add()
        self.leaky_relu_out = tf.keras.layers.LeakyReLU()

    def __call__(self, x):
        temp = x
        x = super().__call__(x)
        x = self.add([temp, x])
        return self.leaky_relu_out(x)


def CNN(l2_strength=0.001):
    """
    Construct and return an (uncompiled) conv2d model out of Conv2DBlocks.
    """
    inp = tf.keras.Input(shape=(None, None, 3))
    x = inp
    x = tf.keras.layers.Conv2D(
        filters=32, kernel_size=(7, 7), strides=(1, 1),
        dilation_rate=(1, 1), padding='valid',
        kernel_regularizer=tf.keras.regularizers.l2(l2_strength))(x)
    x = tf.keras.layers.BatchNormalization(renorm=True)(x)
    x = tf.keras.layers.LeakyReLU()(x)
    x = tf.keras.layers.MaxPool2D(pool_size=(2, 2), strides=(2, 2))(x)

    for i in range(3):
        x = Conv2DBlock(
            n_channels=2**(i+5), n_layers=1, kernel_size=(1, 1),
            kernel_regularizer=tf.keras.regularizers.l2(l2_strength))(x)
        x = ResidualBlock(
            n_channels=2**(i+5), n_layers=2, kernel_size=(3, 3),
            kernel_regularizer=tf.keras.regularizers.l2(l2_strength))(x)
        x = ResidualBlock(
            n_channels=2**(i+5), n_layers=2, kernel_size=(3, 3),
            kernel_regularizer=tf.keras.regularizers.l2(l2_strength))(x)
        x = tf.keras.layers.MaxPool2D(pool_size=(2, 2), strides=(2, 2))(x)

    x = tf.keras.layers.GlobalAveragePooling2D()(x)
    for i in range(1):
        x = tf.keras.layers.Dense(1024//(i+1))(x)
        x = tf.keras.layers.LeakyReLU()(x)
        x = tf.keras.layers.Dropout(0.5)(x)

    x = tf.keras.layers.Dense(3)(x)
    out = tf.keras.layers.Softmax()(x)
    return tf.keras.Model(inputs=inp, outputs=out)


def Baseline(input_shape):
    """
    Simple model that just takes the image and flattens it for a feedforward
    neural network.
    """
    inp = tf.keras.Input(shape=input_shape)
    x = tf.keras.layers.Flatten()(inp)
    x = tf.keras.layers.Dense(1024, activation='relu')(x)
    out = tf.keras.layers.Dense(3, activation='softmax')(x)
    return tf.keras.Model(inputs=inp, outputs=out)
