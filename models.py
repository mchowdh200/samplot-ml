import tensorflow as tf


class Conv2DBlock:
    """
    Composition of 2D convolutions
    """
    def __init__(self, n_channels=32, n_layers=1, kernel_regularizer=None,
                 kernel_size=(3, 3), dilation_rate=(1, 1),
                 batch_norm=False, dropout_rate=0.0):
        self.n_channels = n_channels
        self.n_layers = n_layers
        self.kernel_regularizer = kernel_regularizer
        self.kernel_size = kernel_size
        self.dilation_rate = dilation_rate
        self.batch_norm = batch_norm
        self.dropout_rate=dropout_rate

        self.conv_layers = [
            tf.keras.layers.Conv2D(
                filters=self.n_channels, kernel_size=self.kernel_size,
                dilation_rate=self.dilation_rate,
                kernel_regularizer=self.kernel_regularizer, padding='same')
            for i in range(self.n_layers)]


    def __call__(self, x):

        for layer in self.conv_layers:
            x = layer(x)
            if self.batch_norm:
                x = tf.keras.layers.BatchNormalization()(x)
            if self.dropout_rate > 0:
                x = tf.keras.layers.SpatialDropout2D(rate=self.dropout_rate)(x)

            x = tf.keras.layers.LeakyReLU()(x)
        return x


class ResidualBlock(Conv2DBlock):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __call__(self, x):
        temp = x
        x = super().__call__(x)
        return tf.keras.layers.LeakyReLU()(temp + x)



def CNN(dropout_rate=0.0):
    """
    Construct and return an (uncompiled) conv2d model out of Conv2DBlocks.
    """
    inp = tf.keras.Input(shape=(None, None, 3))
    x = inp
    x = tf.math.reduce_max(tf.image.sobel_edges(x), axis=4)

    x = tf.keras.layers.Conv2D(
        filters=32, kernel_size=(7, 7), strides=(1, 1),
        dilation_rate=(2, 2), padding='valid')(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.LeakyReLU()(x)
    x = tf.keras.layers.MaxPool2D(pool_size=(2, 2))(x)

    for i in range(4):
        x = Conv2DBlock(
            n_channels=32*(i+1), n_layers=1, kernel_size=(1, 1),
            batch_norm=True, dropout_rate=dropout_rate)(x)
        x = ResidualBlock(
            n_channels=32*(i+1), n_layers=2, kernel_size=(3, 3),
            batch_norm=True, dropout_rate=dropout_rate)(x)
        x = ResidualBlock(
            n_channels=32*(i+1), n_layers=2, kernel_size=(3, 3),
            batch_norm=True, dropout_rate=dropout_rate)(x)
        x = ResidualBlock(
            n_channels=32*(i+1), n_layers=2, kernel_size=(3, 3),
            batch_norm=True, dropout_rate=dropout_rate)(x)
        x = tf.keras.layers.MaxPool2D(pool_size=(2, 2), strides=(2, 2))(x)

    x = tf.keras.layers.GlobalAveragePooling2D()(x)
    for i in range(1):
        x = tf.keras.layers.Dense(1024//(i+1))(x)
        x = tf.keras.layers.LeakyReLU()(x)
        # x = tf.keras.layers.Dropout(0.5)(x)

    out = tf.keras.layers.Dense(3, activation='softmax')(x)
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
    

