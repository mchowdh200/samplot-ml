import functools
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
                # TODO compare with batch renormalization
                x = tf.keras.layers.BatchNormalization(renorm=True)(x)
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
        x = tf.keras.layers.LeakyReLU()(x)
        # this seems to cause problems when
        # trying to save the model             
        # return tf.keras.layers.LeakyReLU()(temp + x) 
        x = tf.keras.layers.Add()([temp, x])
        return tf.keras.layers.LeakyReLU()(x)


class GroupedConvolution:
    """
    This "layer" will take an input tensor, slice it up into equal pieces,
    apply a conv2d layer to each, and concatenate the results.
    """
    def __init__(self, channels, cardinality, kernel_regularizer=None,
                 kernel_size=(3, 3), dilation_rate=(1, 1), strides=(1, 1)):
        assert channels % cardinality == 0, "cardinality does not evenly divide channels."
        self.channels = channels
        self.cardinality = cardinality
        self.kernel_regularizer = kernel_regularizer
        self.kernel_size = kernel_size
        self.dilation_rate = dilation_rate
        self.strides=strides

        # channels per group
        n = self.channels // self.cardinality

        # slices up the input by the channel dimension,
        self.split = [
            tf.keras.layers.Lambda(
                functools.partial(
                    lambda x, i: x[:, :, :, i*n : i*n + n], i=i))
            for i in range(self.cardinality)]
        # self.split = []
        # for i in range(self.cardinality):
        #     self.split.append(
        #         tf.keras.layers.Lambda(
        #             lambda x: x[:, :, :, i*n : i*n + n]
        #         )
        #     )

        # then applies conv2d to each slice independently
        self.group_conv = [
            tf.keras.layers.Conv2D(
                filters=n,
                kernel_size=self.kernel_size,
                strides=self.strides,
                padding='same')
            for _ in range(self.cardinality)]

        self.batch_norm = tf.keras.layers.BatchNormalization()
        self.leaky_relu = tf.keras.layers.LeakyReLU()

    def __call__(self, x):
        groups = [f(x) for f in self.split]
        groups = [conv(group) for (conv, group) in zip(self.group_conv, groups)]
        y = tf.keras.layers.concatenate(groups)
        return self.leaky_relu(self.batch_norm(y))
        


class ResNeXtBlock:
    def __init__(self, channels_in, channels_out, cardinality=32, 
                 downsample=True, project_shortcut=False):
        self.channels_in = channels_in
        self.channels_out = channels_out
        self.cardinality = cardinality
        self.downsample = downsample
        self.project_shortcut = project_shortcut

        self.proj_down = tf.keras.layers.Conv2D(
            self.channels_in,
            kernel_size=(1, 1),
            strides=(1, 1),
            padding='same')
        self.batch_norm1 = tf.keras.layers.BatchNormalization()
        self.leaky_relu1 = tf.keras.layers.LeakyReLU()

        self.grouped_conv = GroupedConvolution(
            channels=self.channels_in,
            cardinality=self.cardinality,
            kernel_size=(3, 3),
            strides=(3, 3) if self.downsample else (1, 1))

        # used to project the shortcut if the input does not have the same
        # number of channels as the output.
        self.shortcut_proj = tf.keras.layers.Conv2D(
            self.channels_out,
            kernel_size = (1, 1),
            strides=(3, 3) if self.downsample else (1, 1),
            padding='same')

        self.proj_up = tf.keras.layers.Conv2D(
            self.channels_out,
            kernel_size=(1, 1),
            strides=(1, 1),
            padding='same')
        self.batch_norm2 = tf.keras.layers.BatchNormalization()
        self.leaky_relu2 = tf.keras.layers.LeakyReLU()

        # applied after residual connection
        self.leaky_relu3 = tf.keras.layers.LeakyReLU()

    def __call__(self, x):
        shortcut = x

        x = self.proj_down(x)
        x = self.batch_norm1(x)
        x = self.leaky_relu1(x)

        # if x is not going to be the same dimension as the output,
        # let the shortcut start after the first projection layer
        if self.project_shortcut or self.downsample:
            shortcut = self.shortcut_proj(shortcut)

        x = self.grouped_conv(x)
        x = self.proj_up(x)
        x = self.batch_norm2(x)
        x = self.leaky_relu2(x)
        x = tf.keras.layers.add([x, shortcut])
        return self.leaky_relu3(x)


def ResNeXt(base_filters=32):

    filters = base_filters

    # input stages
    inp = tf.keras.Input(shape=(None, None, 3))
    x = tf.keras.layers.Conv2D(
        filters=filters,
        kernel_size=(7, 7),
        strides=(2, 2))(inp)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.LeakyReLU()(x)
    x = tf.keras.layers.MaxPool2D(pool_size=(3, 3), strides=(2, 2))(x)

    filters *= 2

    # residual blocks
    for i in range(3):
        for j in range(3):
            x = ResNeXtBlock(channels_in=filters,
                             channels_out=filters*2,
                             project_shortcut=(j == 0 and i == 0),
                             downsample=(j == 0),
                             cardinality=32)(x)
        filters *= 2

    # for i in range(3):
    #     x = ResNeXtBlock(channels_in=FILTERS,
    #                      channels_out=FILTERS*4,
    #                      project_shortcut=False,
    #                      downsample=i == 0,
    #                      cardinality=32)(x)

    # for i in range(3):
    #     x = ResNeXtBlock(channels_in=FILTERS,
    #                      channels_out=FILTERS*2,
    #                      project_shortcut=False,
    #                      downsample=i == 0,
    #                      cardinality=32)(x)

    # output stages
    x = tf.keras.layers.GlobalAveragePooling2D()(x)
    x = tf.keras.layers.Dense(3)(x)
    out = tf.keras.layers.Softmax()(x)
    return tf.keras.Model(inputs=inp, outputs=out)

            
def CNN(dropout_rate=0.0):
    """
    Construct and return an (uncompiled) conv2d model out of Conv2DBlocks.
    """
    inp = tf.keras.Input(shape=(None, None, 3))
    x = inp
    # x = tf.math.reduce_max(tf.image.sobel_edges(x), axis=4)

    x = tf.keras.layers.Conv2D(
        filters=32, kernel_size=(7, 7), strides=(1, 1),
        dilation_rate=(2, 2), padding='valid')(x)
    # TODO compare with batch renormalization
    x = tf.keras.layers.BatchNormalization(renorm=True)(x)
    x = tf.keras.layers.LeakyReLU()(x)
    x = tf.keras.layers.MaxPool2D(pool_size=(2, 2))(x)

    for i in range(3):
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
    

