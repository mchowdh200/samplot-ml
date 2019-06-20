import tensorflow as tf
import utils


class Conv2DBlock(tf.keras.Model):
    """
    Block containing a conv1d layer followed by dropout, batchnorm, and pooling
    """
    def __init__(self,
                 filters=128, kernel_size=(6, 6),
                 strides=(1, 1), dilation_rate=(1, 1),
                 dropout_rate=0, pool_size=(2, 2),
                 data_format='channels_last',
                 padding='valid',
                 batch_norm=False):
        super().__init__()

        self.batch_norm = batch_norm

        self.conv2d = tf.keras.layers.Conv2D(
            filters=filters, kernel_size=kernel_size,
            strides=strides, dilation_rate=dilation_rate,
            data_format=data_format,
            padding=padding,
            kernel_initializer='glorot_uniform')

        self.leaky_relu = tf.keras.layers.Activation(
            tf.keras.layers.LeakyReLU())

        self.dropout = tf.keras.layers.Dropout(rate=dropout_rate)

        if batch_norm:
            self.normalization = tf.keras.layers.BatchNormalization()

        self.pool = tf.keras.layers.MaxPool2D(
            pool_size=pool_size,
            data_format=data_format)

    def call(self, input_tensor):
        x = self.conv2d(input_tensor)
        x = self.leaky_relu(x)
        x = self.dropout(x)
        if self.batch_norm:
            x = self.normalization(x)
        return self.pool(x)

    def compute_output_shape(self, input_shape):
        x = self.conv2d.compute_output_shape(input_shape)
        return self.pool.compute_output_shape(x)


class CNN(tf.keras.Model):
    def __init__(self):
        super().__init__()
        
        # TODO remove the first batch norm once I've gotten mean/stddev of the training set.
        self.batch_norm = tf.keras.layers.BatchNormalization()
        self.conv_layers = [Conv2DBlock(filters=32*(i+2), kernel_size=(3, 3),
                                        pool_size=(2, 2), batch_norm=True, 
                                        padding='same', dropout_rate=0.1) 
                            for i in range(5)]
        self.global_avg_pool = tf.keras.layers.GlobalAveragePooling2D()
        self.dense = tf.keras.layers.Dense(256)
        self.leaky_relu = tf.keras.layers.LeakyReLU()
        self.dense_out = tf.keras.layers.Dense(3, activation='softmax')

    def call(self, x):
        x = self.batch_norm(x)
        for layer in self.conv_layers:
            x = layer(x)
        x = self.global_avg_pool(x)
        x = self.dense(x)
        x = self.leaky_relu(x)
        return self.dense_out(x)





