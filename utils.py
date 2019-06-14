import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split


def get_filenames(df):
    """
    Get filenames of images which are composed of fields in the dataframe
    """
    return f'./data/crop/{df.SVTYPE}_{df.chrom}_{df.start}-{df.end}.crop.png'

def get_labels(df, threshold=0.5):
    """
    Apply labels to the dataframe. only consider DELs with a confidence 
    greater than the threshold.
    """
    if df.SVTYPE == 'DEL' and df.score > threshold:
        return 1
    return 0


def load_image(path, label):
    """
    Used to load image from provided filepath for use with a tensorflow dataset.
    Most images we have are 3 channel, but there are some that are 1/4 channels,
    so we just make all 3 channel then normalize to 0-1 range.  We pass through
    the image label so we can yield it in a tensorflow datset
    """
    image = tf.io.read_file(path)
    return tf.image.decode_png(image, channels=3)/255, label


def get_dataset(data_dir='./data', training=True, batch_size=32):
    """
    Create a tensorflow dataset that yields image/label pairs to a model
    """

    if training:
        train_test = 'train'
    else:
        train_test = 'test'


    df = pd.read_csv(f"{data_dir}/{train_test}.bed", sep='\t',
                     names=['chrom', 'start', 'end', 'SVTYPE', 'score'])

    n_images = len(df)

    df['filename'] = df.apply(get_filenames, axis=1)
    df['label'] = df.apply(get_labels, axis=1)

    # construct the dataset
    AUTOTUNE = tf.data.experimental.AUTOTUNE

    ds = tf.data.Dataset.from_tensor_slices(
        (df.filename.values, tf.cast(df.label.values, tf.int64)))

    ds = ds.map(load_image)
    ds = ds.cache()
    ds = ds.apply(tf.data.experimental.shuffle_and_repeat(buffer_size=n_images))
    ds = ds.batch(batch_size)
    ds = ds.prefetch(buffer_size=AUTOTUNE)

    return ds, n_images
    

if __name__ == '__main__':
    # test out the dataset code with a barebones keras model
    BATCH_SIZE = 32
    ds, n_images = get_dataset(batch_size=BATCH_SIZE)
    # it = iter(ds)
    # print(next(it)[0].shape)

    model = tf.keras.Sequential([
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(32, activation='relu'),
        tf.keras.layers.Dense(2, activation='sigmoid'),
    ])

    model.compile(loss='categorical_crossentropy', optimizer='Adam', metrics=['Accuracy'])
    model.fit(ds, 
              steps_per_epoch=np.ceil(n_images/BATCH_SIZE),
              epochs=10)

