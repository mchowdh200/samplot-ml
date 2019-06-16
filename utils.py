import pathlib
import numpy as np
import pandas as pd
import tensorflow as tf


def get_filenames(data_dir='./data/cropped'):
    """
    Get filenames of images from provided root directory.
    """
    data_root = pathlib.Path(data_dir)
    return [str(path) for path in list(data_root.glob('*.png'))]


def get_labels(filenames):
    """
    Filnames are in the following format:
        <chrm>_<start>_<end>_<sample>_<genotype>.png
    We will use the genotypes as the labels.
    """
    labels = [f.split('_')[4].split('.')[0] for f in filenames]
    label_to_index = {'ref': 0, 'het': 1, 'alt': 2}

    return [label_to_index(l) for l in labels]


def load_image(path, label):
    """
    Used to load image from provided filepath for use with a tensorflow dataset.
    Most images we have are 3 channel, but there are some that are 1/4 channels,
    so we just make all 3 channel then normalize to 0-1 range.  We pass through
    the image label so we can yield it in a tensorflow datset
    """
    image = tf.io.read_file(path)
    return tf.image.decode_png(image, channels=3)/255, label


def get_dataset(data_dir='./data/cropped', batch_size=32):
    """
    Create a tensorflow dataset that yields image/label pairs to a model
    """

    # df = pd.read_csv(f"{data_dir}/mean_score.bed", sep='\t',
    #                  names=['chrom', 'start', 'end', 'SVTYPE', 'score'])

    # n_images = len(df)

    # df['filename'] = df.apply(get_filenames, axis=1)
    # df['label'] = df.apply(get_labels, axis=1)

    AUTOTUNE = tf.data.experimental.AUTOTUNE

    filenames = get_filenames(data_dir)
    labels = get_labels(filenames)

    ds = tf.data.Dataset.from_tensor_slices(
        (filenames, tf.cast(labels, tf.int64)))

    ds = ds.map(load_image)
    ds = ds.cache()
    ds = ds.apply(tf.data.experimental.shuffle_and_repeat(buffer_size=n_images))
    ds = ds.batch(batch_size)
    ds = ds.prefetch(buffer_size=AUTOTUNE)

    return ds, n_images
    

if __name__ == '__main__':
    # # test out the dataset code with a barebones keras model
    # BATCH_SIZE = 32
    # ds, n_images = get_dataset(batch_size=BATCH_SIZE)
    # # it = iter(ds)
    # # print(next(it)[0].shape)

    # model = tf.keras.Sequential([
    #     tf.keras.layers.Flatten(),
    #     tf.keras.layers.Dense(32, activation='relu'),
    #     tf.keras.layers.Dense(2, activation='sigmoid'),
    # ])

    # model.compile(loss='categorical_crossentropy', optimizer='Adam', metrics=['Accuracy'])
    # model.fit(ds, 
    #           steps_per_epoch=np.ceil(n_images/BATCH_SIZE),
    #           epochs=10)

    filenames = get_filenames()
    print(get_labels(filenames)[:5])



