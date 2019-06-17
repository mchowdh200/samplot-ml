import os
# import pathlib
# import numpy as np
# import pandas as pd
import tensorflow as tf


# shape that we want to resize images to.  Original images are 2090 x 575
IMAGE_SHAPE = [262, 72, 3]

def get_filenames(data_dir='./data', training='train'):
    """
    Get filenames of images from provided root directory.
    """
    with open(f'{data_dir}/{training}.txt', 'r') as file:
        return [f'{data_dir}/cropped/{fname.rstrip()}' for fname in file]


def get_labels(filenames):
    """
    Filnames are in the following format:
        <chrm>_<start>_<end>_<sample>_<genotype>.png
    We will use the genotypes as the labels.
    """
    labels = [f.split('_')[4].split('.')[0] for f in filenames]
    label_to_index = {'ref': 0, 'het': 1, 'alt': 2}

    return [label_to_index[l] for l in labels]


def load_image(path):
    """
    Used to load image from provided filepath for use with a tensorflow dataset.
    Most images we have are 3 channel, but there are some that are 1/4 channels,
    so we just make all 3 channel then normalize to 0-1 range.
    """
    image = tf.io.read_file(path)
    image = tf.image.decode_png(image, channels=3)
    image = tf.image.resize(image, IMAGE_SHAPE[:2]) # original was 2090 x 575
    return image/255


def parse(x):
    """
    Used to reformat serialzed tensor to original shape
    """
    result = tf.io.parse_tensor(x, out_type=tf.float32)
    result = tf.reshape(result, IMAGE_SHAPE)
    return result


def get_dataset(data_dir='./data', batch_size=32, training='train'):
    """
    Create a tensorflow dataset that yields image/label pairs to a model
    """
    print('-'*80)
    print(f'getting {training} dataset...')
    print('-'*80)

    AUTOTUNE = tf.data.experimental.AUTOTUNE

    filenames = get_filenames(data_dir, training)
    labels = get_labels(filenames)
    n_images = len(filenames)

    # basic datasets from filenames and their labels
    image_ds = tf.data.Dataset.from_tensor_slices(filenames)
    image_ds = image_ds.map(load_image)
    image_ds = image_ds.map(tf.io.serialize_tensor)
    label_ds = tf.data.Dataset.from_tensor_slices(tf.cast(labels, tf.int64))

    # create the tfrecord file from the image dataset.  This takes a while
    # so I only want to do this if the file is not already present.
    if not os.path.isfile(f"{data_dir}/{training}.tfrec"):
        tfrec = tf.data.experimental.TFRecordWriter(f"{data_dir}/{training}.tfrec")
        tfrec.write(image_ds)

    # now load the dataset/parse the serialized tensors
    tfrds = tf.data.TFRecordDataset(f"{data_dir}/{training}.tfrec")
    tfrds = tfrds.map(parse, num_parallel_calls=AUTOTUNE)

    # combine with labels
    ds = tf.data.Dataset.zip((tfrds, label_ds))
    ds = ds.shuffle(buffer_size=n_images).repeat()
    ds = ds.batch(batch_size).prefetch(buffer_size=AUTOTUNE)

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
    print(filenames[:5])
    print(get_labels(filenames)[:5])


