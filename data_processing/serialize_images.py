import argparse
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from itertools import zip_longest, takewhile


# Original images are 2090 x 575
ORIG_SHAPE = [575, 2090, 3]

# amount to crop away during random cropping augmentation
RAND_CROP = np.array([5, 31, 0])
# RAND_CROP = np.array([0, 35, 0])

# we down scale each dimension by constant factor
SCALE_FACTOR = 8
IMAGE_SHAPE = np.array([np.ceil(ORIG_SHAPE[0]/SCALE_FACTOR).astype(int), 
                        np.ceil(ORIG_SHAPE[1]/SCALE_FACTOR).astype(int), 
                        3])


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks,
    Taken from python itertools docs
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def _int64_feature(value):
    if not isinstance(value, list):
        value = [value]
    return tf.train.Feature(int64_list=tf.train.Int64List(value=value))


def _bytes_feature(value):
    # If the value is an eager tensor BytesList 
    # won't unpack a string from an EagerTensor.
    if isinstance(value, tf.Tensor):
        value = value.numpy() 
    elif not isinstance(value, bytes):
        value = value.encode()
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))


def get_filenames(data_dir, training='train'):
    """
    Get filenames of images from provided root directory.
    """
    with open(f'{data_dir}/{training}.txt', 'r') as file:
        return [f'{data_dir}/crop/{fname.rstrip()}' for fname in file]


def get_labels(filenames, num_classes=3):
    """
    Filnames are in the following format:
        <chrm>_<start>_<end>_<sample>_<genotype>.png
    We will use the genotypes as the labels.
    Additionally, we apply label smoothing to the categorical labels
    controlled by the parameter eps.
    """
    labels = [f.split('_')[-1].split('.')[0].lower() for f in filenames]

    if num_classes == 3:
        label_to_index = {'ref': 0, 'het': 1, 'alt': 2}
    else:
        label_to_index = {'ref': 0, 'del': 1}
    return [label_to_index[l] for l in labels]


def load_image(path):
    """
    Used to load image from provided filepath for use with a tensorflow dataset.
    Most images we have are 3 channel, but there are some that are 1/4 channels,
    so we just make all 3 channel then normalize to 0-1 range.
    """
    image = tf.io.read_file(path)
    image = tf.image.decode_png(image, channels=3)
    image = tf.image.convert_image_dtype(image, tf.float32)
    image = tf.image.resize(image, IMAGE_SHAPE[:2])
    return image


def parse_image(x):
    """
    Used to reformat serialzed images to original shape
    """
    result = tf.io.parse_tensor(x, out_type=tf.float32)
    result = tf.reshape(result, IMAGE_SHAPE)
    return result


def serialize_example(filename, label):
    """
    given a filepath to an image (label contained in filename),
    serialize the (image, label)
    """
    example = tf.train.Example(
        features=tf.train.Features(
            feature={
                'filename': _bytes_feature(filename),
                'image': _bytes_feature(
                    tf.io.serialize_tensor(load_image(filename))),
                'label': _int64_feature(label),
            }
        )
    )
    return example.SerializeToString()


def to_tfrecords(data_dir, imgs_per_record=1000,
                 training='train', num_classes=3):
    """
    Write train and/or val set to a set of TFRecords
    """
    filenames = get_filenames(data_dir, training)
    labels = get_labels(filenames, num_classes=num_classes)
    n_images = len(filenames)
    assert len(filenames) == len(labels)

    Parallel(n_jobs=-1)(
        delayed(write_batch)(data_dir=data_dir,
                             batch_index=i,
                             file_label_pairs=takewhile(
                                 lambda x: x is not None,
                                 file_label_pairs),
                             training=training) 
        for i, file_label_pairs in enumerate(
            grouper(zip(filenames, labels), imgs_per_record))
    )


def write_batch(data_dir, batch_index, file_label_pairs, training='train'):
    """
    Write a single batch of images to TFRecord format
    """
    with tf.io.TFRecordWriter(
        f"{data_dir}/{training}/{training}_{batch_index:05d}.tfrec") as writer:
        for file_label in file_label_pairs:
            filename, label = file_label
            serialized_example = serialize_example(filename, label)
            writer.write(serialized_example)


def parse_serialized_example(serialized_example):
    """
    Given a serialized example with the below format, extract/deserialize
    the image and (one-hot) label
    """
    features = {
        'filename': tf.io.FixedLenFeature((), tf.string),
        'image': tf.io.FixedLenFeature((), tf.string),
        'label': tf.io.FixedLenFeature((), tf.int64)
    }

    example = tf.io.parse_single_example(serialized_example, features)
    image = parse_image(example['image'])
    label = example['label']
    return image, tf.one_hot(label, depth=3, dtype=tf.int64)
        

if __name__ == "__main__":
    # Create the sharded dataset
    # to_tfrecords(
    #     data_dir='/home/murad/Desktop/1kgtest',
    #     training='train',
    #     imgs_per_record=1000,
    #     num_classes=3)

    # Test loading an image to make sure we can recover it
    files = tf.data.Dataset.list_files("/home/murad/Desktop/1kgtest/train/train_*.tfrec")
    dataset = tf.data.TFRecordDataset(files)
    dataset = dataset.map(parse_serialized_example).shuffle(buffer_size=1000)

    count=0
    for image, label in dataset:
        count += 1
        print(label.numpy())
        plt.imshow(image.numpy())
        plt.title(f'{np.argmax(label.numpy())}')
        plt.show()
        exit()
    print(count)
