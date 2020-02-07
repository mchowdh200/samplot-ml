# import functools
import os
import sys
import functools
from itertools import zip_longest, takewhile

import numpy as np
import tensorflow as tf
from joblib import Parallel, delayed


# Original images are 2090 x 575
ORIG_SHAPE = [575, 2090, 3]

# amount to crop away during random cropping augmentation
# RAND_CROP = np.array([5, 31, 0])
RAND_CROP = np.array([0, 31, 0])

# we down scale each dimension by constant factor
SCALE_FACTOR = 8
IMAGE_SHAPE = np.array([np.ceil(ORIG_SHAPE[0]/SCALE_FACTOR).astype(int), 
                        np.ceil(ORIG_SHAPE[1]/SCALE_FACTOR).astype(int), 
                        3])


class DataWriter:
    def __init__(self, data_list, out_dir, training, num_classes=3):
        self.out_dir = out_dir
        self.training = training
        self.num_classes=num_classes
        self.filenames = [fname.rstrip() for fname in open(data_list)]
        self.labels = DataWriter._get_labels(self.filenames, num_classes)
        assert len(self.filenames) == len(self.labels)

    @staticmethod
    def _grouper(iterable, n, fillvalue=None):
        """
        Collect data into fixed-length chunks or blocks,
        Taken from python itertools docs
        """
        # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
        args = [iter(iterable)] * n
        return zip_longest(fillvalue=fillvalue, *args)

    @staticmethod
    def _int64_feature(value):
        if not isinstance(value, list):
            value = [value]
        return tf.train.Feature(int64_list=tf.train.Int64List(value=value))

    @staticmethod
    def _bytes_feature(value):
        # If the value is an eager tensor BytesList 
        # won't unpack a string from an EagerTensor.
        if isinstance(value, tf.Tensor):
            value = value.numpy() 
        elif not isinstance(value, bytes):
            value = value.encode()
        return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

    @staticmethod
    def _get_labels(filenames, num_classes=3):
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

    @staticmethod
    def _load_image(path):
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

    @staticmethod
    def get_basic_dataset(image_list, num_processes):
        """
        return a tensorflow dataset consisting of just filenames
        and images (no labels)
        """
        filename_ds = tf.data.Dataset.from_tensor_slices(
            [filename.rstrip() for filename in open(image_list, 'r')])
        image_ds = filename_ds.map(
            DataWriter._load_image,
            num_parallel_calls=num_processes).map(
                tf.image.per_image_standardization,
                num_parallel_calls=num_processes
            )

        return tf.data.Dataset.zip((filename_ds, image_ds))


    @staticmethod
    def _serialize_example(filename, label):
        """
        given a filepath to an image (label contained in filename),
        serialize the (image, label)
        """
        example = tf.train.Example(
            features=tf.train.Features(
                feature={
                    'filename': DataWriter._bytes_feature(filename),
                    'image': DataWriter._bytes_feature(
                        tf.io.serialize_tensor(DataWriter._load_image(filename))),
                    'label': DataWriter._int64_feature(label),
                }
            )
        )
        return example.SerializeToString()

    @staticmethod
    def _write_batch(out_dir, batch_index, file_label_pairs, training):
        """
        Write a single batch of images to TFRecord format
        """
        with tf.io.TFRecordWriter(
            f"{out_dir}/{training}/{training}_{batch_index:05d}.tfrec") as writer:
            for file_label in file_label_pairs:
                filename, label = file_label
                serialized_example = DataWriter._serialize_example(filename, label)
                writer.write(serialized_example)

    def to_tfrecords(self, imgs_per_record=1000):
        """
        Write train and/or val set to a set of TFRecords
        """
        Parallel(n_jobs=-1)(
            delayed(self._write_batch)(self.out_dir,
                                 batch_index=i,
                                 file_label_pairs=takewhile(
                                     lambda x: x is not None,
                                     file_label_pairs),
                                 training=self.training) 
            for i, file_label_pairs in enumerate(
                DataWriter._grouper(zip(self.filenames, self.labels), imgs_per_record))
        )


class DataReader:
    def __init__(self, 
                 data_list,  # list of original images in the dataset
                 tfrec_list, # list of tfrecords in the dataset (s3 or local)
                 num_processes,
                 batch_size,
                 augmentation=False):

        self.data_list = data_list
        self.tfrec_list = tfrec_list
        self.num_processes = num_processes
        self.batch_size = batch_size
        self.augmentation = augmentation

    @staticmethod
    def _parse_image(x):
        """
        Used to reformat serialzed images to original shape
        """
        result = tf.io.parse_tensor(x, out_type=tf.float32)
        result = tf.reshape(result, IMAGE_SHAPE)
        return result

    @staticmethod
    def _parse_serialized_example(serialized_example, augmentation=False):
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

        # read image and perform transormations
        image = DataReader._parse_image(example['image'])
        image = tf.image.per_image_standardization(image)

        # if augmentation:
        #     image = tf.image.random_flip_left_right(image)
        #     image = tf.image.random_crop(image, size=IMAGE_SHAPE - RAND_CROP)


        label = example['label']
        return image, tf.one_hot(label, depth=3, dtype=tf.int64)
    
    def get_dataset(self):
        n_images = len(open(self.data_list).readlines())

        # we don't need examples to be loaded in order (better speed)
        options = tf.data.Options()
        options.experimental_deterministic = False

        with open(self.tfrec_list) as f:
            files = [filename.rstrip() for filename in f]
            dataset = tf.data.Dataset.from_tensor_slices(files) \
                    .shuffle(len(files)) \
                    .with_options(options)


        dataset = dataset.interleave(
            tf.data.TFRecordDataset,
            cycle_length=self.num_processes,
            num_parallel_calls=self.num_processes) \

        dataset = dataset.map(
            functools.partial(DataReader._parse_serialized_example, augmentation=self.augmentation),
            num_parallel_calls=self.num_processes) \
                .repeat() \
                .shuffle(buffer_size=1000) \
                .batch(self.batch_size, drop_remainder=False) \
                .prefetch(buffer_size=tf.data.experimental.AUTOTUNE)


        # TODO incorporate into _parse_serialized_example so that I can
        # operate on just the image and don't have to do anything funky
        # with the labels.
        # if self.augmentation:
        #     dataset = dataset.map(
        #         tf.image.random_flip_left_right,
        #         num_parallel_calls=tf.data.experimental.AUTOTUNE)
        #     dataset = dataset.map(
        #         functools.partial(
        #             tf.image.random_crop, size=IMAGE_SHAPE - RAND_CROP),
        #         num_parallel_calls=tf.data.experimental.AUTOTUNE)

        return dataset, n_images
                         
        

if __name__ == "__main__":

    train_list = sys.argv[1]
    val_list = sys.argv[2]
    out_dir = sys.argv[3]

    # we're just writing datasets, so don't use the gpu
    # (runs out of gpu memory otherwise)
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

    print('Writing training set to tfrecords')
    DataWriter(
        data_list = train_list,
        out_dir = out_dir,
        training = 'train').to_tfrecords(imgs_per_record=1000)

    print('Writing validation set to tfrecords')
    DataWriter(
        data_list = val_list,
        out_dir = out_dir,
        training = 'val').to_tfrecords(imgs_per_record=1000)
