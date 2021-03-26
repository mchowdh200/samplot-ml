import argparse
import numpy as np
import tensorflow as tf
import tensorflow_addons as tfa

# data loading, CNN models
from utils import datasets, models

def main(args):
    model = tf.keras.models.load_model(args.model_path)
    dataset = datasets.DataWriter.get_basic_dataset(
        args.image_list, args.processes)
    dataset = dataset.batch(args.batch_size, drop_remainder=False)
    dataset = dataset.prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    
    for filenames, images in dataset:
        predictions = model(images)
        for file, pred in zip(filenames, predictions):
            f = os.path.splitext(os.path.basename(file.numpy()))[0]
            region = f.decode().split(args.delimiter)
            print(*region[:3], sep='\t', end='\t')
            print(*pred.numpy(), sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--model-path', dest='model_path', type=str, required=True,
        help='Path of trained model')
    parser.add_argument(
        '--image-list', dest='image_list', type=str, required=True,
        help='list of image file paths.')
    parser.add_argument(
        '--processes', dest='processes', type=int, required=True,
        help='number of simultaneous processes.')
    parser.add_argument(
        '--batch-size', dest='batch_size', type=int, required=True,
        help='number of images per patch.')
    parser.add_argument(
        '--delimiter', dest='delimiter', type=int, required=True,
        help='delimiter within image file name.')
    args = parser.parse_args()
    main(args)
