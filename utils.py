import os
import numpy as np
import tensorflow as tf
from sklearn.metrics import classification_report, confusion_matrix


# -----------------------------------------------------------------------------
# Data utils
# -----------------------------------------------------------------------------

# shape that we want to resize images to.  Original images are 2090 x 575
ORIG_SHAPE = [2090, 575, 3]
SCALE_FACTOR = 8
IMAGE_SHAPE = [np.ceil(ORIG_SHAPE[0]/SCALE_FACTOR).astype(int), 
               np.ceil(ORIG_SHAPE[1]/SCALE_FACTOR).astype(int), 
               3]

def get_filenames(data_dir, training='train'):
    """
    Get filenames of images from provided root directory.
    """
    with open(f'{data_dir}/{training}.txt', 'r') as file:
        return [f'{data_dir}/crop/{fname.rstrip()}' for fname in file]


def get_labels(filenames):
    """
    Filnames are in the following format:
        <chrm>_<start>_<end>_<sample>_<genotype>.png
    We will use the genotypes as the labels.
    Additionally, we apply label smoothing to the categorical labels
    controlled by the parameter eps.
    """
    labels = [f.split('_')[5].split('.')[0] for f in filenames]
    label_to_index = {'ref': 0, 'het': 1, 'alt': 2}


    return [tf.keras.utils.to_categorical(label_to_index[l], num_classes=3) 
            for l in labels]


def load_image(path):
    """
    Used to load image from provided filepath for use with a tensorflow dataset.
    Most images we have are 3 channel, but there are some that are 1/4 channels,
    so we just make all 3 channel then normalize to 0-1 range.
    """
    image = tf.io.read_file(path)
    image = tf.image.decode_png(image, channels=3)
    image = tf.image.resize(image, IMAGE_SHAPE[:2])
    return image/255


def parse(x):
    """
    Used to reformat serialzed tensor to original shape
    """
    result = tf.io.parse_tensor(x, out_type=tf.float32)
    result = tf.reshape(result, IMAGE_SHAPE)
    return result


def get_dataset(data_dir, batch_size=32,
                training='train', shuffle=True,
                augmentation=False,
                return_labels=False):
    """
    Create a tensorflow dataset that yields image/label pairs to a model.
    If return_labels is True, then we also return a dataset consisting of just labels
    (note: these will not be shuffled, so the labels will only be valid if shuffle is False)
    """
    print(f'GETTING {training.upper()} DATASET')

    AUTOTUNE = tf.data.experimental.AUTOTUNE

    filenames = get_filenames(data_dir, training)
    labels = get_labels(filenames)
    n_images = len(filenames)

    assert len(filenames) == len(labels)

    # basic datasets from filenames and their labels
    image_ds = tf.data.Dataset.from_tensor_slices(filenames)
    image_ds = image_ds.map(load_image, num_parallel_calls=AUTOTUNE)
    image_ds = image_ds.map(tf.image.per_image_standardization, 
                            num_parallel_calls=AUTOTUNE)

    label_ds = tf.data.Dataset.from_tensor_slices(tf.cast(labels, tf.int64))

    if augmentation:
        image_ds = image_ds.concatenate(
            image_ds.map(tf.image.flip_left_right,
                         num_parallel_calls=AUTOTUNE))
        label_ds = label_ds.concatenate(label_ds)
        n_images *= 2
    image_ds = image_ds.map(tf.io.serialize_tensor)


    # create the tfrecord file from the image dataset.  This takes a while
    # so I only want to do this if the file is not already present.
    if not os.path.isfile(f"{data_dir}/{training}.tfrec"):
        print("generating TFRecord...")
        tfrec = tf.data.experimental.TFRecordWriter(f"{data_dir}/{training}.tfrec")
        tfrec.write(image_ds)

    # now load the dataset/parse the serialized tensors
    tfrds = tf.data.TFRecordDataset(
        f"{data_dir}/{training}.tfrec",
        # num_parallel_reads=os.cpu_count(),
    )
    tfrds = tfrds.map(parse, num_parallel_calls=AUTOTUNE)



    # combine with labels
    ds = tf.data.Dataset.zip((tfrds, label_ds))
    if shuffle:
        ds = ds.shuffle(buffer_size=10000)
    ds = ds.repeat()
    ds = ds.batch(batch_size, drop_remainder=True)
    # ds = ds.prefetch(buffer_size=n_images//(4*batch_size))
    ds = ds.prefetch(buffer_size=AUTOTUNE)

    if return_labels:
        if shuffle:
            raise ValueError('Cannot have shuffle==True when returning labels')
        return ds, label_ds, n_images

    return ds, n_images


# -----------------------------------------------------------------------------
# Model utils
# -----------------------------------------------------------------------------
def display_prediction(path, model):
    """
    Given path of an image and a trained model, returns prediction
    probability distribution over classes class.
    """
    # load image and make a prediction from the saved model.
    img = load_image(path)
    img = tf.image.per_image_standardization(img)
    img = tf.expand_dims(img, axis=0)
    pred = model(img).numpy()

    return pred[0]


def evaluate_model(model, data_dir, batch_size=80):
    """
    Take a dataset with (image, label) pairs along with a trained model and
    evaluate per class metrics.
    """
    # tfmodel = tf.keras.models.load_model(model)

    test_ds, label_ds, n = get_dataset(batch_size=batch_size, data_dir=data_dir,
                                       training='val', shuffle=False, 
                                       return_labels=True)

    if n % 2:
        test_ds = test_ds.take(n-1)
        label_ds = label_ds.take(n-1)
        n-=1


    assert n % batch_size == 0, \
        f'Batch size of {batch_size} does not evenly divide into size of data ({n}).'

    y_true = np.argmax(np.array([x.numpy() for x in label_ds.take(n)]), axis=1)
    y_pred = np.argmax(model.predict(test_ds, steps=np.ceil(n/batch_size)), axis=1)
    print(confusion_matrix(y_true, y_pred))
    print(classification_report(y_true, y_pred))


