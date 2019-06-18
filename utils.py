import os
import numpy as np
import tensorflow as tf
from sklearn.metrics import classification_report, confusion_matrix


# -----------------------------------------------------------------------------
# Data utils
# -----------------------------------------------------------------------------

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


def get_dataset(data_dir='./data', batch_size=32,
                training='train', shuffle=True,
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
    if shuffle:
        ds = ds.shuffle(buffer_size=n_images//4)
    ds = ds.repeat()
    ds = ds.batch(batch_size).prefetch(buffer_size=AUTOTUNE)

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
    Given path of an image and a trained model, plot the image along with its
    predicted class.
    """
    tfmodel = tf.keras.models.load_model(model)
    
    # load image and make a prediction from the saved model.
    img = load_image(path)
    img = tf.expand_dims(img, axis=0)
    pred = tfmodel(img).numpy()

    # TODO plot the original image and display prediction distribution
    return list(pred[0])


def evaluate_model(model, batch_size=80):
    """
    Take a dataset with (image, label) pairs along with a trained model and
    evaluate per class metrics.
    """
    tfmodel = tf.keras.models.load_model(model)

    test_ds, label_ds, n = get_dataset(batch_size=batch_size, training='test', 
                                       shuffle=False, return_labels=True)

    assert n % batch_size == 0, \
        f'Batch size of {80} does not evenly divide into size of data ({n}).'

    y_true = np.array([x.numpy() for x in label_ds.take(n)])
    y_pred = np.argmax(tfmodel.predict(test_ds, steps=np.ceil(n/batch_size)), axis=1)
    print(confusion_matrix(y_true, y_pred))
    print(classification_report(y_true, y_pred))




# if __name__ == '__main__':

    # just picked one of each type of image
    # test_ref = './data/cropped/10_100043404_100047636_NA07000_ref.png'
    # test_het = './data/cropped/10_100375789_100378989_NA18505_het.png'
    # test_alt = './data/cropped/10_100688639_100702031_HG00638_alt.png'

    # model = tf.keras.models.load_model('./saved_models/my_model.h5')
    
    # filenames = get_filenames(training='test')
    # with open('file_predictions.txt', 'w') as f:
    #     f.write('Filename\tPrediction Distribution\n')
    #     for fname in filenames:
    #         f.write(f'{fname}\t{display_prediction(fname, model)}\n')
            
    # evaluate_model(model)


