import os
import functools
# import cv2
import numpy as np
import tensorflow as tf
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# -----------------------------------------------------------------------------
# Data utils
# -----------------------------------------------------------------------------

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


    return [tf.keras.utils.to_categorical(label_to_index[l], num_classes=num_classes) 
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
    return image


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
                return_labels=False,
                num_classes=3,
                processes=1):
    print("PROCESSES=", processes)
    """
    Create a tensorflow dataset that yields image/label pairs to a model.
    If return_labels is True, then we also return a dataset consisting of just labels
    (note: these will not be shuffled, so the labels will only be valid if shuffle is False)
    """
    print(f'GETTING {training.upper()} DATASET')

    #TODO there is a bug in tensorflow 2.0 where using AUTOTUNE
    # in map() causes linear increase in memory useage until the
    # program gets killed.  Check back later when the bug is fixed.
    AUTOTUNE = tf.data.experimental.AUTOTUNE

    filenames = get_filenames(data_dir, training)
    labels = get_labels(filenames, num_classes=num_classes)
    n_images = len(filenames)

    assert len(filenames) == len(labels)

    # basic datasets from filenames and their labels
    image_ds = tf.data.Dataset.from_tensor_slices(filenames)
    image_ds = image_ds.map(load_image, num_parallel_calls=processes)
    image_ds = image_ds.map(tf.image.per_image_standardization, 
                            num_parallel_calls=processes)

    label_ds = tf.data.Dataset.from_tensor_slices(tf.cast(labels, tf.int64))

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
    tfrds = tfrds.map(parse, num_parallel_calls=processes)
    
    # create flipped/cropped versions of images on the fly
    if augmentation:
        tfrds = tfrds.map(tf.image.random_flip_left_right)
        tfrds = tfrds.map(
            functools.partial(
                tf.image.random_crop, size=IMAGE_SHAPE - RAND_CROP))

    # combine with labels
    ds = tf.data.Dataset.zip((tfrds, label_ds))
    if shuffle:
        ds = ds.shuffle(buffer_size=10000)
    ds = ds.repeat()
    ds = ds.batch(batch_size, drop_remainder=False)
    # ds = ds.prefetch(buffer_size=n_images//(4*batch_size))
    ds = ds.prefetch(buffer_size=1000)

    if return_labels:
        if shuffle:
            raise ValueError('Cannot have shuffle==True when returning labels')
        return ds, label_ds, n_images

    return ds, n_images


# -----------------------------------------------------------------------------
# Model utils
# -----------------------------------------------------------------------------
def display_prediction(path, model, augmentation=False, n_aug=15):
    """
    Given path of an image and a trained model, returns prediction
    probability distribution over classes.  Optionally computes prediction
    with test time prediction where we average the predictions of n_aug derived
    images with random lr flips and random crops performed.
    """
    # load image and make a prediction from the saved model.
    img = load_image(path)
    img = tf.image.per_image_standardization(img)

    if augmentation:
        transformed_images = []
        for _ in range(n_aug):
            transformed_images.append(
                tf.image.random_crop(tf.image.random_flip_left_right(img),
                                     size=IMAGE_SHAPE - RAND_CROP))
        img = tf.stack(transformed_images, axis=0)
    else:
        img = tf.expand_dims(img, axis=0)

    pred = model(img)

    if augmentation:
        # majority vote prediction -- 
        x = tf.argmax(pred, axis=-1)
        x = tf.one_hot(x, depth=3)
        x = tf.reduce_sum(x, axis=0)
        return x.numpy()


    return pred.numpy()[0]

def grad_cam(path, model, final_conv, pre_softmax, out_dir=None):
    """
    Grad cam visualization technique for convolutional neural networks.
    https://arxiv.org/pdf/1610.02391.pdf
    - path: image path.
    - model: keras model.
    - final_conv: name of the final convolutional layer.
    - pre_softmax: name of the Dense layer before softmax

    Variables named to match the procedure in the paper
    - y: class prediction score (before softmax)
    - A: Activations of final conv layer
    - dy_dA: derivative of prediction (for a given class) w.r.t. the activations
    - alpha: per channel weights for the activations (global average pooled from dy_dA)
    - L: grad cam weights
    """
    img = tf.io.read_file(path) #load_image(path)
    img = tf.image.decode_png(img, channels=3)

    normed_img = tf.image.resize(img, IMAGE_SHAPE[:2])
    normed_img=tf.image.per_image_standardization(normed_img)
    normed_img=tf.expand_dims(normed_img, axis=0)

    submodel = tf.keras.Model(inputs=model.inputs,
                              outputs=[model.get_layer(pre_softmax).output,
                                       model.get_layer(final_conv).output,
                                       model.layers[-1].output])

    # get the gradient of the pre softmax score w.r.t. the final conv layer
    # feature maps and create the grad-cam
    with tf.GradientTape(persistent=True) as tape:
        y, A , pred = submodel(normed_img)
        dy_dA = [tape.gradient(y[:, i], A) for i in range(3)]


    # convert the heatmap and image to 0-255 range
    img = np.uint8(img)
    plt.close()
    plt.rcParams["figure.figsize"] = (10, 20)
    gs = matplotlib.gridspec.GridSpec(4, 1)
    plt.subplot(gs[0, 0])
    plt.imshow(img)
    plt.tick_params(top=False, bottom=False, 
                    left=False, right=False, 
                    labelleft=False, labelbottom=False)
    plt.xlabel('Original image')
    plt.title(f"[ref, het, alt] = {pred.numpy()[0]}")

    label_dict = {0: 'Homozygous Reference',
                  1: 'Heterozygous',
                  2: 'Homozygous Alternate'}

    for i, x in enumerate(dy_dA):
        alpha = tf.keras.layers.GlobalAveragePooling2D()(x)
        L = tf.keras.layers.ReLU()(tf.math.reduce_sum(alpha*A, axis=3)).numpy()

        # scale the weights and remove the batch dimension
        heatmap = L[0]/np.max(L)

        # image dimensions are swapped in the resize function (why???)
        heatmap = cv2.resize(heatmap, (img.shape[1], img.shape[0]))
        
        heatmap = np.uint8(heatmap*255)
        heatmap = cv2.applyColorMap(heatmap, cv2.COLORMAP_JET)
        cam = cv2.addWeighted(img, 0.6, heatmap, 0.4, 0)

        plt.subplot(gs[i + 1, 0])
        plt.imshow(cam)
        plt.tick_params(top=False, bottom=False, 
                        left=False, right=False, 
                        labelleft=False, labelbottom=False)
        plt.xlabel(label_dict[i])

    plt.tight_layout()
    if out_dir:
        plt.savefig(
            out_dir + '/' + os.path.basename(path).split('.')[0] + '.heatmap.png')
    else:
        plt.show()


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



# if __name__ == '__main__':

#     model = tf.keras.models.load_model('./saved_models/CNN_7_16.h5')
#     grad_cam(
#         path='data/giab/crop/8_146042946_146043066_DEL.png',
#         model=model, 
#         # final_conv='leaky_re_lu_20',
#         final_conv='leaky_re_lu_39',
#         pre_softmax='dense_1')




