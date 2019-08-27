import sys
import tensorflow as tf
import utils

model = tf.keras.models.load_model('./saved_models/CNN_7_16.h5')

for line in sys.stdin: # cat image list
    utils.grad_cam(
        path=line.rstrip(),
        out_dir='./analysis/comparison/images',
        model=model,
        final_conv='leaky_re_lu_39',
        pre_softmax='dense_1'
    )

