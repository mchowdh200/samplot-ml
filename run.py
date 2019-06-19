#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import tensorflow as tf
import utils
import models

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# TODO refactor evaluate and predict to pick the model type

model_index = {
    'CNN' : models.CNN
}


def get_model(model_type, model_params=None, compile_params=None):
    """
    Returns a compiled model with given model/compile parameters
    """
    if model_params:
        model = model_type(**model_params)
    else:
        model = model_type()

    # model.compile(loss='sparse_categorical_crossentropy', 
    #               optimizer=tf.optimizers.Adam(lr=5e-4, amsgrad=True), 
    #               metrics=['SparseCategoricalAccuracy'])
    model.compile(**compile_params)

    return model



# -----------------------------------------------------------------------------
# Subcommand functions
# -----------------------------------------------------------------------------
def predict(args):
    print(utils.display_prediction(args.image, args.model))

def evaluate(args):
    utils.evaluate_model(args.model, args.batch_size)

def train(args):

    # load data
    training_set, n_train = utils.get_dataset(
        batch_size=args.batch_size, training='train')
    test_set, n_test = utils.get_dataset(
        batch_size=args.batch_size, training='test')

    # setup training
    callbacks = [tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss',patience=3),
                 tf.keras.callbacks.EarlyStopping(monitor='val_loss',patience=5,
                                                  restore_best_weights=True, verbose=1)]

    # TODO eventually I'd like to optionally be able to load these from a JSON
    model_params = None
    compile_params = dict(
        loss='sparse_categorical_crossentropy',
        optimizer=tf.optimizers.Adam(lr=5e-4, amsgrad=True),
        metrics=['SparseCategoricalAccuracy'])
    model = get_model(model_index[args.model_type], model_params, compile_params)
    model.fit(training_set,
              steps_per_epoch=np.ceil(n_train/args.batch_size),
              validation_data=test_set,
              validation_steps=np.ceil(n_test/args.batch_size),
              epochs=args.epochs,
              callbacks=callbacks)

    print(model.summary())

    if args.save_to:
        model.save_weights(f"./saved_models/{args.save_to}_{type(model).__name__}",
                           save_format='tf')


# -----------------------------------------------------------------------------
# Get arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(title='Subcommands')

# prediction subcommand -------------------------------------------------------
predict_parser = subparsers.add_parser(
    'predict', help='Use a trained model to classify an image.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
predict_parser.add_argument(
    '--model', '-m', dest='model', type=str, required=True,
    help='Path of trained model')
predict_parser.add_argument(
    '--image', '-i', dest='image', type=str, required=True,
    help='Path of image')
predict_parser.set_defaults(
    func=predict)

# evaluation subcommand -------------------------------------------------------
eval_parser = subparsers.add_parser(
    'evaluate', help='Evaluate a trained model on a labelled dataset',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
eval_parser.add_argument(
    '--model', '-m', dest='model', type=str, required=True,
    help='Path of trained model')
eval_parser.add_argument(
    '--batch-size', '-b', dest='batch_size', type=int, required=False,
    default=80, help='Number of images to feed to model at a time.')
eval_parser.set_defaults(
    func=evaluate)

# training subcommand ---------------------------------------------------------
train_parser = subparsers.add_parser(
    'train', help='Train a new model.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
train_parser.add_argument(
    '--batch-size', '-b', dest='batch_size', type=int, required=False,
    default=80, help='Number of images to feed to model at a time.')
train_parser.add_argument(
    '--epochs', '-e', dest='epochs', type=int, required=False,
    default=100, help='Max number of epochs to train model.')
train_parser.add_argument(
    '--model-type', '-m', dest='model_type', type=str, required=False,
    default='CNN', help='Type of model to train.'
)
train_parser.add_argument(
    '--save-to', '-s', dest='save_to', type=str, required=False,
    default=None, help='filename if you want to save your trained model.')
train_parser.set_defaults(
    func=train)

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args.func(args)
