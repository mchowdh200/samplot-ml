usage: run.py [-h] {predict,evaluate,train} ...

optional arguments:
  -h, --help            show this help message and exit

Subcommands:
  {predict,evaluate,train}
    predict             Use a trained model to classify an image.
    evaluate            Evaluate a trained model on a labelled dataset
    train               Train a new model.


usage: run.py train [-h] [--batch-size BATCH_SIZE] [--epochs EPOCHS]
                    [--model-type MODEL_TYPE] [--learning-rate LR]
                    [--save-to SAVE_TO]

optional arguments:
  -h, --help            show this help message and exit
  --batch-size BATCH_SIZE, -b BATCH_SIZE
                        Number of images to feed to model at a time. (default:
                        80)
  --epochs EPOCHS, -e EPOCHS
                        Max number of epochs to train model. (default: 100)
  --model-type MODEL_TYPE, -m MODEL_TYPE
                        Type of model to train. (default: CNN)
  --learning-rate LR, -lr LR
                        Learning rate for optimizer. (default: 0.0001)
  --save-to SAVE_TO, -s SAVE_TO
                        filename if you want to save your trained model.
                        (default: None)

usage: run.py evaluate [-h] --model-path MODEL_PATH --model-type {CNN}
                       [--use-h5] [--batch-size BATCH_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  --model-path MODEL_PATH, -mp MODEL_PATH
                        Path of trained model (before the dot '.') (default:
                        None)
  --model-type {CNN}, -mt {CNN}
                        Type of model to load. (default: CNN)
  --use-h5, -h5
  --batch-size BATCH_SIZE, -b BATCH_SIZE
                        Number of images to feed to model at a time. (default:
                        80)

usage: run.py predict [-h] --model-path MODEL_PATH --model-type {CNN}
                      [--use-h5] --image IMAGE

optional arguments:
  -h, --help            show this help message and exit
  --model-path MODEL_PATH, -mp MODEL_PATH
                        Path of trained model (default: None)
  --model-type {CNN}, -mt {CNN}
                        Type of model to load. (default: CNN)
  --use-h5, -h5
  --image IMAGE, -i IMAGE
                        Path of image (default: None)
