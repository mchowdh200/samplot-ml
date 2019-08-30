import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--title', dest='title', type=str, default="")
parser.add_argument('-l', '--legend', dest='legend', nargs='+')
parser.add_argument('-f', '--filenames', dest='filenames', nargs='+')
parser.add_argument('-s', '--saveto', dest='saveto', 
                    required=False, default=None)
args = parser.parse_args()

for filename in args.filenames:
    df = pd.read_csv(
        filename, sep=' ', 
        names=['threshold', 'precision', 'recall', 'f1'])
    x = df.recall.values
    y = df.precision.values
    plt.plot(x, y , drawstyle='steps-post')

plt.title(args.title)
plt.xlabel('recall')
plt.ylabel('precision')
plt.legend(args.legend)

if args.saveto:
    plt.savefig(args.saveto)
else:
    plt.show()



