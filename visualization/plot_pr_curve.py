import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

title = sys.argv[1]
legend = [
    sys.argv[2],
    sys.argv[3],
]

for filename in sys.argv[4:]:
    df = pd.read_csv(
        filename, sep=' ', 
        names=['threshold', 'precision', 'recall', 'f1']
    )
    # x = df[df.precision > 0].recall.values
    # y = df[df.precision > 0].precision.values
    x = df.recall.values
    y = df.precision.values
    plt.plot(x, y , drawstyle='steps-post')

plt.xlabel('recall')
plt.ylabel('precision')
plt.legend(legend)
plt.show()



