import sys
import numpy as np
import pandas as pd
from sklearn.metrics import (precision_recall_curve, 
                             accuracy_score,
                             precision_score, 
                             recall_score, 
                             f1_score)
df = pd.read_csv(sys.stdin, sep='\t',
                 usecols=[6, 7, 8],
                 names=['pref', 'phet', 'palt'])

X = np.argmax(df.values, axis=1) > 0

print(sum(X)/len(X))

