import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (precision_recall_curve, 
                             accuracy_score,
                             precision_score, 
                             recall_score, 
                             f1_score)

# Input is output of bedtools -wb with svplaudit 
# mean_scores bed and CNN predictions bed.
# 
# The format of the input is (tab delimitted):
# $chr $start $end DEL $svp_score $chr $start $end $P_ref $P_het $P_alt
# 0    1      2    3   4          5    6      7    8      9      10
df = pd.read_csv(sys.stdin, sep='\t',
                 usecols=[4, 8, 9, 10],
                 names=['svplaudit', 'pref', 'phet', 'palt'])

X = df.values

# keep unambiguous svplaudit scores as per the paper's methods
X = X[(X[:, 0] < 0.2) | (X[:, 0] > 0.8)]

svp_labels = (X[:, 0] > 0.5).astype(int)
CNN_pred = (np.argmax(X[:, 1:], axis=1) > 0).astype(int)

print(len(svp_labels))
print("accuracy: ", accuracy_score(svp_labels, CNN_pred))
print("precision: ", precision_score(svp_labels, CNN_pred))
print("recall: ", recall_score(svp_labels, CNN_pred))
print("f1: ", f1_score(svp_labels, CNN_pred))

