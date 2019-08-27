#!/bin/bash

python plot_pr_curve.py -t "HG002 Precision-Recall Curves" \
    -l "CNN" "Duphold DHFFC" \
    -f "~/data/giab/VCF/CNN_8_6/pr_curve/pr_curve.txt" \
       "~/data/giab/VCF/DHFFC/pr_curve/pr_curve.txt" \
    -s "../figures/HG002-notier1-pr.png"
python plot_pr_curve.py -t "HG00514 Precision-Recall Curves" \
    -l "CNN" "Duphold DHFFC" \
    -f "~/data/HG00514/VCF/CNN_8_6/pr_curve/pr_curve.txt" \
       "~/data/HG00514/VCF/duphold/pr_curve/pr_curve.txt" \
    -s "../figures/HG00514-pr.png"
python plot_pr_curve.py -t "HG00733 Precision-Recall Curves" \
    -l "CNN" "Duphold DHFFC" \
    -f "~/data/HG00733/VCF/CNN_8_6/pr_curve/pr_curve.txt" \
       "~/data/HG00733/VCF/duphold/pr_curve/pr_curve.txt" \
    -s "../figures/HG00733-pr.png"
python plot_pr_curve.py -t "NA19240 Precision-Recall Curves" \
    -l "CNN" "Duphold DHFFC" \
    -f "~/data/NA19240/VCF/CNN_8_6/pr_curve/pr_curve.txt" \
       "~/data/NA19240/VCF/duphold/pr_curve/pr_curve.txt" \
    -s "../figures/NA19240-pr.png"
