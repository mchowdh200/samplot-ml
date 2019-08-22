
BED=$1
cat $1 | awk '{print $3-$2}' | (stats -b 50 -th 5000 -t "BED interval length distribution.")
