
DATASET=$1

echo REF = $(grep ref $DATASET | wc -l)
echo HET = $(grep het $DATASET | wc -l)
echo ALT = $(grep alt $DATASET | wc -l)
