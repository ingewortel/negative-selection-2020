#!/bin/bash

# DESCRIPTION:
# Create "unseen foreign" file. It takes as input a file of training chunks and a file of
# foreign chunks. It then ensures that none of the output foreign chunks also occur in the
# training set file, so that they are 'unseen' no matter which sample from the training
# set we draw.

TRAINFILE=$1
FOREIGNFILE=$2
OUTFILE=$3

# comm line 1: lines unique to file1, line 2 : lines unique to file2, line 3: lines in both. 
# We need lines unique to foreign file, so comm -13 to suppress the other columns

comm -13 <( sort $TRAINFILE ) <( sort $FOREIGNFILE ) > $OUTFILE
