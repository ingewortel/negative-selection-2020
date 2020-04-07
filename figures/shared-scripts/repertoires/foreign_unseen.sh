#!/bin/bash

#DESCRIPTION:
#create "unseen foreign" file.

TRAINFILE=$1
FOREIGNFILE=$2
OUTFILE=$3

# comm line 1: lines unique to file1, line 2 : lines unique to file2, line 3: lines in both. 
# We need lines unique to foreign file, so comm -13 to suppress the other columns

comm -13 <( sort $TRAINFILE ) <( sort $FOREIGNFILE ) > $OUTFILE
