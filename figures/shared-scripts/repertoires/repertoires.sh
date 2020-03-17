#!/bin/bash

# This script performs negative selection with a model of choice and
# stores the automaton describing the post-selection repertoire in .fst.gz format.
#
# INPUT:
# The first argument after the command is the name of the .fst.gz file to
# produce. The second argument contains folder where the model programs
# are located. All other options are set with flags as indicated below.
# 
# OPTIONS:
# 	-	-s 	The name of the file with self epitopes to train on
#			(Default is empty training set)
#	-	-n	The length of epitopes
#			(Default 6 AAs)
#	-	-r	The value of r in the r-contiguous/r-pattern model.
#			This parameter is called 't' in the paper. (Default 6)
#	-	-m	The type of matching to use: "contiguous", 
#			"pattern", "pattern-exact", etc.
#			(The script calls ./MATCHTYPE-fa, so argument
#			for MATCHTYPE can be any string for which there
#			exists this code. Default is "contiguous".
#			We use 'contiguous' throughout the paper, but have obtained
#			similar results using other rules (not published).
#	-	-l	With argument "-lang" for using the language alphabet 
#			(abcdefghijklmnopqrstuvwxyz_) instead of the amino
#			acids (default). 
# 
# OUTPUT: .fst.gz file with the surviving T cell repertoire. 



# ---------------------------------------------------------------------
# INPUT ARGUMENTS
# ---------------------------------------------------------------------

# Preparation: make a directory for temporary files

mkdir -p data/tmp

# The main input argument is the set of test epitopes

OUTFILENAME=$1
MODELDIR=$2


# Some defaults:

echo -n > data/tmp/emptyself.txt
TRAINSELF=data/tmp/emptyself.txt
MOTIFLENGTH=6
R=6
#FILEID=$(echo $(date | sed 's/ /_/g' | sed 's/://g')$RANDOM)
FILEID=$(echo $OUTFILENAME | sed 's@\/@@g' | sed 's/fst.gz/fst/g')-$RANDOM
MATCHTYPE="pattern-exact"
L=""

OPTIND=3
while getopts ":s:n:r:m:l:" opt; do
	case $opt in
		s)
			TRAINSELF=$OPTARG >&2
			;;
		n)
			MOTIFLENGTH=$OPTARG >&2
			;;
		r)
			R=$OPTARG >&2
			;;
		m)
			MATCHTYPE=$OPTARG >&2
			;;
		l)
			L=$OPTARG >&2
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			;;
	esac
done



# temporary files

echo -n > data/tmp/all-$FILEID #$FILEID-$R.fst
echo -n > data/tmp/self-$FILEID #$FILEID-$R.fst
echo -n > data/tmp/neg-$FILEID #$FILEID-$R.fst
TEMP_FILES="data/tmp/all-$FILEID data/tmp/self-$FILEID data/tmp/neg-$FILEID data/tmp/neg-$FILEID.gz"

function clean_up {

	# Perform program exit housekeeping
	rm $TEMP_FILES

}

trap "clean_up" EXIT


# ---------------------------------------------------------------------
# AUTOMATON GENERATION
# ---------------------------------------------------------------------
# Step 1 : Make an automaton for positively selected repertoire
# 
# Based on the training epitopes, generate an automaton to describe a
# "positively selected repertoire": all TCRs that recognize any of these
# sequences. 
# (It can make sense to put this first step into a different script or put 
# it into the makefile, because this can take long.)


cat $TRAINSELF | $MODELDIR/$MATCHTYPE-fa$L $MOTIFLENGTH $R | fstcompile --acceptor > data/tmp/self-$FILEID



# Step 2 : Make an automaton that contains all possible TCRs.

$MODELDIR/makerep-$MATCHTYPE-fa$L $MOTIFLENGTH $R  | fstcompile --acceptor > data/tmp/all-$FILEID




# ---------------------------------------------------------------------
# GET REPERTOIRE 
# ---------------------------------------------------------------------
# Step 3: Subtract positively selected from all TCR to obtain negatively selected TCR.

fstdifference data/tmp/all-$FILEID data/tmp/self-$FILEID | fstminimize > data/tmp/neg-$FILEID

# Now zip the repertoire file
gzip -k data/tmp/neg-$FILEID

# Move/rename to outfilename.
cp data/tmp/neg-$FILEID.gz $OUTFILENAME



