#!/bin/bash

# This script generates training epitope sets of the correct size. 
#
# INPUT:
# First argument: 	file with training epitopes to select from
# Second argument: 	NTRAIN number of training epitopes.
# 
# OPTIONS:
# 	-	-T 	The type of epitope selection: R (RANDOM), O (OPTIMAL), 
#			or B (BIASED). Default: Random.
#			Note that if T=B, the file with training epitopes should
#			have a second column with the (biased) probabilities of
#			drawing each epitope. 
#	-	-b	Optional argument for biased training: the power to raise
#			the probability to (increasing bias strength). Default: 1. 
#			This argument is ignored if T is set to R or O.
#	-	-r	The value of r in r-contiguous model (this parameter is called
#			't' in the paper). Important when T = O.
# 
# OUTPUT: chosen training epitopes to STDOUT



# ---------------------------------------------------------------------
# PREPARATION
# ---------------------------------------------------------------------

# Input arguments:
TRAINSELF=$1
NTRAIN=$2

# Defaults:
TYPE="R"
pow=1

# OPTIONS:
OPTIND=3
while getopts ":T:b:r:" opt; do
	case $opt in
		T)
			TYPE=$OPTARG >&2
			;;
		b)	
			pow=$OPTARG >&2
			;;
		r)
			R=$OPTARG >&2
			;;
    		\?)
      			echo "Invalid option: -$OPTARG" >&2
      			;;
  	esac
done

# ---------------------------------------------------------------------
# CODE
# ---------------------------------------------------------------------

# Epitope choice dependent on type:
	
if [ $TYPE == "R" ] ; then

	# Random epitope selection: shuffle lines in trainfile and select
	# the correct number of epitopes.

	cat $TRAINSELF | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $NTRAIN

elif [ $TYPE == "O" ] ; then
	
	# Selection from the optimal set: take lines from the top of the 
	# sorted file without shuffling

	cat $TRAINSELF-r$R.txt | head -n $NTRAIN 
	
elif [ $TYPE == "B" ] ; then

	# Biased epitope selection: use Rscript scripts/biased-sample-n.R
		
	Rscript scripts/biased-sample-n.R $TRAINSELF $NTRAIN $pow

else 
	echo "Unknown argument to TYPE (-T). Please choose either R (RANDOM),
		O (OPTIMIZED) or B (BIASED)."
	exit 1

fi
	



