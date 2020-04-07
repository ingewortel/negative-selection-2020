#!/bin/bash

# This script takes a repertoire .fst file and a set of test epitopes, and
# calculates precursor frequencies for each epitope.
#
# INPUT:
# The first argument after the command contains the name of the file with
# test epitopes. This file should contain the epitope sequence in the first
# column. The second argument
# contains the name of the repertoire .fst file. The third is the directory
# where model code is located. 
# 
# Other arguments are:
#	-	-n	The length of epitopes
#			(Default 6 AAs)
#	-	-r	The value of r in the r-contiguous/r-pattern
#			(Default 6)
#	-	-m	The type of matching to use: "contiguous", 
#			"pattern", "pattern-exact", etc.
#			(The script calls ./MATCHTYPE-fa, so argument
#			for MATCHTYPE can be any string for which there
#			exists this code. Default is "contiguous".
#	-	-l	With argument "-lang" for using the language alphabet 
#			(abcdefghijklmnopqrstuvwxyz_) instead of the amino
#			acids (default).
# 
# OUTPUT: STDOUT with the following columns:
#	1	type of test epitope (e.g. "self" or "hiv", as specified in
#		the input epitope file. 
#	2	The test epitope sequence (as in the input epitope file)
#	3	The total number of survivors in the simulated repertoire
#	4	Precursor count for current value of r


# ---------------------------------------------------------------------
# INPUT ARGUMENTS
# ---------------------------------------------------------------------

# Arguments:
TESTEPFILE=$1
REPERTOIREFILE=$2
MODELDIR=$3


# Some defaults:
MOTIFLENGTH=6
R=6
MATCHTYPE="pattern-exact"
L=""

# Flags and arguments:
OPTIND=4
while getopts ":n:r:m:l:" opt; do
	case $opt in
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



# ---------------------------------------------------------------------
# CODE
# ---------------------------------------------------------------------

# Total number of T cells in the current repertoire: 

totalsurvivors=$(cat $REPERTOIREFILE | fstprint | $MODELDIR/countpaths)

# Compute precursor frequencies for all foreign epitopes.

$MODELDIR/$MATCHTYPE-negative-selection$L $TESTEPFILE $MOTIFLENGTH $R \
	< $REPERTOIREFILE | awk -v t=$totalsurvivors '{print t,$1,""}'


#| \
#	awk -v t=$totalsurvivors '{print t,$1,""}'

#for lineno in $(seq 1 $(wc -l < $TESTEPFILE)) ; do

	# Store current line as 'peptide'.  
#	peptide=$(cat $TESTEPFILE | awk "NR==$lineno"'{print $2}')

	# count total number of survivors and add to file, followed by space
#	echo -n $totalsurvivors" "

	# count matches for given R and echo it to the file.
	#survivors=$(echo $peptide | $MODELDIR/$MATCHTYPE-fa$L $MOTIFLENGTH $R | fstcompile --acceptor | fstintersect $REPERTOIREFILE - | fstprint | $MODELDIR/countpaths)

	# echo the number of survivors to the file.
#	echo -n $survivors" "

#	echo ""
#done
