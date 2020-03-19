#!/bin/bash

# This script computes precursor frequencies for a single simulation.
#
# INPUT:
# Arg 1:	directory where repertoires are located
# Arg 2: 	"ntrain" - current training set size.
# Arg 3:	"R" - current value of r to use
# Arg 4: 	file with epitopes to sample test epitopes from.
#		By default this sampling is done randomly. See
#		option -u to select unseen test epitopes.
# Arg 5: 	string describing type of test epitopes (eg "self", "hiv", "xh")
# Arg 6:	directory where model code is located
# Arg 7:	Number of test epitopes to select
# Arg 8: 	current simulation number
# Arg 9:	Type of training set selection (eg "R","O","B1",...).
# 
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
#			exists this code. Default is "pattern-exact".
#	-	-l	With argument "-lang" for using the language alphabet 
#			(abcdefghijklmnopqrstuvwxyz_) instead of the amino
#			acids (default). 
#	-	-u	Use this flag if test epitopes should be "unseen" (ie not in training set).
# 
# OUTPUT: precursor counts to OUTFILE.

# Arguments
REPERTOIREDIR=$1
ntrain=$2
R=$3
TESTFILE=$4
MODELDIR=$5
sim=$6
TYPE=$7
OUTFILE=$8

# Defaults:
MOTIFLENGTH=6
MATCHTYPE="pattern-exact"
L=""
UNSEEN=0

# Options and flags:
OPTIND=9
while getopts ":n:m:l:u" opt; do
	case $opt in
		n)
			MOTIFLENGTH=$OPTARG >&2
			;;
		m)
			MATCHTYPE=$OPTARG >&2
			;;
		l)
			L=$OPTARG >&2
			;;
		u)
			UNSEEN=1 >&2
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			;;
	esac
done


FILEID=$RANDOM-pf-$MATCHTYPE-$sim-$ntrain-$R

mkdir -p data/fixtest-pfout/$MATCHTYPE/$TYPE
tmtesteps=data/tmp-fixtest/testeps-$FILEID.txt
tmrepertoire=data/tmp-fixtest/rep-$FILEID.fst

mkdir -p data/tmp-fixtest
echo -n > $tmtesteps


# ---------------------------------------------------------------------
# CODE:
# ---------------------------------------------------------------------


REP=$REPERTOIREDIR/$MATCHTYPE/repertoire-r$R-$TYPE-n$ntrain-sim$sim.fst.gz
TRAINFILE=$REPERTOIREDIR/trainsets/traineps-r$R-$TYPE-n$ntrain-sim$sim.txt

# Check if the repertoire exists
if [ ! -f $REP ] ; then
	echo "Error: Cannot find repertoire $REP. Please generate repertoires before continuing."
	exit 1
fi

# copy and unpack repertoire file
cp $REP $tmrepertoire.gz
gunzip -k $tmrepertoire.gz
trap "rm $tmrepertoire; rm $tmrepertoire.gz; rm $tmtesteps" EXIT

# Use supplied test set
cat $TESTFILE | awk '{print $2}' > $tmtesteps

# Compute precursor frequencies and store in OUTFILE
bash ../shared-scripts/analysis/precursor-frequencies-c.sh $tmtesteps $tmrepertoire $MODELDIR -n $MOTIFLENGTH -r $R -m $MATCHTYPE -l $L | paste $TESTFILE - | awk -v sim="$sim" '{print $0, sim}' > $OUTFILE


