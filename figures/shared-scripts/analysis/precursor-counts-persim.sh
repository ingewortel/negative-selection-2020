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
EPFILE=$4
EPITOPEID=$5
MODELDIR=$6
NTEST=$7
sim=$8
TYPE=$9
OUTFILE=${10}

# Defaults:
MOTIFLENGTH=6
MATCHTYPE="pattern-exact"
L=""
UNSEEN=0

# Options and flags:
OPTIND=11
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


FILEID=$RANDOM-pf-$EPITOPEID-$MATCHTYPE-$sim-$ntrain-$R

mkdir -p data/pfout/$MATCHTYPE/$TYPE
tmtesteps=data/tmp/testeps-$FILEID.txt
tmrepertoire=data/tmp/rep-$FILEID.fst

mkdir -p data/tmp
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

# Create a test set
# If option -u, UNSEEN=1 and test epitopes should be unseen:
if [ $UNSEEN -eq 1 ] ; then
	comm -13 <( sort $TRAINFILE ) <( sort $EPFILE ) | perl -MList::Util -e 'print List::Util::shuffle <>' | head -$NTEST > $tmtesteps

# If not, UNSEEN=0 and test epitopes are selected randomly:		
else
	cat $EPFILE | perl -MList::Util -e 'print List::Util::shuffle <>' | head -$NTEST > $tmtesteps
fi

# Compute precursor frequencies and store in OUTFILE
bash ../shared-scripts/analysis/precursor-frequencies-c.sh $tmtesteps $tmrepertoire $MODELDIR -n $MOTIFLENGTH -r $R -m $MATCHTYPE -l $L | paste $tmtesteps - | awk -v ID="$EPITOPEID" '{print ID, $0}' > $OUTFILE





