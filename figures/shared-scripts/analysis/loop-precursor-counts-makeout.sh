#!/bin/bash

# This script loops over values of t and trainset sizes and creates a makefile to compute precursor frequencies.
#
# INPUT:
# Arg 1:	directory where repertoires are located
# Arg 2: 	"selfsizes" - a vector of training set sizes to use.
# Arg 3:	"rvec" - vector of r values to loop over
# Arg 4: 	file with epitopes to sample test epitopes from.
#		By default this sampling is done randomly. See
#		option -u to select unseen test epitopes.
# Arg 5: 	string describing type of test epitopes (eg "self", "hiv", "xh")
# Arg 6:	directory where model code is located
# Arg 7:	Number of test epitopes to select
# Arg 8: 	Number of simulations
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
# OUTPUT: Makefile that generates files : data/pfout/MATCHTYPE/TYPE/output-r[r]-[lang/peps]-[name]-n[ntrain]-sim[sim].txt

# Arguments
REPERTOIREDIR=$1
SELFSIZES=$2
RVEC=$3
EPFILE=$4
EPITOPEID=$5
MODELDIR=$6
NTEST=$7
NSIM=$8
TYPE=$9



# Defaults:
MOTIFLENGTH=6
MATCHTYPE="pattern-exact"
L=""
UNSEEN=0

# Options and flags:
OPTIND=10
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




FILEID=$RANDOM-pf-$EPITOPEID-$MATCHTYPE

mkdir -p data/pfout/$MATCHTYPE/$TYPE
tmtesteps=data/tmp/testeps-$FILEID.fst
tmrepertoire=data/tmp/rep-$FILEID.fst

mkdir -p data/tmp
echo -n > $tmtesteps


# ---------------------------------------------------------------------
# CODE:
# ---------------------------------------------------------------------

echo "all : "
DIR=data/pfout/$MATCHTYPE/$TYPE

nr=$( echo ${RVEC[@]} | wc -w )
nn=$( echo ${SELFSIZES[@]} | wc -w )

repnum=$( echo "$nr * $nn * $NSIM" | bc )
simfiles=""
outfiles=""

for R in ${RVEC[@]} ; do

	for ntrain in ${SELFSIZES[@]} ; do


		# Output target filename:
		if [ "$L" == "-lang" ] ; then
			OUTFILE=$DIR/output-r$R-lang-$EPITOPEID-n$ntrain
		else
			OUTFILE=$DIR/output-r$R-peps-$EPITOPEID-n$ntrain
		fi
	
		simfiles1=""
		for sim in $(seq 1 $NSIM) ; do

			# write intermediate targets to makefile
			echo "$OUTFILE-sim$sim.txt : $EPFILE"
			
			if [ "$UNSEEN" == "0" ] ; then
				echo -e "\t@"bash ../shared-scripts/analysis/precursor-counts-persim.sh $REPERTOIREDIR $ntrain $R $EPFILE $EPITOPEID $MODELDIR $NTEST $sim $TYPE \$\@ -n $MOTIFLENGTH -m $MATCHTYPE -l $L 
			else
				echo -e "\t@"bash ../shared-scripts/analysis/precursor-counts-persim.sh $REPERTOIREDIR $ntrain $R $EPFILE $EPITOPEID $MODELDIR $NTEST $sim $TYPE \$\@ -n $MOTIFLENGTH -m $MATCHTYPE -u -l $L
			fi
			echo ""
			simfiles=$simfiles" "$OUTFILE-sim$sim.txt
			simfiles1=$simfiles1" "$OUTFILE-sim$sim.txt
		done

		# Write final target to makefile
		echo "$OUTFILE-r$R.txt : $simfiles1"
		echo -e "\t@"cat \$\^ " > \$@ &&"rm \$\^
		echo ""

		outfiles=$outfiles" "$OUTFILE-r$R.txt
		#echo "all : $OUTFILE-r$R.txt"

	done

done

#echo "final : $simfiles"
#echo -e "\t@"rm $simfiles
echo ""
echo "all : $outfiles"

echo "message : "
echo -e "\t@"echo "Computing precursor counts for $EPITOPEID for $repnum files, recipe: bash ../shared-scripts/analysis/precursor-counts-persim.sh $REPERTOIREDIR [ntrain] [R] $EPFILE $EPITOPEID $MODELDIR $NTEST [sim] $TYPE \$\@ -n $MOTIFLENGTH -m $MATCHTYPE -l $L"








