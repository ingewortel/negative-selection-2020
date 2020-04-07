#!/bin/bash

# This script loops over values of r and trainset sizes and computes precursor frequencies
# for test epitopes. Multiple simulations are stored in one outputfile. 
# Repertoirefiles must already be generated beforehand and present in the repertoire directory.
#
# INPUT:
# Arg 1:	directory where repertoires are located
# Arg 2: 	"selfsizes" - a vector of training set sizes to use.
# Arg 3:	"rvec" - vector of r values to loop over
# Arg 4:	directory where model code is located
# Arg 5: 	Number of simulations
# Arg 6:	Type of training set selection (eg "R","O","B1",...).
# Arg 7:	Directory of the test epitope files.
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
# 
# OUTPUT: Files in data/fixtest-pfout/MATCHTYPE/TYPE/output-r[r]-[lang/peps]-[name]-n[ntrain]-r[r].txt
# with simulations underneath each other. 


# Arguments
REPERTOIREDIR=$1
SELFSIZES=$2
RVEC=$3
MODELDIR=$4
NSIM=$5
TYPE=$6

# Defaults:
MOTIFLENGTH=6
MATCHTYPE="pattern-exact"
L=""
UNSEEN=0

# Options and flags:
OPTIND=7
while getopts ":n:m:l:" opt; do
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
		\?)
			echo "Invalid option: -$OPTARG" >&2
			;;
	esac
done

# ---------------------------------------------------------------------
# CODE:
# ---------------------------------------------------------------------

echo "all : "
DIR=data/fixtest-pfout/$MATCHTYPE/$TYPE

nr=$( echo ${RVEC[@]} | wc -w )
nn=$( echo ${SELFSIZES[@]} | wc -w )

repnum=$( echo "$nr * $nn * $NSIM" | bc )
simfiles=""
outfiles=""

for R in ${RVEC[@]} ; do

	for ntrain in ${SELFSIZES[@]} ; do


		# Output target filename:
		if [ "$L" == "-lang" ] ; then
			OUTFILE=$DIR/output-r$R-lang-all-n$ntrain
		else
			OUTFILE=$DIR/output-r$R-peps-all-n$ntrain
		fi
	
		simfiles1=""
		for sim in $(seq 1 $NSIM) ; do

			# write intermediate targets to makefile
			echo "$OUTFILE-sim$sim.txt : $EPFILE"
			echo -e "\t@"bash ../shared-scripts/analysis/precursor-counts-persim-fixtest.sh $REPERTOIREDIR $ntrain $R data/testsets/all-sim$sim.txt $MODELDIR $sim $TYPE \$\@ -n $MOTIFLENGTH -m $MATCHTYPE -l $L 
		
			simfiles=$simfiles" "$OUTFILE-sim$sim.txt
			simfiles1=$simfiles1" "$OUTFILE-sim$sim.txt
		done
		
		# Write final target to makefile
		echo "$OUTFILE-r$R.txt : $simfiles1"
		echo -e "\t@"cat \$\^ " > \$@ &&"rm \$\^
		echo ""

		outfiles=$outfiles" "$OUTFILE-r$R.txt

	done

done

#echo "final : $simfiles"
#echo -e "\t@"rm $simfiles
echo ""
echo "all : $outfiles"

echo "message : "
echo -e "\t@"echo "Computing precursor counts for $EPITOPEID for $repnum files, recipe: ../scripts/analysis/precursor-counts-persim-fixtest.sh $REPERTOIREDIR [ntrain] [R] data/testsets/all-sim[sim].txt $MODELDIR [sim] $TYPE \$\@ -n $MOTIFLENGTH -m $MATCHTYPE -l $L"








