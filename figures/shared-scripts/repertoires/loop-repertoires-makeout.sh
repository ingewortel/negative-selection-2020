#!/bin/bash

# This script generates post-selection repertoires for given training set sizes,
# selection method and value of the parameter r (which we call t in the paper). 
# It first calls epitope-sets.sh to generate the epitope training sets, and then 
# uses them to generate the repertoire .fst.gz files.
#
# INPUT:
# First argument: 	TRAINSELF -- the file with self epitopes to choose from.
#			!!!	For random selection: the file with 6mers.
#			 	For optimal selection: "data/greedytol-r$$R.txt"
#				  (the right file will be chosen within the loop).
#				For biased selection: the file with 6mers + 
#				  a probability.
# Second argument: 	SELFSIZES -- vector with numbers of self epitopes to select on.
# Third argument:	Rvec -- vector of r-values to make repertoires for (this parameter
#						is called t in the paper).
# Fourth argument: 	MODELDIR -- the directory where model code is stored. 
# Fifth argument:	REPERTOIREDIR -- the directory to store the repertoire files in. 
# 
# OPTIONS:
# 	-	-T 	TYPE -- The type of epitope selection: R (RANDOM), 
#			O (OPTIMAL), or B (BIASED). Default: Random.
#	-	-b	pow -- optional argument for biased training: the power 
# 			to raise the probability to (increasing bias strength). 
#			Default: 1. This argument is ignored if T is set to R or O.
#	-	-S	SIMS -- The range of simulation numbers to use (different repertoires)
#			Default: "1 1" (simulation 1 to 1).
#	-	-n	The length of epitopes
#			(Default 6 AAs or letters)
#	-	-m	The type of matching to use: "contiguous", 
#			"pattern", "pattern-exact", etc.
#			(The script calls ./MATCHTYPE-fa, so argument
#			for MATCHTYPE can be any string for which there
#			exists this code. Default is "pattern-exact". We use 'contiguous' throughout
#			the paper but have obtained similar results using other matching rules
#			(not published).
#	-	-l	With argument "-lang" for using the language alphabet 
#			(abcdefghijklmnopqrstuvwxyz_) instead of the amino
#			acids (default). 
# 
# OUTPUT: .fst.gz files with the surviving T cell repertoires. Naming of the file:
# REPERTOIREDIR/repertoire-r[R]-[TYPE]-n[NTRAIN]-sim[sim].fst.gz
# If they don't exist yet: also makes REPERTOIREDIR/../trainsets/traineps-r[R]-[TYPE]-n[NTRAIN]-sim[sim].txt


# ---------------------------------------------------------------------
# PREPARATION
# ---------------------------------------------------------------------

TRAINSELF=$1
SELFSIZES=$2
Rvec=$3
MODELDIR=$4
REPERTOIREDIR=$5

# OPTIONS:

# Defaults
TYPE="R"
pow=1
SIMS="1 1"
MATCHTYPE="pattern-exact"
MOTIFLENGTH=6
L=""

# Flags
OPTIND=6
while getopts ":T:b:S:m:n:l:" opt; do
	case $opt in
		T)
			TYPE=$OPTARG >&2
			;;
		b)	
			pow=$OPTARG >&2
			;;
		S)
			SIMS=$OPTARG >&2
			;;
		m)
			MATCHTYPE=$OPTARG >&2
			;;
		n)
			MOTIFLENGTH=$OPTARG >&2
			;;
		l)
			L=$OPTARG >&2
			;;
    		\?)
      			echo "Invalid option: -$OPTARG" >&2
      			;;
  	esac
done

SIMS=($SIMS)

if [ $TYPE == "B" ] ; then
	TYPE2=$(echo $TYPE$pow)
else
	TYPE2=$TYPE
fi



# ---------------------------------------------------------------------
# MAKE REPERTOIRES
# ---------------------------------------------------------------------

echo "all : all-repertoires"
echo ""

allreps=""

# Loop over the simulations, selfsizes and rvalues

for sim in $(seq ${SIMS[0]} ${SIMS[1]}) ; do

	# Print progress:
	#echo .. $TYPE2 $MATCHTYPE -- Simulation $sim of ${SIMS[1]}


	for NTRAIN in $SELFSIZES ; do

		# Print progress:
		#echo ......... $TYPE2 $MATCHTYPE -- Size $NTRAIN of $SELFSIZES

		
		for R in $Rvec ; do

			trap "exit 1" SIGINT SIGTERM

			# Step 1: Check if the right training file already exists.
			# If not, make it by selecting training epitopes from TRAINSELF file.
			TRAINFILE=$REPERTOIREDIR/../trainsets/traineps-r$R-$TYPE2-n$NTRAIN-sim$sim.txt
			if [ -e $TRAINFILE ] ; then

				sleep 0
				#echo ......... $TRAINFILE exists.

			else
			
				# Trainfile does not exist yet. Make it:
				echo -e "$TRAINFILE:\n\t@"bash ../shared-scripts/repertoires/train-epitopes.sh $TRAINSELF $NTRAIN -T $TYPE -b $pow -r $R "> \$@"
				echo ""
				#sleep 0

			fi
	
			# Step 2: Create the negatively selected repertoire and save it to the
			# fst file

			OUT=$REPERTOIREDIR/repertoire-r$R-$TYPE2-n$NTRAIN-sim$sim.fst.gz
			#bash scripts/repertoires.sh $OUT $MODELDIR -s $TRAINFILE -n $MOTIFLENGTH -r $R -m $MATCHTYPE -l $L
			
			echo "$OUT : $TRAINFILE"
			echo -e "\t@"bash ../shared-scripts/repertoires/repertoires.sh \$@ $MODELDIR -s $TRAINFILE -n $MOTIFLENGTH -r $R -m $MATCHTYPE -l $L  #\&\& \\
			echo ""
			allreps="$allreps $OUT"

			#echo "all : $OUT"

			

		done

	done

done

echo "all-repertoires : $allreps"
echo ""

echo "message : "
echo -e "\t@"echo "Generating repertoires in $REPERTOIREDIR, recipe: scripts/repertoires.sh \$@ $MODELDIR -s $REPERTOIREDIR/../trainsets/traineps-r[r]-$TYPE2-n[NTRAIN]-sim[sim].txt -n $MOTIFLENGTH -r [r] -m $MATCHTYPE -l $L"


