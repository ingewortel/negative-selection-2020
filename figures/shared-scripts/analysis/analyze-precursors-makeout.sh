# This script loops over values of t, trainset sizes, and languages. 
# It creates a makefile to analyze each file with precursor frequencies, and
# concatenates those to two large dataframes: $FILE-calc.txt and $FILE-sim1.txt
# 
#
# INPUT:
# Arg 1:	directory where precursor counts are located
# Arg 2: 	"selfsizes" - a vector of training set sizes to use.
# Arg 3:	"rvec" - vector of r values to loop over. This is the model parameter t.
# Arg 4:	"languages" - a vector of the languages to analyze.
# Arg 5:	"nsim" - the number of simulations.
# Arg 6:	"file" - prefix for the name of the output files to generate.
# 
# 
# Other arguments are:
#	-	-l	With argument "-lang" for using the language alphabet 
#			(abcdefghijklmnopqrstuvwxyz_) instead of the amino
#			acids (default). 
# 
# OUTPUT: Makefile that generates files : data/$FILE-calc.txt and data/$FILE-sim1.txt

# Arguments
DIR=$1
SELFSIZES=$2
RVEC=$3
LANGUAGES=$4
NSIM=$5
FILE=$6


# This code allows the script to take options with flags.
# Here, the option -l lang allows the user to specify that the language version
# of the model should be used (not the peptide version).
OPTIND=7
while getopts ":l:" opt; do
	case $opt in
		l)
			L=$OPTARG >&2
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			;;
	esac
done

# Handle input arguments
if [ "$L" == "-lang" ] ; then
	STRINGID="lang"
else
	STRINGID="peps"
fi


# Generate the directory for analysis if it is not already there
mkdir -p $DIR/analysis

# Count how many repertoires will be analyzed
nr=$( echo ${RVEC[@]} | wc -w )
nn=$( echo ${SELFSIZES[@]} | wc -w )
nl=$( echo ${LANGUAGES[@]} | wc -w )
repnum=$( echo "$nr * $nn * $nl * 2" | bc )

# These two variables track the files that are generated.
filelist1=""
filelist2=""

echo "all : $FILE-calc.txt $FILE-sim1.txt"

# For each value of t, ntrain, and language, a single file is generated first.
for R in ${RVEC[@]} ; do
	for ntrain in ${SELFSIZES[@]} ; do
		for lang in ${LANGUAGES[@]} ; do

			OUTFILE=$DIR/analysis/analysis-r$R-$STRINGID-$lang-n$ntrain-calc.txt
			OUTFILE2=$DIR/analysis/analysis-r$R-$STRINGID-$lang-n$ntrain-sim1.txt
			echo "$OUTFILE : ../shared-scripts/analysis/analyze-precursors3.R"
			echo -e "\t@"Rscript ../shared-scripts/analysis/analyze-precursors3.R $R $ntrain $lang $NSIM $STRINGID $DIR

			filelist1=$filelist1" "$OUTFILE
			filelist2=$filelist2" "$OUTFILE2

		done
	done
done

# This concatenates all the small files from the previous loop into one large dataframe
# for each type of output.
echo "$FILE-calc.txt : $filelist1"
echo -e "\t@cat \$< | head -1  > \$@ && cat \$^ | grep -v ntrain  >> \$@"

echo "$FILE-sim1.txt : $filelist1"
echo -e "\t@cat $OUTFILE2 | head -1  > \$@ && cat $filelist2 | grep -v ntrain  >> \$@"
