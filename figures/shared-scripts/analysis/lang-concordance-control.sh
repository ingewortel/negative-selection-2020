#!/bin/bash

trainfile=$1
sim=$2
modeldir=$3

mkdir -p data/lang-concordance/tmp
mkdir -p data/lang-concordance/repertoires
mkdir -p data/lang-concordance/trainsets

# temporary files
tmself=data/lang-concordance/tmp/tmself-$RANDOM-sim$sim.txt
tmtest=data/lang-concordance/tmp/tmtest-$RANDOM-sim$sim.txt
trap "rm $tmself; rm $tmtest" EXIT


# Step 1: create a set of 1000 strings, randomly divide into train and self/nonself test set

	# Create a set of 1000 strings from the English training sample, and shuffle
	python ../shared-scripts/get6mers/chunkify.py $trainfile | awk '{ print $1, length}' | grep 6 | awk '{print $1}' | perl -MList::Util -e 'print List::Util::shuffle <>' > $tmself

	# First 800 strings: "train set"
	trainfile=data/lang-concordance/trainsets/train-sim$sim.txt
	cat $tmself | head -n 800 > $trainfile

	# Last 100 strings: "foreign" test
	cat $tmself | tail -n 100 | awk '{print "en_t", $1}' > $tmtest

	# remaining 100: "self" test set
	cat $tmself | tail -n 200 | head -n 100 | awk '{print "self", $1}' >> $tmtest

# Step 2: negative selection on trainset

	# Create a repertoire using the training set
	repfile=data/lang-concordance/repertoires/rep-sim$sim.fst
	bash ../shared-scripts/repertoires/repertoires.sh $repfile.gz $modeldir -s $trainfile -n 6 -r 3 -m "contiguous" -l "-lang"

# Step 3: compute precursor frequencies
	# Unpack the repertoire
	gunzip -k $repfile.gz

	# Compute precursor frequencies for the test epitopes. 
	outfile=data/lang-concordance/out-sim$sim.txt
	bash ../shared-scripts/analysis/precursor-frequencies-c.sh <( cat $tmtest | awk '{print $2}' ) $repfile $modeldir -n 6 -r 3 -m "contiguous" -l "-lang" | paste $tmtest - > $outfile


