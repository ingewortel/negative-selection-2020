#!/bin/bash

selffile=$1
sims=$2

mkdir -p data/peps-concordance/tmp


# temporary files
tmself=data/peps-concordance/tmp/tmself-$RANDOM-sim$sim.txt
trap "rm $tmself" EXIT


# Step 1: randomly pick 100 self/nonself test set

for sim in $(seq 1 $sims) ; do
	trap "exit 1" SIGINT SIGTERM

	outfile=data/peps-concordance/out-sim$sim.txt

	# Shuffle the self strings
	cat $selffile | sort | uniq | perl -MList::Util -e 'print List::Util::shuffle <>' > $tmself

	# First 100 strings: "foreign" test
	cat $tmself | head -n 500 | awk '{print "foreign", $1}' > $outfile

	# Second 100: "self" test set
	cat $tmself | head -n 1000 | tail -n 500 | awk '{print "self", $1}' >> $outfile

done
