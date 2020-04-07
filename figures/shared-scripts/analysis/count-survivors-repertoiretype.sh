REPERTOIREDIR=$1
MODELDIR=$2
TYPES=$3
RVEC=$4
SIZES=$5
NSIM=$6



for R in ${RVEC[@]} ; do

	for ntrain in ${SIZES[@]} ; do

		for tt in ${TYPES[@]} ; do

			for s in $(seq 1 $NSIM) ; do

				repfile=$REPERTOIREDIR/repertoire-r$R-$tt-n$ntrain-sim$s.fst.gz
				surv=$(gunzip -c $repfile | fstprint | $MODELDIR/countpaths)
				echo $R $ntrain $tt $s $surv

			done
		done
	done
done
