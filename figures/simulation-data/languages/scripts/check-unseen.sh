
languages=$1

for l in ${languages[@]} ; do
	echo $l :
	echo "       " $( cat data/$l.txt | wc -l ) lines
	echo "       " Overlap with en_t.txt : $( comm -12 <(sort data/en_t.txt) <(sort data/$l.txt) | wc -l ) lines
	echo "       " Unseen : $( cat data/$l-unseen.txt | wc -l ) lines
	echo "       " Overlap with en_t.txt : $( comm -12 <(sort data/en_t.txt) <(sort data/$l-unseen.txt) | wc -l ) lines
done
