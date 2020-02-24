# This folder contains code to generate the language strings from the downloaded bible
# texts in the 'samples/' folder. We include the strings used for this manuscript in the
# 'data/' folder for more exact replicability, but these can be overwritten by new, 
# randomly generated strings if desired: type 'make clean' followed by 'make' to generate
# new sets.


.SECONDARY:
.DELETE_ON_ERROR:

# A small number needed to compute the biased sampling probability.
epsilon=0.01

# Names of the languages in shorthand and full
languages=en en2 hi la me pd ta xh en_t
languages_long=english english hiligaynon latin middle-english plautdietsch tagalog xhosa english-train english-train
language_i=$(shell seq 1 $(words $(languages)))

all: $(foreach l,$(languages), data/$(l)-unseen.txt)

# ------------ LANGUAGE DATA
# This script doesn't really do anything except show that the data/*-unseen.txt files
# indeed do not overlap with data/en_t.txt anymore (type 'make check' to see this).
check : scripts/check-unseen.sh
	bash $< "$(languages)"

# Generate files of unseen foreign by comparing each file of foreign chunks against the
# file of training set chunks.
data/%-unseen.txt : scripts/foreign_unseen.sh data/%.txt data/en_t.txt | data
	bash $< data/en_t.txt data/$*.txt data/$*-unseen.txt

# This code takes the text from samples/[language-file].txt, and chops it into chunks
# of the correct length (in this case, 6). 
define CHUNKS
data/$(1).txt : scripts/chunkify.py samples/$(2).txt | data
	python $$< samples/$(2).txt | awk 'length==6' > $$@
endef
$(foreach i,$(language_i),$(eval $(call CHUNKS,$(word $(i),$(languages)),$(word $(i),$(languages_long)))))

data : 
	mkdir -p $@

clean:
	rm -rf data

