#CURRENT PATH: 2020-negative-selection-repo/figures/simulation-data/repertoires

MODELDIR=../../shared-scripts/negselmodel/src

# This file contains some settings in variables used below (eg. 'nsim' is defined here)
include ../../settings.env
SIMS=1 $(nsim)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: check-lang-contiguous data/new-lang


data : 
	mkdir -p data

clean:
	rm -rf data && rm -rf languages


# --------------------------------------------------------------------
#				LANGUAGES
# --------------------------------------------------------------------

# print progress to console
check-lang-contiguous : 
	bash scripts/repertoire-status.sh "." "$(lang_selfsizes)" "$(rvec)" "R" "languages" "contiguous" $(nsim)

# Make the repertoires using the makefile 'data/rand-contiguous-lang.mk', which is
# generated below. This Makefile will first 
# generate a training set for each simulation and each training set size using the 
# script 'scripts/train-epitopes.sh'. For the languages, we use random training set 
# selection so these sequences are sampled randomly from the file 
# '../languages/data/en_t.txt'. The training sets are then used to produce the 
# corresponding repertoire using the code in $(MODELDIR).
data/rand-contiguous-lang : data/rand-contiguous-lang.mk data/add-lang-contiguous \
	scripts/train-epitopes.sh scripts/repertoires.sh ../languages/data/en_t.txt \
	../../shared-scripts/analysis/print-progress.sh $(MODELDIR)/contiguous-fa-lang \
	$(MODELDIR)/makerep-contiguous-fa-lang
	$(MAKE) -f $< message && \
	$(MAKE) -f $< &&\
	touch $@ && touch data/add-lang-contiguous


# Generate a makefile to auto-generate all the repertoires. See data/rand-contiguous-lang
# above for a description.
data/rand-contiguous-lang.mk : scripts/loop-repertoires-makeout.sh data/add-lang-contiguous data/total-lang.txt
	mkdir -p languages/contiguous && \
	mkdir -p languages/trainsets && \
    	bash $< ../languages/data/en_t.txt "$(lang_selfsizes)" "$(rvec)" $(MODELDIR) languages/contiguous data/total-lang.txt -T "R" -S "$(SIMS)" -n 6 -m "contiguous" -l "-lang" > $@

# This is just a helper to tell the program that new repertoires should be added.
data/add-lang-contiguous : | data
	touch $@

# Compute the total number of repertoires that will be generated at the current settings
# (number of training set sizes, number of sims, number of values of the model parameter t)	
data/total-lang.txt : | data
	nn=$$( echo $(lang_selfsizes) | wc -w ) &&\
	nr=$$( echo $(rvec) | wc -w) &&\
	nreps=$$( echo "$(nsim)*$$nn*$$nr" | bc) &&\
	echo $$nreps > $@ && touch data/add-lang-contiguous

data/new-lang : | data
	rm -rf languages/trainsets ;\
	rm -rf languages/contiguous ;\
	rm -rf data/rand-contiguous-lang