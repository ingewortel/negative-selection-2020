.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=5
NTEST=50
trainset=../input-data/peptides/self-6mers.txt
pathogens=hiv hepb hepc hcmv ebola zika mal lis vac
rvec=3 4


all : figure6.pdf

# ==================== FIGURE 6 ================================================
# This code combines all panels into a single pdf for figure 5. 
figure6.pdf : latex/figure6.pdf
	@cp $< $@

latex/figure6.pdf : latex/figure6.tex ../shared-scripts/latex/figure-preamble.tex \
	cartoon-panels simulation-panels
	@cd latex && pdflatex figure6.tex#> texlog.log


# These are the simulation plots to generate
panels=panelB panelC panelD panelE
simulation-panels :$(foreach p, $(panels), plots/F6$(p).pdf)


# ==================== PANEL A ================================================
# This figure does not have to be generated.
cartoon-panels : cartoons/clusters.pdf




# ==================== PANEL B ================================================
# For this panel, we first compute an "optimal" set of self-peptides that can
# delete as many different self-reactive TCRs as possible, and then repeat the
# analysis from Figure 5E (although here we plot the % of TCRs removed rather than
# the % of TCRs surviving selection).

# Step 1: Find the optimal set of self-peptides.
data/peptide-sets/greedy-r%.txt : ../shared-scripts/trainset-sampling/greedytol.R \
	../shared-scripts/trainset-sampling/problem-functions.R $(trainset) | data/peptide-sets
	@echo "...Computing optimal self peptides for t = $*" && \
	Rscript $< $(trainset) $* $$(cat $(trainset) | wc -l) > $@

progress/greedy : $(foreach g,$(rvec),data/peptide-sets/greedy-r$(g).txt)  | progress

# Step 2: Use the fixed test sets from Figure 5E.
progress/panelB-testsets : | data/testsets progress
	@echo "...Copying test sets from Figure 5E" &&\
	cd ../figure5 && $(MAKE) progress/panelE-test-all &&\
	cp data/testsets/all-* ../figure6/data/testsets/ && touch $@

# Step 3: Make trainsets of peptides not in the test set
progress/panelB-train-peps-O-r% : ../shared-scripts/trainset-sampling/greedy-train.r progress/greedy \
	progress/panelB-testsets ../shared-scripts/repertoires/train-epitopes.sh | \
	data/fixtest-repertoires/peptides/trainsets
	@echo "...Generating training sets" && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		Rscript ../shared-scripts/trainset-sampling/greedy-train.r data/peptide-sets/greedy-r$*.txt data/testsets/all-sim$$sim.txt data/fixtest-repertoires/peptides/trainsets/gtmp-sim$$sim-r$*.txt &&\
		for n in $(pep_selfsizes) ; do \
			bash ../shared-scripts/repertoires/train-epitopes.sh \
				data/fixtest-repertoires/peptides/trainsets/gtmp-sim$$sim $$n -T O -r \
				$* > data/fixtest-repertoires/peptides/trainsets/traineps-r$*-O-n$$n-sim$$sim.txt ;\
		done ;\
	done && touch $@

progress/panelB-train-peps-O : $(foreach r,$(rvec), progress/panelB-train-peps-O-r$(r) )
	@touch $@


# Step 4: Generate repertoires using these generated test sets.
data/mkfiles/panelB-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	progress/panelB-train-peps-O | data/fixtest-repertoires/peptides/contiguous \
	data/fixtest-repertoires/peptides/trainsets data/mkfiles
	@bash $< data/peptide-sets/greedy "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/fixtest-repertoires/peptides/contiguous -T "O" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

progress/panelB-repertoires :  data/mkfiles/panelB-repertoires.mk progress/panelB-train-peps-O
	@echo "...Generating repertoires for Figure6B" && $(MAKE) -f $< && touch $@


# Step 5: Analyze repertoires using generated test sets
data/mkfiles/frequencies-panelB.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh \
	progress/panelB-repertoires | data/mkfiles
	@bash $< data/fixtest-repertoires/peptides "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		$(NSIM) "O" -n 6 -m "contiguous" > $@

# Use the makefile above to generate the actual files with precursor frequencies.
# They will be stored in data/fixtest-pfout/contiguous/O/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/fixtest-pfout/contiguous/O/output-r[t value]-peps-[pathogen]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/fixtest-pfout/contiguous/O/output-r[t value]-peps-[pathogen]-n[ntrain]-r[t value].txt
progress/panelB-frequencies-O : data/mkfiles/frequencies-panelB.mk \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection progress/panelB-repertoires
	@echo "...Analyzing precursor frequencies for Figure 6B..." && $(MAKE) -f $< && touch $@


# The frequencies for random selection have already been computed in Figure 5, so copy
# them here.
progress/panelB-frequencies-R : | progress data/fixtest-pfout/contiguous/R
	@echo "...Getting precursor frequencies from Figure 5 for random trainsets" &&\
	cd ../figure5 && $(MAKE) progress/panelE-frequencies && \
	cp data/fixtest-pfout/contiguous/R/* ../figure6/data/fixtest-pfout/contiguous/R/ &&\
	touch $@

progress/panelB-frequencies : progress/panelB-frequencies-R progress/panelB-frequencies-O
	@touch $@

# Step 6: TCR survival rates are computed from the precursor frequencies computed in the 
# previous step. 
# Compute from the random sets
data/tsurv-R-r%.txt : progress/panelB-frequencies-R
	@echo "...Analyzing TCR survival for Figure 5E" && for n in $(pep_selfsizes) ; do \
		cat data/fixtest-pfout/contiguous/R/output-r$*-peps-all-n$$n-r$*.txt | awk -v n="$$n" -v r="$*" '{print $$0, n, r}' >> $@ ;\
	done && touch $@

# Compute from the optimal sets
data/tsurv-O-r%.txt : progress/panelB-frequencies-O
	@echo "...Analyzing TCR survival for Figure 5E" && for n in $(pep_selfsizes) ; do \
		cat data/fixtest-pfout/contiguous/O/output-r$*-peps-all-n$$n-r$*.txt | awk -v n="$$n" -v r="$*" '{print $$0, n, r}' >> $@ ;\
	done && touch $@

# Combine results
data/tsurv-RvsG-r%.txt : data/tsurv-O-r%.txt data/tsurv-R-r%.txt
	cat data/tsurv-R-r$*.txt | awk '$$1 ~ /^self/ {print "R", $$0}' > $@ &&\
	cat data/tsurv-O-r$*.txt | awk '$$1 ~ /^self/ {print "O", $$0}' >> $@


# Step 7: plot the result.
plots/F6panelB.pdf : ../shared-scripts/plotting/plot-tsurv-selectiontype.R \
	data/tsurv-RvsG-r4.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 6B" && Rscript $< data/tsurv-RvsG-r4.txt $(NSIM) $(NTEST) "peps" $@






# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	@cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	@mkdir -p data 
	
data/peptide-sets :
	@mkdir -p $@
	
data/peptides-fixtest/trainsets :
	@mkdir -p $@

data/fixtest-pfout/contiguous/R :
	@mkdir -p $@

data/fixtest-repertoires : 
	@mkdir -p $@

data/fixtest-repertoires/peptides/contiguous :
	@mkdir -p $@
	
data/fixtest-repertoires/peptides/trainsets :
	@mkdir -p $@
	
data/repertoires : 
	@mkdir -p $@

data/repertoires/peptides/contiguous :
	@mkdir -p $@
	
data/repertoires/peptides/trainsets :
	@mkdir -p $@
	
data/testsets :
	@mkdir -p $@
	
data/mkfiles :
	@mkdir -p $@
	
progress :
	@mkdir -p $@
	
plots :
	@mkdir -p plots

clean: latex-clean
	@rm -rf data && rm -rf plots && rm -rf progress
