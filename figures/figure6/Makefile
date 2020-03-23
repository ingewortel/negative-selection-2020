.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=30
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
	@cd latex && pdflatex figure6.tex> texlog.log


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

progress/panelB-train-peps-R-r% : progress/panelE-test-all | data/fixtest-repertoires/peptides/trainsets progress
	@echo "...Making train sets.." && for sim in $$(seq 1 $(NSIM) ) ; do \
		bash ../shared-scripts/repertoires/foreign_unseen.sh data/testsets/all-sim$$sim.txt $(trainset) data/fixtest-repertoires/peptides/trainsets/tmp-sim$$sim.txt && \
		for n in $(pep_selfsizes) ; do \
			bash ../shared-scripts/repertoires/train-epitopes.sh data/fixtest-repertoires/peptides/trainsets/tmp-sim$$sim.txt $$n > data/fixtest-repertoires/peptides/trainsets/traineps-r$*-R-n$$n-sim$$sim.txt ;\
		done ;\
	done && touch $@

progress/panelB-train-peps-O : $(foreach r,$(rvec), progress/panelB-train-peps-O-r$(r) )
	@touch $@
	
progress/panelB-train-peps-R : $(foreach r,$(rvec), progress/panelB-train-peps-R-r$(r) )
	@touch $@


# Step 4: Generate repertoires using these generated test sets.
data/mkfiles/panelB-repertoires-O.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	progress/panelB-train-peps-O | data/fixtest-repertoires/peptides/contiguous \
	data/fixtest-repertoires/peptides/trainsets data/mkfiles
	@bash $< data/peptide-sets/greedy "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/fixtest-repertoires/peptides/contiguous -T "O" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

data/mkfiles/panelB-repertoires-R.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	progress/panelB-train-peps-R | data/fixtest-repertoires/peptides/contiguous \
	data/fixtest-repertoires/peptides/trainsets data/mkfiles
	@bash $< $(trainset) "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/fixtest-repertoires/peptides/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@


progress/panelB-repertoires-O :  data/mkfiles/panelB-repertoires-O.mk progress/panelB-train-peps-O
	@echo "...Generating repertoires for Figure6B" && $(MAKE) -f $< && touch $@

progress/panelB-repertoires-R :  data/mkfiles/panelB-repertoires-R.mk progress/panelB-train-peps-R
	@echo "...Generating repertoires for Figure6B" && $(MAKE) -f $< && touch $@

# Step 5: Analyze repertoires using generated test sets
data/mkfiles/frequencies-panelB-O.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh \
	progress/panelB-repertoires-O | data/mkfiles
	@bash $< data/fixtest-repertoires/peptides "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		$(NSIM) "O" -n 6 -m "contiguous" > $@
		
data/mkfiles/frequencies-panelB-R.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh \
	progress/panelB-repertoires-R | data/mkfiles
	@bash $< data/fixtest-repertoires/peptides "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		$(NSIM) "R" -n 6 -m "contiguous" > $@

# Use the makefile above to generate the actual files with precursor frequencies.
# They will be stored in data/fixtest-pfout/contiguous/O/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/fixtest-pfout/contiguous/O/output-r[t value]-peps-[pathogen]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/fixtest-pfout/contiguous/O/output-r[t value]-peps-[pathogen]-n[ntrain]-r[t value].txt
progress/panelB-frequencies-O : data/mkfiles/frequencies-panelB-O.mk \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection progress/panelB-repertoires-O
	@echo "...Analyzing precursor frequencies for Figure 6B..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies-R : data/mkfiles/frequencies-panelB-R.mk \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection progress/panelB-repertoires-R
	@echo "...Analyzing precursor frequencies for Figure 6B..." && $(MAKE) -f $< && touch $@


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


# ==================== PANEL C ================================================

# Step 1: Build graph to analyze exchangeability
data/exch-all.txt : ../shared-scripts/analysis/graph-exch.py $(trainset) ../input-data/peptides/hiv-6mers.txt | data
	@echo "Building graph to get exchangeability" && python3 $^ 4 | sort > $@

# Step 2: Plot figure.
plots/F6panelC.pdf : ../shared-scripts/plotting/plot-exchangeability-compare.R \
	data/exch-all.txt progress/greedy ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 6C" &&  Rscript $< data/exch-all.txt $(trainset) data/peptide-sets/greedy-r4.txt $@

# ==================== PANEL D ================================================

# Self-hiv discrimination after negative selection on different sizes of trainsets
# and with t = 3 or t = 4 - random or optimal set. The data for the random sets was
# already generated in the previous figure.

# *** STEP 1: Make the repertoires just like in figure 5, but now with optimal trainsets.
# We first make a Makefile in data/mkfiles/panelD-repertoires.mk
# for generating the repertoires. This ensures that we can then
# make all required repertoires in parallel. 
data/mkfiles/panelD-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	progress/greedy| \
	data/repertoires/peptides/contiguous data/repertoires/peptides/trainsets data/mkfiles
	@bash $< data/peptide-sets/greedy "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/peptides/contiguous -T "O" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

# Ensure that the repertoires for panelB are made first, as these also are used for panelC.
progress/panelD-repertoires : data/mkfiles/panelD-repertoires.mk | progress
	@echo "...Generating repertoires for Figure6D..." && $(MAKE) -f $< && touch $@

# *** STEP 2: After generating the repertoires, we measure TCRs/test peptide.
# We first generate makefiles to count TCRs per test peptide for all simulations,
# for all test peptides of the given pathogen. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test peptides of size NTEST (which is set on top of this page). 
# For this figure, we use r=3/r=4 but take all values of ntrain.
data/mkfiles/panelD-frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelD-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "O" \
		-n 6 -m "contiguous" -u > $@

data/mkfiles/panelD-frequencies-hiv.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/hiv-6mers.txt progress/panelD-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/hiv-6mers.txt "hiv" $(MODELDIR) $(NTEST) $(NSIM) "O" \
		-n 6 -m "contiguous" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/panelD-frequencies-O-% : data/mkfiles/panelD-frequencies-%.mk progress/panelD-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection \
	../shared-scripts/negselmodel/src/countpaths
	@echo "...Analyzing precursor frequencies for Figure 6D $* ..." && $(MAKE) -f $< && touch $@

# The frequencies for the random trainsets were computed in figure 5, copy them here.
progress/panelD-frequencies-R : | progress data/pfout/contiguous/R
	@echo "...Getting precursor frequencies from Figure 5 for random trainsets" &&\
	cd ../figure5 && $(MAKE) progress/panelD-frequencies && \
	cp data/pfout/contiguous/R/output* ../figure6/data/pfout/contiguous/R/ &&\
	touch $@

progress/panelD-frequencies : progress/panelD-frequencies-O-hiv progress/panelD-frequencies-O-self \
	progress/panelD-frequencies-R
	@touch $@
	
# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/panelD-analysis-%.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelD-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/$* "$(pep_selfsizes)" "$(rvec)" "hiv" "$(NSIM)" \
		"data/panelD-precursors-$*" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/panelD-precursors-%-calc.txt : data/mkfiles/panelD-analysis-%.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelD-frequencies
	@echo "...Analyzing repertoires for Figure 5D" && $(MAKE) -f $< $@

# Combine the results with random (R) and optimal (O) sets into a single file:
data/panelD-precursors-greedy-calc.txt : data/panelD-precursors-O-calc.txt data/panelD-precursors-R-calc.txt
	cat $< | head -1 > $@ &&\
	cat $^ | grep -v type >> $@


# *** STEP 4: Make the plot.
# For this figure only t = 4.
plots/F6panelD.pdf : ../shared-scripts/plotting/plot-top.R \
	data/panelD-precursors-greedy-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/panelD-precursors-greedy-calc.txt 4 "hiv" "sem" "peps" $@


# ==================== PANEL E ================================================

# We only need to make the plot, since we use the data generated for panel C.
plots/F6panelE.pdf : ../shared-scripts/plotting/plot-exch-vs-foreign.R \
	data/exch-all.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/exch-all.txt plots/F6panelE


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

data/pfout/contiguous/R:
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
	@rm -rf data && rm -rf plots && rm -rf progress && rm -f figure6.pdf
