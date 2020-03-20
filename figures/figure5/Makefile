.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=5
NTEST=50
trainset=../input-data/peptides/self-6mers.txt
pathogens=hiv hepb hepc hcmv ebola zika mal lis vac
rvec=3 4

nsim_control=100

all : figure5.pdf

# ==================== FIGURE 5 ================================================
# This code combines all panels into a single pdf for figure 5. 
figure5.pdf : latex/figure5.pdf
	@cp $< $@

latex/figure5.pdf : latex/figure5.tex ../shared-scripts/latex/figure-preamble.tex \
	cartoon-panels simulation-panels
	@cd latex && pdflatex figure5.tex#> texlog.log


# These are the simulation plots to generate
panels=panelB-lang panelB-peps-r4 panelC panelD-r3 panelD-r4 panelE-r3 panelE-r4
simulation-panels :$(foreach p, $(panels), plots/F5$(p).pdf)


# ==================== PANEL A ================================================
# This is a cartoon that is actually made for figure2, so we generate the plot
# in ../figure2/ and copy it here.
cartoon-panels : plots/F5panelA.pdf

plots/F5panelA.pdf : | plots
	@echo "Building cartoon for Figure 5A" && cd ../figure2/cartoons && $(MAKE) plots/tcr-peps.pdf && \
	cp plots/tcr-peps.pdf ../../figure5/$@

# ==================== PANEL B ================================================
# Comparing the concordances for peptides and languages. 
# 
# Step 1: concordances for languages. They already were generated in figure 4, so
# copy the data here.
# Concordances for the tested languages
data/concordances/lang-r3.txt : | data/concordances
	@echo "...Getting language concordances" && cd ../figure4 && $(MAKE) data/concordances/r3.txt &&\
	cp data/concordances/r3.txt ../figure5/$@
# Concordances for the trainset
data/concordances/lang-control-data.txt : | data/concordances
	cd ../figure4 && $(MAKE) data/concordances/control-data.txt && cp data/concordances/control-data.txt ../figure5/$@

# Step 2: plot language concordances.
plots/F5panelB-lang.pdf : ../shared-scripts/plotting/plot-concordance.R \
	data/concordances/lang-r3.txt data/concordances/lang-control-data.txt \
	@echo "...Getting language control concordances" && ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 5B for languages" Rscript $< data/concordances/lang-r3.txt data/concordances/lang-control-data.txt "lang" $@ && rm -f Rplots.pdf


# Step 3: Make a single file with 6mers from all pathogens to compute concordance on. 
# This is just the input strings used for all the figures, so this does not involve any 
# simulation of negative selection or any sampling of "training strings".
data/concordances/all-6mers.txt : $(trainset) $(foreach p,$(pathogens),../input-data/peptides/$(p)-6mers.txt) | data/concordances
	@cat $< | awk '{print "self", $$1}' > $@ ;\
	for p in $(pathogens); do \
		cat ../input-data/peptides/$$p-6mers.txt | awk -v p=$$p '{print p, $$1}' >> $@ ;\
	done

# Step 4: Compute concordances for all these pathogens at t = 4
data/concordances/peps-r4.txt : ../shared-scripts/analysis/compute-concordance.R \
	data/concordances/all-6mers.txt ../shared-scripts/analysis/concordance-functions2.R
	@echo "....Computing all peptide concordances at t = 4" && Rscript $< data/concordances/all-6mers.txt 4 1 $@


# Step 5: Compute concordances for the "self"
data/control-data-r4.txt : ../shared-scripts/analysis/peps-analyze-control.R \
	../shared-scripts/analysis/concordance-functions2.R progress/panelB-pep-control 
	@Rscript $< "data/peps-concordance/out-sim" $(nsim_control) 4 $@

progress/panelB-pep-control : ../shared-scripts/analysis/peps-concordance-control.sh $(trainset)
	@echo "...Computing peptide control concordances" && bash $< $(trainset) $(nsim_control) && touch $@

# Step 6: Make the plot.
plots/F5panelB-peps-r4.pdf : ../shared-scripts/plotting/plot-concordance.R \
	data/concordances/peps-r4.txt data/control-data-r4.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 5B for peptides" && Rscript $< data/concordances/peps-r4.txt data/control-data-r4.txt "peps" $@



# ==================== PANEL C ================================================
# Graph visualization of human/hiv peptides.

# Step 1: generate the graph data (nodes and edges)
# This does not actually involve simulations, so no dependencies except folders.
data/graphs/ghuman-hiv.dot : ../shared-scripts/analysis/concgraph-proteome.py \
	$(trainset) ../input-data/peptides/hiv-6mers.txt | data/graphs
	@echo "...Building human-hiv graph" && python3 $^ 4 "off" > $@


# Step 2: Build graph svg graphic from .dot file using graphviz program.
# some options for graphviz:
#GRAPHVIZ_OPTS=-Nlabel="" -Nstyle=filled -Npenwidth=0.0 \
#	   -Nfixedsize=true -Nwidth=0.2 -Epenwidth=4\
#	   -Gsize=14,14 -Grotate=90 -Goutputorder=edgesfirst -Gsplines=line
	   
GRAPHVIZ_OPTS=-Nlabel="" -Nstyle=filled -Nshape=circle -Npenwidth=0.0 \
	-Ncolor=lightgray -Nfixedsize=true -Nwidth=0.2 \
	-Epenwidth=5 -Ecolor=lightgray \
	-Gsize=200,200 -Goutputorder=edgesfirst -Gsplines=line

plots/F5panelC.svg : data/graphs/ghuman-hiv.dot | plots
	@echo "Plotting Figure 5C" && sfdp $(GRAPHVIZ_OPTS) -Nwidth=0.5 -Tsvg -o $@ $<
	
# Step 3: convert svg graph to pdf
plots/F5panelC.pdf : plots/F5panelC.svg | plots
	@echo "...Converting svg to pdf" && rsvg-convert $< -f pdf -o $@


# ==================== PANEL D ================================================
# Self-hiv discrimination after negative selection on different sizes of trainsets
# and with t = 3 or t = 4.

# *** STEP 1: Make the repertoires just like in figure 3.
# We first make a Makefile in data/mkfiles/panelD-repertoires.mk
# for generating the repertoires. This ensures that we can then
# make all required repertoires in parallel. 
data/mkfiles/panelD-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | \
	data/repertoires/peptides/contiguous data/repertoires/peptides/trainsets data/mkfiles
	@bash $< $(trainset) "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/peptides/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

# Ensure that the repertoires for panelB are made first, as these also are used for panelC.
progress/panelD-repertoires : data/mkfiles/panelD-repertoires.mk | progress
	@echo "...Generating repertoires for Figure5D..." && $(MAKE) -f $< && touch $@


# *** STEP 2: After generating the repertoires, we measure TCRs/test peptide.
# We first generate makefiles to count TCRs per test peptide for all simulations,
# for all test peptides of the given pathogen. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test peptides of size NTEST (which is set on top of this page). 
# For this figure, we use r=3/r=4 but take all values of ntrain.
data/mkfiles/panelD-frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelD-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@

data/mkfiles/panelD-frequencies-hiv.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/hiv-6mers.txt progress/panelD-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/hiv-6mers.txt "hiv" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/panelD-frequencies-% : data/mkfiles/panelD-frequencies-%.mk progress/panelD-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection \
	../shared-scripts/negselmodel/src/countpaths
	@echo "...Analyzing precursor frequencies for Figure 5D $* ..." && $(MAKE) -f $< && touch $@

progress/panelD-frequencies : progress/panelD-frequencies-hiv progress/panelD-frequencies-self
	@touch $@

# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/panelD-analysis.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelD-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(pep_selfsizes)" "$(rvec)" "hiv" "$(NSIM)" \
		"data/panelD-precursors" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/panelD-precursors-calc.txt : data/mkfiles/panelD-analysis.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelD-frequencies
	@echo "...Analyzing repertoires for Figure 5D" && $(MAKE) -f $< $@


# *** STEP 4: Make the plot.
plots/F5panelD-r%.pdf : ../shared-scripts/plotting/plot-top.R data/panelD-precursors-calc.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 5D r$*" && Rscript $< data/panelD-precursors-calc.txt $* "hiv" "sd" "peps" $@


# ==================== PANEL E ================================================
# For the TCR survival plot, the workflow is slightly different than for panel D.
# Here, we need to see how the number of TCRs per peptide drop as the training set
# increases, which means we need the *same* test peptides for simulations of
# different training set sizes. We therefore select those test sets first (one test
# set for each replicate simulation number, but the same test set is used across
# different values for ntrain). We then select training peptides from all peptides not
# in the test set, and generate the repertoire on those. Frequencies can then be
# analyzed in the same way as before, but using the fixed test sets chosen in step 1.
	

# Step 1: Generate test sets of NTEST peptides each, for self/hiv and each simulation
# that is needed for the output panels.
progress/panelE-test-hiv : ../input-data/peptides/hiv-6mers.txt | data/testsets progress
	@echo Making test set for hiv && \
	for sim in $$(seq 1 $(NSIM)) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk -v id="hiv" '{print id, $$1}' > data/testsets/hiv-sim$$sim.txt ;\
	done && touch $@

progress/panelE-test-self : $(trainset) | data/testsets progress
	@echo Making test set for self && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk '{print "self", $$1}' > data/testsets/self-sim$$sim.txt ;\
	done && touch $@

progress/panelE-test-all : progress/panelE-test-hiv progress/panelE-test-self | progress
	@echo merging test sets && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat data/testsets/self-sim$$sim.txt data/testsets/hiv-sim$$sim.txt > data/testsets/all-sim$$sim.txt ;\
	done && touch $@
	

# Step 2: Generate repertoires using peptides *not* in those predefined test sets.
# this is done in the same way as before, but now we make the training sets first
# (so they are not generated automatically anymore).
# First create the training sets for t=1, then just use the same ones for other t values:
progress/panelE-train-r1 : progress/panelE-test-all | data/fixtest-repertoires/peptides/trainsets progress
	@echo "...Making train sets for Figure 5E" && for sim in $$(seq 1 $(NSIM) ) ; do \
		bash ../shared-scripts/repertoires/foreign_unseen.sh data/testsets/all-sim$$sim.txt $(trainset) data/fixtest-repertoires/peptides/trainsets/tmp-sim$$sim.txt && \
		for n in $(pep_selfsizes) ; do \
			bash ../shared-scripts/repertoires/train-epitopes.sh data/fixtest-repertoires/peptides/trainsets/tmp-sim$$sim.txt $$n > data/fixtest-repertoires/peptides/trainsets/traineps-r1-R-n$$n-sim$$sim.txt ;\
		done ;\
	done && touch $@

progress/panelE-train : progress/panelE-train-r1
	@for r in $$(seq 3 4) ; do \
		for f in $$( ls data/fixtest-repertoires/peptides/trainsets/traineps-r1-R-*.txt); do \
			cp $$f $$( echo $$f | sed "s/r1/r$$r/g" ); \
		done ;\
	done && touch $@

data/mkfiles/panelE-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	progress/panelE-train | data/fixtest-repertoires/peptides/contiguous \
	data/fixtest-repertoires/peptides/trainsets data/mkfiles
	@bash $< $(trainset) "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/fixtest-repertoires/peptides/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

progress/panelE-repertoires :  data/mkfiles/panelE-repertoires.mk progress/panelE-train
	@echo "...Generating repertoires for Figure5E" && $(MAKE) -f $< && touch $@


# Step 3: Analyze precursor frequencies for the fixed self and hiv test sets.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test peptides from step1. It uses the corresponding negatively 
# selected repertoires from the previous step.
data/mkfiles/frequencies-panelE.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh \
	progress/panelE-repertoires | data/mkfiles
	@bash $< data/fixtest-repertoires/peptides "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		$(NSIM) "R" -n 6 -m "contiguous" > $@

# Use the makefile above to generate the actual files with precursor frequencies.
# They will be stored in data/fixtest-pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/fixtest-pfout/contiguous/R/output-r[t value]-peps-[pathogen]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/fixtest-pfout/contiguous/R/output-r[t value]-peps-[pathogen]-n[ntrain]-r[t value].txt
progress/panelE-frequencies : data/mkfiles/frequencies-panelE.mk \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection progress/panelE-repertoires
	@echo "...Analyzing precursor frequencies for Figure 5E..." && $(MAKE) -f $< && touch $@

# Step 4: TCR survival rates are computed from the precursor frequencies computed in the 
# previous step.
# For a given value of t, merge the data of self and xh. 
data/tsurv-r%.txt : progress/panelE-frequencies
	@echo "...Analyzing TCR survival for Figure 5E" && for n in $(pep_selfsizes) ; do \
		cat data/fixtest-pfout/contiguous/R/output-r$*-peps-all-n$$n-r$*.txt | awk -v n="$$n" -v r="$*" '{print $$0, n, r}' >> $@ ;\
	done && touch $@


# Step 5: plot the result.
plots/F5panelE-r%.pdf : ../shared-scripts/plotting/plot-tsurv2.R data/tsurv-r%.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 5E" && Rscript $< data/tsurv-r$*.txt $(NSIM) $(NTEST) "peps" $@

clean-panelE :
	rm -f figure5.pdf latex/figure5.pdf progress/panelE* data/mkfiles/*panelE* data/testsets/hiv-sim*.txt \
		data/testsets/self-sim*.txt data/testsets/all-sim*.txt data/tsurv-r3.txt plots/F5panelE*.pdf && \
	rm -rf data/fixtest-repertoires



# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	mkdir -p data 
	
data/graphs :
	mkdir -p $@

data/concordances :
	mkdir -p $@
	
data/fixtest-repertoires : 
	mkdir -p $@

data/fixtest-repertoires/peptides/contiguous :
	mkdir -p $@
	
data/fixtest-repertoires/peptides/trainsets :
	mkdir -p $@
	
data/repertoires : 
	mkdir -p $@

data/repertoires/peptides/contiguous :
	mkdir -p $@
	
data/repertoires/peptides/trainsets :
	mkdir -p $@
	
data/testsets :
	mkdir -p $@
	
data/mkfiles :
	mkdir -p $@
	
	
progress :
	mkdir -p $@
	
plots :
	mkdir -p plots

clean: latex-clean
	rm -rf data && rm -rf plots && rm -rf progress
