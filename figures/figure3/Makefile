.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
panelC_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800
panelD_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800
NSIM=30
NTEST=50
trainset=../input-data/languages/en_t.txt

all : figure3.pdf

# ==================== FIGURE 3 ==================================================
# This code combines all panels into a single pdf for figure 3. 
figure3.pdf : latex/figure3.pdf
	@cp $< $@

latex/figure3.pdf : latex/figure3.tex ../shared-scripts/latex/figure-preamble.tex \
	cartoon-panels simulation-panels
	@cd latex && pdflatex figure3.tex > texlog.log

# These are the simulation plots to generate
panels=panelB1 panelB2 panelC panelD
simulation-panels :$(foreach p, $(panels), plots/F3$(p).pdf)

# ==================== CARTOONS (panel A) ========================================
# The first two panels of this figure do not involve simulations but show examples
# of how the model/simulations work.
cartoon-panels : cartoons/negative-selection-2.pdf


# ==================== PANEL B ====================================================

# *** STEP 1: In the first step of the simulation, training strings are compared against 
# the total repertoire to obtain a negatively selected repertoire.
#
# We first make a Makefile in data/mkfiles/panelX-repertoires.mk
# for generating the repertoires for each of the panels. This ensures that we can then
# make all required repertoires in parallel. 
data/mkfiles/panelB1-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | \
	data/repertoires/languages/contiguous data/repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "0" "3" $(MODELDIR) data/repertoires/languages/contiguous \
		-T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@
	
data/mkfiles/panelB2-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | \
	data/repertoires/languages/contiguous data/repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "500" "3" $(MODELDIR) data/repertoires/languages/contiguous \
		-T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

# The progress/xx files are just empty files to indicate that the repertoires have been
# completed. The actual repertoires will end up in an automatically generated directory 
# data/repertoires/languages/contiguous, and the used training sets (also generated 
# automatically) will be saved in data/repertoires/languages/trainsets.
progress/panelB1-repertoires : data/mkfiles/panelB1-repertoires.mk | progress
	@echo ".......Generating repertoires for Figure3B1..." && $(MAKE) -f $< && touch $@

progress/panelB2-repertoires : data/mkfiles/panelB2-repertoires.mk | progress
	@echo ".......Generating repertoires for Figure3B2..." && $(MAKE) -f $< && touch $@

progress/panelB-repertoires : progress/panelB1-repertoires progress/panelB2-repertoires
	@touch $@

# *** STEP 2: After generating the repertoires, we test how many motifs still recognize
# the strings in the test set.
# We first again generate makefiles to count motifs per test string for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
# For this figure, we only need the analyiss for n = 0 or 500, r = 3.
data/mkfiles/panelB-frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelB-repertoires | data/mkfiles
	bash $< data/repertoires/languages "0 500" "3" $(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

data/mkfiles/panelB-frequencies-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt progress/panelB-repertoires | data/mkfiles
	bash $< data/repertoires/languages "0 500" "3" \
		../input-data/languages/xh-unseen.txt "xh" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
progress/panelB-frequencies-% : data/mkfiles/panelB-frequencies-%.mk progress/panelB-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies : progress/panelB-frequencies-xh progress/panelB-frequencies-self
	@touch $@

# *** STEP 3: For this panel, we just take the raw motif counts for the test strings of
# a single simulation, and convert them to a frequency (motifs/(10^6 total motifs left)).
# We plot those here:
plots/F3panelB1.pdf : ../shared-scripts/plotting/plot-pf-dot.R progress/panelB-frequencies \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel B1" && Rscript $< "data/pfout/contiguous/R" "lang" "3" "xh" "0" $(NSIM) $@ 220

plots/F3panelB2.pdf : ../shared-scripts/plotting/plot-pf-dot.R progress/panelB-frequencies \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel B2" && Rscript $< "data/pfout/contiguous/R" "lang" "3" "xh" "500" $(NSIM) $@ 220


# ==================== PANEL C ====================================================

# *** STEP 1: Make the repertoires just like in panel B, but now for different training
# set sizes. Also already make all 30 simulations (for logistic reasons).
# We first make a Makefile in data/mkfiles/panelX-repertoires.mk
# for generating the repertoires for each of the panels. This ensures that we can then
# make all required repertoires in parallel. 
data/mkfiles/panelC-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | \
	data/repertoires/languages/contiguous data/repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "$(panelC_selfsizes)" "3" $(MODELDIR) \
		data/repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

# Ensure that the repertoires for panelB are made first, as these also are used for panelC.
progress/panelC-repertoires : data/mkfiles/panelC-repertoires.mk progress/panelB-repertoires
	@echo ".......Generating repertoires for Figure3C..." && $(MAKE) -f $< && touch $@
	
# *** STEP 2: After generating the repertoires, we again test motifs/test string, but now
# for all trainset sizes. 
# We first again generate makefiles to count motifs per test string for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
# For this figure, we only need the r=3 but take all values of ntrain.
data/mkfiles/panelC-frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelC-repertoires | data/mkfiles
	bash $< data/repertoires/languages "$(panelC_selfsizes)" "3" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

data/mkfiles/panelC-frequencies-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt progress/panelC-repertoires | data/mkfiles
	bash $< data/repertoires/languages "$(panelC_selfsizes)" "3" \
		../input-data/languages/xh-unseen.txt "xh" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/panelC-frequencies-% : data/mkfiles/panelC-frequencies-%.mk progress/panelC-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

progress/panelC-frequencies : progress/panelC-frequencies-xh progress/panelC-frequencies-self
	@touch $@


# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/panelC-analysis.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelC-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(panelC_selfsizes)" "3" "xh" "$(NSIM)" \
		"data/panelC-precursors" -l "-lang" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/panelC-precursors-sim1.txt : data/mkfiles/panelC-analysis.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelC-frequencies
	@echo "Analyzing repertoires for Figure 3C" && $(MAKE) -f $< $@


# *** STEP 4: Make the plot.
plots/F3panelC.pdf : ../shared-scripts/plotting/plot-pf.R data/panelC-precursors-sim1.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel C" && Rscript $< data/panelC-precursors-sim1.txt 3 "xh" 1 $@ 0



# ==================== PANEL D ====================================================


# *** STEP 1: We use the same repertoires as for panel C.
# Ensure that the repertoires for panelC are made first, as these also are used for panelD.
progress/panelD-repertoires :  progress/panelC-repertoires
	@touch $@
	
# *** STEP 2: Use the same frequencies as for panel C.
progress/panelD-frequencies : progress/panelC-frequencies
	@touch $@

# *** STEP 3: Analyze precursor counts.
# Precursor frequencies computed for training sets of different sizes are now analyzed 
# differently to yield the 	%foreign in the best recognized strings. This is done for 
# 30 simulations, yielding a mean and SD/SEM for each simulation.
data/mkfiles/panelD-analysis.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelD-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(panelD_selfsizes)" "3" "xh" "$(NSIM)" \
		"data/panelD-precursors" -l "-lang" > $@

data/panelD-precursors-calc.txt : data/mkfiles/panelD-analysis.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelD-frequencies
	@echo "Analyzing repertoires for Figure 3D" && $(MAKE) -f $< $@

# *** STEP 4: Make the plot.
plots/F3panelD.pdf : ../shared-scripts/plotting/plot-top.R data/panelD-precursors-calc.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel D" && Rscript $< data/panelD-precursors-calc.txt 3 "xh" \
		"sd" "lang" $@ && rm -f Rplots.pdf




# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	
	
	
progress :
	mkdir -p $@

data : 
	mkdir -p data 
	
data/mkfiles :
	mkdir -p $@
	
data/repertoires : 
	mkdir -p $@
	
data/repertoires/languages/contiguous : 
	mkdir -p $@
	
data/repertoires/languages/trainsets : 
	mkdir -p $@
	
data/output :
	mkdir -p $@
	
plots :
	mkdir -p plots

clean: latex-clean
	rm -rf data && rm -rf plots && rm -rf progress && rm -f figure3.pdf
