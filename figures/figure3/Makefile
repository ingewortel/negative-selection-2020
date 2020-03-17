.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
panelC_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800
panelD_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800
NSIM=30
NTEST=50
trainset=../input-data/languages/en_t.txt

all : figure3.pdf

# ==================== FIGURE 3 ================================================
# This code combines all panels into a single pdf for figure 3. 
figure3.pdf : latex/figure3.pdf
	@cp $< $@

latex/figure3.pdf : latex/figure3.tex ../shared-scripts/latex/figure-preamble.tex \
	cartoon-panels simulation-panels
	@cd latex && pdflatex figure3.tex > texlog.log


# ==================== CARTOONS ================================================
# The first two panels of this figure do not involve simulations but show examples
# of how the model/simulations work.
cartoon-panels : cartoons/negative-selection-2.pdf

# ==================== PANELS ================================================
# These are the simulation plots to generate
panels=panelB1 panelB2 panelC panelD
simulation-panels :$(foreach p, $(panels), plots/F3$(p).pdf)


plots/F3panelB1.pdf : ../shared-scripts/plotting/plot-pf-dot.R data/frequencies ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel B1" && Rscript $< "data/pfout/contiguous/R" "lang" "3" "xh" "0" $(NSIM) $@ 220

plots/F3panelB2.pdf : ../shared-scripts/plotting/plot-pf-dot.R data/frequencies ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel B2" && Rscript $< "data/pfout/contiguous/R" "lang" "3" "xh" "500" $(NSIM) $@ 220

plots/F3panelC.pdf : ../shared-scripts/plotting/plot-pf.R data/precursors-sim1.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel C" && Rscript $< data/precursors-sim1.txt 3 "xh" 1 $@ 0

plots/F3panelD.pdf : ../shared-scripts/plotting/plot-top.R data/precursors-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Fig3 panel D" && Rscript $< data/precursors-calc.txt 3 "xh" "sd" "lang" $@ && rm -f Rplots.pdf

# ==================== SIMULATION DATA 3: PLOTTED DATA ============================
# Processed data used for plotting, using the data from SIMULATION DATA 2.
# Panel B: Computed precursor frequencies for the test strings are used directly; so this
# 			step is not required. 
# Panel C: Precursor frequencies are computed for training sets of different sizes, and 
# 			processed further to give median + IQR of a single simulation for each ntrain.
# Panel D: Similar, but now precursor frequencies are analyzed differently to yield the
# 			%foreign in the best recognized strings. This is done for 30 simulations, 
# 			yielding a mean, SD, and SEM for each simulation.

# Analyzed version of the raw data for panel C
data/precursors-sim1.txt : data/analysis.mk ../shared-scripts/analysis/analyze-precursors3.R data/frequencies
	@echo "Analyzing repertoires for Figure 3C" && $(MAKE) -f $< $@

# Analyzed version of the raw data for panel D
data/precursors-calc.txt : data/analysis.mk ../shared-scripts/analysis/analyze-precursors3.R data/frequencies
	@echo "Analyzing repertoires for Figure 3D" && $(MAKE) -f $< $@

data/analysis.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh data/output
	@bash $< data/pfout/contiguous/R "$(panelD_selfsizes)" "3" "xh" "$(NSIM)" "data/precursors" -l "-lang" > $@


# ==================== SIMULATION DATA 2: PRECURSOR FRQUENCIES ======================
# Raw data from the simulations, where the remaining number of motifs for each test 
# string are counted. This uses the repertoires generated in SIMULATION DATA 1. 

# This generates an empty file data/frequencies which just shows that all the 
# simulations have been performed. In turn, the files data/frequencies-xh and
# data/frequencies-self indicate that all files in the corresponding makefiles 
# data/frequencies-self.mk and data/frequencies-xh.mk have been generated.
data/frequencies : data/frequencies-xh data/frequencies-self
	@echo "Done analyzing precursor frequencies for Figure 3." && touch $@

# Use the makefiles below to generate the actual files with precursor frequencies.
# They will be stored in data/pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
data/frequencies-% : data/frequencies-%.mk data/repertoires/all \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
data/frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	data/repertoires/all
	@rm -f data/pfout/contiguous/R/output*-self-*.txt;\
	bash $< data/repertoires/languages "$(panelD_selfsizes)" "3" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

data/frequencies-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt data/repertoires/all 
	@rm -f data/pfout/contiguous/R/output*-xh-*.txt;\
	bash $< data/repertoires/languages "$(panelD_selfsizes)" "3" \
		../input-data/languages/xh-unseen.txt "xh" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@


# ==================== SIMULATION DATA 1: REPERTOIRES ===============================
# In the first step of the simulation, training strings are compared against the total
# repertoire to obtain a negatively selected repertoire.
#
# Below, we first make separate Makefiles in data/repertoires/panelX-repertoires.mk
# for generating the repertoires for each of the panels. This ensures that we can then
# make those repertoires in parallel. The actual repertoires will end up in an 
# automatically generated directory data/repertoires/languages/contiguous, and the used
# training sets (also generated automatically) will be saved in data/repertoires/languages/trainsets.

data/repertoires/all : $(foreach p, $(panels), data/repertoires/$(p))
	@echo "Done generating repertoires for Figure 3." && touch $@

data/repertoires/panelB1 : data/repertoires/panelB1-repertoires.mk
	@echo ".......Generating repertoires for Figure3B1..." && $(MAKE) -f $< && touch $@

data/repertoires/panelB2 : data/repertoires/panelB2-repertoires.mk
	@echo ".......Generating repertoires for Figure3B2..." && $(MAKE) -f $< && touch $@

# Ensure that the repertoires for panelB are made first, as these also are used for panelC.
data/repertoires/panelC : data/repertoires/panelC-repertoires.mk data/repertoires/panelB1 data/repertoires/panelB2
	@echo ".......Generating repertoires for Figure3C..." && $(MAKE) -f $< && touch $@

# Ensure that the repertoires for panelC are made first, as these also are used for panelD.
data/repertoires/panelD :  data/repertoires/panelD-repertoires.mk data/repertoires/panelC
	@echo ".......Generating repertoires for Figure3D..." && $(MAKE) -f $< && touch $@

# The following are written by the script ../shared-scripts/repertoires/loop-repertoires-makeout.sh.
# They generate makefiles for automatically generating all the required repertoires.
data/repertoires/panelB1-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets
	@bash $< $(trainset) "0" "3" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 1" -n 6 -m "contiguous" -l "-lang" > $@
	
data/repertoires/panelB2-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets
	@bash $< $(trainset) "500" "3" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 1" -n 6 -m "contiguous" -l "-lang" > $@

data/repertoires/panelC-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets
	@bash $< $(trainset) "$(panelC_selfsizes)" "3" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 1" -n 6 -m "contiguous" -l "-lang" > $@

data/repertoires/panelD-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets
	@bash $< $(trainset) "$(panelD_selfsizes)" "3" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@


# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	mkdir -p data 
	
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
	rm -rf data && rm -rf plots
