.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=30
NTEST=50
trainset=../input-data/peptides/self-6mers.txt

pathogens=hiv hepb hepc hcmv ebola zika mal lis vac
rvec=4
epsilon=0.01


all : figureS6.pdf

# ==================== FIGURE S6 ================================================
# This code combines all panels into a single pdf for figure S6. 
figureS6.pdf : latex/figureS6.pdf
	@cp $< $@

latex/figureS6.pdf : latex/figureS6.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS6.tex
	
	#> texlog.log


simulation-panels : $(foreach p, $(pathogens), plots/SF6-$(p).pdf)


# ==================== PANEL ================================================
# Discrimination curves for the biased trainsets.
#
# Step 1: Compute the biased sets of self-peptides.
# First AA frequency scores
data/aa-freq-scores.txt : ../shared-scripts/trainset-sampling/aa-freq-score.R $(trainset)
	Rscript $^ $@ "1 1 1 1 1 1"

# Then compute peptide probability based on those.
data/peptide-sets/pep-prob.txt : ../shared-scripts/trainset-sampling/biased-probability.R \
	data/aa-freq-scores.txt | data/peptide-sets
	@echo "...Computing biased self peptides" && Rscript $^ $@ $(epsilon)

# *** STEP 1: Make the repertoires just like in figure 5, but now with biased trainsets.
# We first make a Makefile in data/mkfiles/panelB-repertoires.mk
# for generating the repertoires. This ensures that we can then
# make all required repertoires in parallel. 

data/mkfiles/repertoires-B%.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	data/peptide-sets/pep-prob.txt | \
	data/repertoires/peptides/contiguous data/repertoires/peptides/trainsets data/mkfiles
	@bash $< data/peptide-sets/pep-prob.txt "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/peptides/contiguous -T "B" -S "1 $(NSIM)" -n 6 -m "contiguous" -b $* > $@
		

# Actually make the repertoires for both the weak (B1) and strong (B5) bias.
# For the randomly selected trainsets, copy the repertoires from figure 5.
progress/repertoires-B% : data/mkfiles/repertoires-B%.mk | progress
	@echo "...Generating repertoires for Figure7B..." && $(MAKE) -f $< && touch $@
	
progress/repertoires-R : | progress data/repertoires/peptides/contiguous \
	data/repertoires/peptides/trainsets
	@echo "...Getting repertoires from Figure 5..." && \
	cd ../figure5 && $(MAKE) progress/panelD-repertoires && \
	cp data/repertoires/peptides/contiguous/* ../figureS6/data/repertoires/peptides/contiguous/ && \
	cp data/repertoires/peptides/trainsets/* ../figureS6/data/repertoires/peptides/trainsets/ && \
	touch $@

progress/repertoires : progress/repertoires-B1 progress/repertoires-B5 progress/repertoires-R
	@touch $@

# *** STEP 2: After generating the repertoires, we measure TCRs/test peptide.
# We first generate makefiles to count TCRs per test peptide for all simulations,
# for all test peptides of the given pathogen. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test peptides of size NTEST (which is set on top of this page). 
# For this figure, we use r=3/r=4 but take all values of ntrain.
# self:
data/mkfiles/frequencies-self-R.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/frequencies-self-B1.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "B1" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/frequencies-self-B5.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "B5" \
		-n 6 -m "contiguous" -u > $@

# rest
data/mkfiles/frequencies-%-R.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/%-6mers.txt progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/$*-6mers.txt "$*" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/frequencies-%-B1.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/%-6mers.txt progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/$*-6mers.txt "$*" $(MODELDIR) $(NTEST) $(NSIM) "B1" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/frequencies-%-B5.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/%-6mers.txt progress/repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/$*-6mers.txt "$*" $(MODELDIR) $(NTEST) $(NSIM) "B5" \
		-n 6 -m "contiguous" -u > $@
		
# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/frequencies-%-R : data/mkfiles/frequencies-%-R.mk progress/repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure S6 $* R..." && $(MAKE) -f $< && touch $@

progress/frequencies-%-B1 : data/mkfiles/frequencies-%-B1.mk progress/repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure S6 $* B1..." && $(MAKE) -f $< && touch $@

progress/frequencies-%-B5 : data/mkfiles/frequencies-%-B5.mk progress/repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure S6 $* B5..." && $(MAKE) -f $< && touch $@

progress/frequencies : progress/frequencies-self-R progress/frequencies-self-B1 progress/frequencies-self-B5 \
	$(foreach p, $(pathogens), progress/frequencies-$(p)-R progress/frequencies-$(p)-B1 progress/frequencies-$(p)-B5)
	@touch $@
	
# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/analysis-R.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(pep_selfsizes)" "$(rvec)" "$(pathogens)" "$(NSIM)" \
		"data/precursors-R" > $@
		
data/mkfiles/analysis-B1.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/B1 "$(pep_selfsizes)" "$(rvec)" "$(pathogens)" "$(NSIM)" \
		"data/precursors-B1" > $@
		
data/mkfiles/analysis-B5.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/B5 "$(pep_selfsizes)" "$(rvec)" "$(pathogens)" "$(NSIM)" \
		"data/precursors-B5" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/precursors-R-calc.txt : data/mkfiles/analysis-R.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/frequencies
	@echo "...Analyzing repertoires for Figure S6" && $(MAKE) -f $< $@
data/precursors-B1-calc.txt : data/mkfiles/analysis-B1.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/frequencies
	@echo "...Analyzing repertoires for Figure S6" && $(MAKE) -f $< $@
data/precursors-B5-calc.txt : data/mkfiles/analysis-B5.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/frequencies
	@echo "...Analyzing repertoires for Figure S6" && $(MAKE) -f $< $@
	
# Combine the results with random (R) and biased (B1/B5) sets into a single file:
data/precursors-bias-calc.txt : data/precursors-B5-calc.txt \
	data/precursors-B1-calc.txt data/precursors-R-calc.txt
	cat $< | head -1 > $@ &&\
	cat $^ | grep -v type >> $@

# *** STEP 4: Make the plots.
# For this figure only t = 4.

plots/SF6-%.pdf : ../shared-scripts/plotting/plot-top.R \
	data/precursors-bias-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure S6" && Rscript $< data/precursors-bias-calc.txt 4 "$*" "sem" "peps" $@ &&\
	rm -f Rplots.pdf


# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	@cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	@mkdir -p data 
	
data/peptide-sets :
	@mkdir -p $@
	
	
data/repertoires : 
	@mkdir -p $@

data/repertoires/peptides/contiguous :
	@mkdir -p $@
	
data/repertoires/peptides/trainsets :
	@mkdir -p $@

	
data/mkfiles :
	@mkdir -p $@
	
progress :
	@mkdir -p $@
	
plots :
	@mkdir -p plots

clean: latex-clean
	@rm -rf data && rm -rf plots && rm -rf progress && rm -f figureS6.pdf
