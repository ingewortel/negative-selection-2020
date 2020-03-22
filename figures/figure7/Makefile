.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=30
NTEST=50
trainset=../input-data/peptides/self-6mers.txt
trainset2=../input-data/peptides/self1-6mers.txt# this trainset now contains most of all
												# self-peptides, but some are left out
												# for the "other self" experiment.
pathogens=hiv hepb hepc hcmv ebola zika mal lis vac
rvec=4
epsilon=0.01


all : figure7.pdf

# ==================== FIGURE 7 ================================================
# This code combines all panels into a single pdf for figure 5. 
figure7.pdf : latex/figure7.pdf
	@cp $< $@

latex/figure7.pdf : latex/figure7.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figure7.tex


# These are the simulation plots to generate
panels=panelA panelB-hiv panelB-self2 panelC
simulation-panels :$(foreach p, $(panels), plots/F7$(p).pdf)


# ==================== PANEL A ================================================
# Plot of frequency score versus exchangeability.
#
# Step 1: Get the computed exchangeabilities from figure 6, and filter the self-ones
data/exch-all.txt : | data
	@echo "... Getting exchangeabilities from Figure 6" &&\
	cd ../figure6 && $(MAKE) $@ && cp data/exch-all.txt ../figure7/data/exch-all.txt

data/exch-self.txt : data/exch-all.txt
	cat $< | awk '$$2==0{print}' > $@

# Step 2: Compute the AA frequency scores for the self peptides
data/self-scores.txt : ../shared-scripts/trainset-sampling/aa-freq-score.R $(trainset) | data
	Rscript $^ $@ "1 1 1 1 1 1"

# Step 3: Merge with exchangeabilities into a single file
data/self-degreescore.txt : data/self-scores.txt data/exch-self.txt | data
	cat $< | sort | awk '{print $$2}' | paste data/exch-self.txt - | awk '{print $$1, $$3, $$6}' > $@

# Step 4: Plot the result.
plots/F7panelA.pdf : ../shared-scripts/plotting/plot-degree-vs-score-2.R \
	data/self-degreescore.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/self-degreescore.txt $@ "peps"


# ==================== PANEL B ================================================
# Discrimination curves for the biased trainsets.
#
# Step 1: Compute the biased sets of self-peptides.
# First AA frequency scores
data/aa-freq-scores.txt : ../shared-scripts/trainset-sampling/aa-freq-score.R $(trainset2)
	Rscript $^ $@ "1 1 1 1 1 1"

# Then compute peptide probability based on those.
data/peptide-sets/pep-prob.txt : ../shared-scripts/trainset-sampling/biased-probability.R \
	data/aa-freq-scores.txt | data/peptide-sets
	@echo "...Computing biased self peptides" && Rscript $^ $@ $(epsilon)

# *** STEP 1: Make the repertoires just like in figure 5, but now with biased trainsets.
# We first make a Makefile in data/mkfiles/panelB-repertoires.mk
# for generating the repertoires. This ensures that we can then
# make all required repertoires in parallel. 
#
# Also make new 'random' repertoires with the same slightly smaller 
# training set (for consistency with the self-other self figure). But note that the
# results for the other pathogens are pretty much the same as when the entire self
# trainset is used.
data/mkfiles/panelB-repertoires-B%.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	data/peptide-sets/pep-prob.txt | \
	data/repertoires/peptides/contiguous data/repertoires/peptides/trainsets data/mkfiles
	@bash $< data/peptide-sets/pep-prob.txt "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/peptides/contiguous -T "B" -S "1 $(NSIM)" -n 6 -m "contiguous" -b $* > $@
		
data/mkfiles/panelB-repertoires-R.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh \
	$(trainset2) | \
	data/repertoires/peptides/contiguous data/repertoires/peptides/trainsets data/mkfiles
	@bash $< $(trainset2) "$(pep_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/peptides/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" > $@

# Actually make the repertoires for both the weak (B1) and strong (B5) bias, and for 
# random selection (R)
progress/panelB-repertoires-B% : data/mkfiles/panelB-repertoires-B%.mk | progress
	@echo "...Generating repertoires for Figure7B..." && $(MAKE) -f $< && touch $@
	
progress/panelB-repertoires-R : data/mkfiles/panelB-repertoires-R.mk | progress
	@echo "...Generating repertoires for Figure7B..." && $(MAKE) -f $< && touch $@

progress/panelB-repertoires : progress/panelB-repertoires-B1 \
	progress/panelB-repertoires-B5 progress/panelB-repertoires-R
	@touch $@

# *** STEP 2: After generating the repertoires, we measure TCRs/test peptide.
# We first generate makefiles to count TCRs per test peptide for all simulations,
# for all test peptides of the given pathogen. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test peptides of size NTEST (which is set on top of this page). 
# For this figure, we use r=3/r=4 but take all values of ntrain.
# self:
data/mkfiles/panelB-frequencies-self-R.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset2) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/panelB-frequencies-self-B1.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset2) "self" $(MODELDIR) $(NTEST) $(NSIM) "B1" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/panelB-frequencies-self-B5.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		$(trainset2) "self" $(MODELDIR) $(NTEST) $(NSIM) "B5" \
		-n 6 -m "contiguous" -u > $@

# hiv:
data/mkfiles/panelB-frequencies-hiv-R.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/hiv-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/hiv-6mers.txt "hiv" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/panelB-frequencies-hiv-B1.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/hiv-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/hiv-6mers.txt "hiv" $(MODELDIR) $(NTEST) $(NSIM) "B1" \
		-n 6 -m "contiguous" -u > $@
		
data/mkfiles/panelB-frequencies-hiv-B5.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/hiv-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/hiv-6mers.txt "hiv" $(MODELDIR) $(NTEST) $(NSIM) "B5" \
		-n 6 -m "contiguous" -u > $@
		
# "other self"
data/mkfiles/panelB-frequencies-self2-R.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/self2-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/self2-6mers.txt "self2" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@		
data/mkfiles/panelB-frequencies-self2-B1.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/self2-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/self2-6mers.txt "self2" $(MODELDIR) $(NTEST) $(NSIM) "B1" \
		-n 6 -m "contiguous" -u > $@
data/mkfiles/panelB-frequencies-self2-B5.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/self2-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "$(rvec)" \
		../input-data/peptides/self2-6mers.txt "self2" $(MODELDIR) $(NTEST) $(NSIM) "B5" \
		-n 6 -m "contiguous" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/panelB-frequencies-%-R : data/mkfiles/panelB-frequencies-%-R.mk progress/panelB-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure 7B $* R..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies-%-B1 : data/mkfiles/panelB-frequencies-%-B1.mk progress/panelB-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure 7B $* B1..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies-%-B5 : data/mkfiles/panelB-frequencies-%-B5.mk progress/panelB-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/contiguous-negative-selection $(MODELDIR)/countpaths
	@echo "...Analyzing precursor frequencies for Figure 7B $* B5..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies : progress/panelB-frequencies-hiv-B1 progress/panelB-frequencies-self-B1 \
	progress/panelB-frequencies-hiv-B5 progress/panelB-frequencies-self-B5 \
	progress/panelB-frequencies-self2-B1 progress/panelB-frequencies-self2-B5 \
	progress/panelB-frequencies-self2-R progress/panelB-frequencies-self-R \
	progress/panelB-frequencies-hiv-R
	@touch $@
	
# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/panelB-analysis-R.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelB-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(pep_selfsizes)" "$(rvec)" "hiv self2" "$(NSIM)" \
		"data/panelB-precursors-R" > $@
		
data/mkfiles/panelB-analysis-B1.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelB-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/B1 "$(pep_selfsizes)" "$(rvec)" "hiv self2" "$(NSIM)" \
		"data/panelB-precursors-B1" > $@
		
data/mkfiles/panelB-analysis-B5.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelB-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/B5 "$(pep_selfsizes)" "$(rvec)" "hiv self2" "$(NSIM)" \
		"data/panelB-precursors-B5" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/panelB-precursors-R-calc.txt : data/mkfiles/panelB-analysis-R.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelB-frequencies
	@echo "...Analyzing repertoires for Figure 7B" && $(MAKE) -f $< $@
data/panelB-precursors-B1-calc.txt : data/mkfiles/panelB-analysis-B1.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelB-frequencies
	@echo "...Analyzing repertoires for Figure 7B" && $(MAKE) -f $< $@
data/panelB-precursors-B5-calc.txt : data/mkfiles/panelB-analysis-B5.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelB-frequencies
	@echo "...Analyzing repertoires for Figure 7B" && $(MAKE) -f $< $@
	
# Combine the results with random (R) and biased (B1/B5) sets into a single file:
data/panelB-precursors-bias-calc.txt : data/panelB-precursors-B5-calc.txt \
	data/panelB-precursors-B1-calc.txt data/panelB-precursors-R-calc.txt
	cat $< | head -1 > $@ &&\
	cat $^ | grep -v type >> $@

# *** STEP 4: Make the plots.
# For this figure only t = 4.
plots/F7panelB-hiv.pdf : ../shared-scripts/plotting/plot-top.R \
	data/panelB-precursors-bias-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 7B" && Rscript $< data/panelB-precursors-bias-calc.txt 4 "hiv" "sem" "peps" $@
plots/F7panelB-self2.pdf : ../shared-scripts/plotting/plot-top.R \
	data/panelB-precursors-bias-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure 7B" && Rscript $< data/panelB-precursors-bias-calc.txt 4 "self2" "sem" "peps" $@



# ==================== PANEL C ================================================
# Analyzing all pathogens on the R/B1/B5 training sets. We repeat the above, but now
# using the entire "self" trainset again.
# The data for this figure will be generated in figure S6, since this plot just 
# represents a single point of all those plots.

data/precursors-bias-calc.txt : | data
	@echo "...Getting data from figure S6" && \
	cd ../figureS6 && $(MAKE) $@ && cp $@ ../figure7/data/precursors-bias-calc.txt

plots/F7panelC.pdf : ../shared-scripts/plotting/plot-compare-pathogens.R \
	data/precursors-bias-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< $*




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
	@rm -rf data && rm -rf plots && rm -rf progress
