.SECONDARY:
.DELETE_ON_ERROR:

rvec=1 2 3 4 5 6
MODELDIR=../shared-scripts/negselmodel/src
pep_selfsizes=0 1 100 1000 5000 10000 25000 40000 55000 70000 100000 150000 200000 250000 260000
NSIM=30
NTEST=50
trainset=../input-data/peptides/self-6mers.txt
pathogens=hiv hepb hepc hcmv ebola zika mal lis vac

all : figureS3.pdf

# ==================== FIGURE S3 ================================================
# This code combines all panels into a single pdf for figure S3. 
figureS3.pdf : latex/figureS3.pdf
	@cp $< $@

latex/figureS3.pdf : latex/figureS3.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS3.tex
	
	#> texlog.log


simulation-panels : $(foreach p, $(pathogens), plots/SF3panelA-$(p).pdf plots/SF3panelB-$(p).pdf)


# ==================== PANEL A ================================================
# concordances

# Step 1: get 6mers from figure 5
data/all-6mers.txt : | data
	@echo "...Getting 6mers from Figure 5"
	cd ../figure5 && $(MAKE) data/concordances/all-6mers.txt &&\
	cp data/concordances/all-6mers.txt ../figureS3/data/all-6mers.txt

# Step 2: compute concordances for all thresholds t
data/concordances-r%.txt : ../shared-scripts/analysis/compute-concordance.R \
	../shared-scripts/analysis/concordance-functions2.R data/all-6mers.txt | data/concordances
	Rscript $< data/all-6mers.txt $* 1 $@

progress/concordances : $(foreach r, $(rvec), data/concordances-r$(r).txt) | progress
	@touch $@
	
# Step 3: combine data per pathogen
data/concordances-%.txt : progress/concordances
	for r in $(rvec); do \
		conc=$$(cat data/concordances-r$$r.txt | awk '$$1 == "$*" {print $$2}');\
		echo $$r $$conc >> $@ ;\
	done

# Step 4: make plots
plots/SF3panelA-%.pdf : ../shared-scripts/plotting/plot-concordance-rvalue.R \
	data/concordances-%.txt  ../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure S3A"  && \
	Rscript $< data/concordances-$*.txt $@


# ==================== PANEL B ================================================
# pathogens: discrimination curves for pathogens other than hiv.

# Step 1: get repertoires from figure 5D.
progress/panelB-repertoires : | progress data/repertoires/peptides/contiguous \
	data/repertoires/peptides/trainsets
	@echo "...Copying repertoires from figure 5D" &&\
	cd ../figure5 && $(MAKE) progress/panelD-repertoires && \
	cp data/repertoires/peptides/contiguous/* ../figureS3/data/repertoires/peptides/contiguous/ && \
	cp data/repertoires/peptides/trainsets/* ../figureS3/data/repertoires/peptides/trainsets/ && \
	touch $@

# Step 2: analyze precursor frequencies
# *** STEP 2: After generating the repertoires, we measure TCRs/test peptide.
# We first generate makefiles to count TCRs per test peptide for all simulations,
# for all test peptides of the given pathogen. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test peptides of size NTEST (which is set on top of this page). 
data/mkfiles/panelB-frequencies-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "4" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@

data/mkfiles/panelB-frequencies-%.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/peptides/%-6mers.txt progress/panelB-repertoires | data/mkfiles
	@bash $< data/repertoires/peptides "$(pep_selfsizes)" "4" \
		../input-data/peptides/$*-6mers.txt "$*" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -u > $@

# Use these makefiles to generate the actual files with motifs/test string. The progress/X
# files just indicate that the step has been completed, while the actual data
# will be stored in data/pfout/contiguous/R/. See also previous panel.
progress/panelB-frequencies-% : data/mkfiles/panelB-frequencies-%.mk progress/panelB-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection \
	../shared-scripts/negselmodel/src/countpaths
	@echo "...Analyzing precursor frequencies for Figure S3B $* ..." && $(MAKE) -f $< && touch $@

progress/panelB-frequencies : progress/panelB-frequencies-self \
	$(foreach p, $(pathogens), progress/panelB-frequencies-$(p))
	@touch $@

# *** STEP 3: Further analyze the data; Precursor frequencies computed for training sets 
# of different sizes (previous step) are here processed further to give median + IQR over 
# the recognition of different test strings (but still in a single simulation for each ntrain). 
# First again generate a makefile:
data/mkfiles/panelB-analysis.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh \
	progress/panelB-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "$(pep_selfsizes)" "4" "$(pathogens)" "$(NSIM)" \
		"data/panelB-precursors" > $@

# This file is produced by the makefile, by analyzing all the separate trainsets and merging
# the results:
data/panelB-precursors-calc.txt : data/mkfiles/panelB-analysis.mk \
	../shared-scripts/analysis/analyze-precursors3.R progress/panelB-frequencies
	@echo "...Analyzing repertoires for Figure S3B" && $(MAKE) -f $< $@


# *** STEP 4: Make the plot.
plots/SF3panelB-%.pdf : ../shared-scripts/plotting/plot-top.R data/panelB-precursors-calc.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure S3B $*" && Rscript $< data/panelB-precursors-calc.txt 4 "$*" "sd" "peps" $@



# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	@cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	@mkdir -p data 
	
data/repertoires/peptides/contiguous :
	@mkdir -p $@
	
data/repertoires/peptides/trainsets :
	@mkdir -p $@
	
data/mkfiles :
	@mkdir -p $@
	
	
data/concordances :
	@mkdir -p $@

progress : 
	@mkdir -p progress
	
	
plots :
	@mkdir -p plots

clean: latex-clean
	@rm -rf data && rm -rf plots && rm -rf progress
