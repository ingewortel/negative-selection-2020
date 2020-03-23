.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
NTEST=50
NSIM=30
languages=xh en en2 me la pd ta hi
trainset=../input-data/languages/en_t.txt
rvec=1 2 3 4 5 6
lang_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800

all : figureS1.pdf

# ==================== FIGURE S1 ==================================================
# This code combines all panels into a single pdf for figure S1. 
figureS1.pdf : latex/figureS1.pdf
	@cp $< $@

latex/figureS1.pdf : latex/figureS1.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS1.tex > texlog.log

# These are the simulation plots to generate
panels=panelA panelB panelC
simulation-panels :$(foreach p, $(panels), plots/SF1$(p).pdf)


# ==================== PANEL A ====================================================
# This is a replotted version of Figure 3E.

# Step 1 : Copy the data from figure 3.

data/panelA-precursors-calc.txt : | data 
	@echo "Getting figure 3 data for SFigure 1A" && \
	cd ../figure3 $(MAKE) data/panelD-precursors-calc.txt && \
	cp data/panelD-precursors-calc.txt ../figureS1/data/panelA-precursors-calc.txt

# *** STEP 4: Make the plot.
plots/SF1panelA.pdf : ../shared-scripts/plotting/plot-top.R data/panelA-precursors-calc.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting SFig1 panel A" && Rscript $< data/panelA-precursors-calc.txt 3 "xh" \
		"sem" "lang" $@ && rm -f Rplots.pdf


# ==================== PANEL B ====================================================
# This is a replotted version of Figure 3E.

# Step 1: Same repertoires as figure 3, so copy them here
progress/panelB-repertoires : | data/repertoires progress
	cd ../figure3 && $(MAKE) progress/panelD-repertoires && cd ../figureS1 && \
	cp -r ../figure3/data/repertoires/* data/repertoires/ && touch $@

# Step 2: Compute affinities, first make a makefile and then run these to analyze
# all data
data/affn-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt progress/panelB-repertoires
	@echo "...Computing affinities " && \
	bash $< data/repertoires/languages "500" "3" ../input-data/languages/xh-unseen.txt \
		"xh" $(MODELDIR) $(NTEST) $(NSIM) "R" -n 6 -m "contiguous" \
		-l "-lang" -u | sed 's:persim:persim-allR:g' | sed 's:output:affn:g' > $@

data/affn-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	$(trainset) progress/panelB-repertoires
	@echo "...Computing affinities " && \
	bash $< data/repertoires/languages "500" "3" $(trainset) "self" $(MODELDIR) \
		$(NTEST) $(NSIM) "R" -n 6 -m "contiguous" -l "-lang" \
		-u | sed 's:persim:persim-allR:g' | sed 's:output:affn:g' > $@
		
		
# Ensure that the affinites are computed using the makefiles
progress/affn-% : data/affn-%.mk progress/panelB-repertoires ../shared-scripts/analysis/precursor-counts-persim-allR.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh $(MODELDIR)/contiguous-negative-selection-lang \
	$(MODELDIR)/countpaths | progress
	$(MAKE) -f $< && touch $@

progress/panelB-affn : progress/affn-xh progress/affn-self

plots/SF1panelB.pdf : ../shared-scripts/plotting/plot-affn-lang-example.R progress/panelB-affn | plots
	Rscript $< data/pfout/contiguous/R/affn-r3-lang-self-n500-r3.txt \
		data/pfout/contiguous/R/affn-r3-lang-xh-n500-r3.txt 30 3 $@


# ==================== PANEL C ====================================================

# Step 1: add repertoires for other values of the threshold t

# Makefile for generating the repertoires:
data/mkfiles/panelC-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | \
	data/repertoires/languages/contiguous data/repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "$(lang_selfsizes)" "$(rvec)" $(MODELDIR) \
		data/repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

# Ensure that the repertoires for panelB are made first, as these also are used for panelC.
progress/panelC-repertoires : data/mkfiles/panelC-repertoires.mk progress/panelB-repertoires
	@echo ".......Generating repertoires for Figure S1C..." && $(MAKE) -f $< && touch $@
	

# Step 2: analyze deletion in the repertoires of step 1
data/deletion-perc.txt : ../shared-scripts/analysis/count-survivors-repertoiretype.sh \
	progress/panelC-repertoires | data
	@echo "...Analyzing deletion of motifs in repertoire" && \
	bash $< data/repertoires/languages/contiguous $(MODELDIR) "R" "$(rvec)" "0 1 10 50 100 250 500 800" $(NSIM) > $@


# Step 2: plot
plots/SF1panelC.pdf : ../shared-scripts/plotting/plot-deletion-versus-t.R \
	data/deletion-perc.txt | plots
	Rscript $^ 500 R $@


# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	

data/repertoires :
	mkdir -p $@
	
data/mkfiles :
	mkdir -p $@
	
data/repertoires/languages/contiguous :
	@mkdir -p $@

data/repertoires/languages/trainsets :
	@mkdir -p $@
		
progress :
	mkdir -p $@

data : 
	mkdir -p data 
	

plots :
	mkdir -p plots

clean: latex-clean
	rm -rf data && rm -rf plots && rm -rf progress && rm -f figureS1.pdf
