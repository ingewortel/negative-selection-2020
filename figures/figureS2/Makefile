.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
NTEST=50
NSIM=30
languages=xh en en2 me la pd ta hi
trainset=../input-data/languages/en_t.txt
rvec=1 2 3 4 5 6
lang_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800

all : figureS2.pdf

# ==================== FIGURE S2 ==================================================
# This code combines all panels into a single pdf for figure S2. 
figureS2.pdf : latex/figureS2.pdf
	@cp $< $@

latex/figureS2.pdf : latex/figureS2.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS2.tex

# These are the simulation plots to generate
panels=panelA panelB panelC-r3 panelC-r4
simulation-panels :$(foreach p, $(panels), plots/SF2$(p).pdf)


# ==================== PANEL A ====================================================

# Step 1: compute the graph for the peptide "alphabet"
plots/Fig-specificity-peps.pdf : ../shared-scripts/plotting/recognition.R \
	../shared-scripts/plotting/mytheme.R | plots
	Rscript $< 20 "-peps"


# Step 2: annotate in latex
plots/SF2panelA.pdf : latex/figure-features-peps.tex plots/Fig-specificity-peps.pdf | plots
	cd latex && pdflatex figure-features-peps.tex && cd .. && \
	mv latex/figure-features-peps.pdf plots/SF2panelA.pdf


# ==================== PANEL B ====================================================

# The data needed for this figure was already generated in figure 5D, 
# so get the data from there
data/panelB-precursors-sim1.txt : | data
	@echo "...Analyzing data for Figure S2B" && cd ../figure5 && \
	$(MAKE) -f data/mkfiles/panelD-analysis.mk data/panelD-precursors-sim1.txt &&\
	cp data/panelD-precursors-sim1.txt ../figureS2/data/panelB-precursors-sim1.txt
	
plots/SF2panelB.pdf : ../shared-scripts/plotting/plot-pf.R data/panelB-precursors-sim1.txt \
	../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/precursors-sim1.txt 4 "hiv" 1 $@ 20 "peps"

# ==================== PANEL C ====================================================
# This is the same data as in Figure 5D, but now plotted with SEM instead of SD.

# Step 1: get the data from Figure 5
data/panelC-precursors-calc.txt : | data
	@echo "...Getting analysis from Figure 5" && \
	cd ../figure 5 $(MAKE) data/panelD-precursors-calc.txt && \
	cp data/panelD-precursors-calc.txt ../figureS2/data/panelC-precursors-calc.txt


# Step 2: make the plot
plots/SF2panelC-r%.pdf : ../shared-scripts/plotting/plot-top.R data/panelC-precursors-calc.txt \
	../shared-scripts/plotting/mytheme.R | plots
	@echo "Plotting Figure S2C r$*" && Rscript $< data/panelC-precursors-calc.txt $* "hiv" "sem" "peps" $@


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
	rm -rf data && rm -rf plots && rm -rf progress
