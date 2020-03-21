.SECONDARY:
.DELETE_ON_ERROR:


all : figureS4.pdf

# ==================== FIGURE S4 ================================================
# This code combines all panels into a single pdf for figure S4. 
figureS4.pdf : latex/figureS4.pdf
	@cp $< $@

latex/figureS4.pdf : latex/figureS4.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS4.tex
	
	#> texlog.log


simulation-panels : plots/SF4.pdf


# ==================== PANEL ================================================

# Use the data already generated in figure 6.
data/precursors-greedy-calc.txt : | data
	@echo "Getting data from figure 6" &&\
	cd ../figure6 && $(MAKE) data/panelD-precursors-greedy-calc.txt &&\
	cp data/panelD-precursors-greedy-calc.txt ../figureS4/data/precursors-greedy-calc.txt


# Make the plot.
plots/SF4.pdf : ../shared-scripts/plotting/plot-top.R \
	data/precursors-greedy-calc.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/precursors-greedy-calc.txt 3 "hiv" "sem" "peps" $@



# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	@cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	@mkdir -p data 
	
	
plots :
	@mkdir -p plots

clean: latex-clean
	@rm -rf data && rm -rf plots && rm -rf progress
