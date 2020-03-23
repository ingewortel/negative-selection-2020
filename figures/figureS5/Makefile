.SECONDARY:
.DELETE_ON_ERROR:




trainset=../input-data/peptides/self-6mers.txt


all : figureS5.pdf

# ==================== FIGURE S5 ================================================
# This code combines all panels into a single pdf for figure S5. 
figureS5.pdf : latex/figureS5.pdf
	@cp $< $@

latex/figureS5.pdf : latex/figureS5.tex ../shared-scripts/latex/figure-preamble.tex \
	simulation-panels
	@cd latex && pdflatex figureS5.tex
	
	#> texlog.log


simulation-panels : plots/SF5.pdf


# ==================== PANEL ================================================

# Step 1 get the optimal training set for t = 4, computed in figure 6
data/greedy-r4.txt : | data
	@echo "Getting optimal training set from figure 6" && \
	cd ../figure6 && $(MAKE) data/peptide-sets/greedy-r4.txt && \
	cp data/peptide-sets/greedy-r4.txt ../figureS5/data/

# Step 2; plot amino acid enrichment of this set versus the full "self"
plots/SF5.pdf : ../shared-scripts/plotting/plot-aa-enrichment.R $(trainset) \
	data/greedy-r4.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< $(trainset) data/greedy-r4.txt $@


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
	@rm -rf data && rm -rf plots && rm -rf progress && rm -f figureS5.pdf
