.SECONDARY:
.DELETE_ON_ERROR:

all : figure2.pdf

figure2.pdf : latex/figure2.pdf
	cp $< $@

latex/figure2.pdf : latex/figure2.tex ../shared-scripts/latex/figure-preamble.tex \
	cartoon-panels simulation-panels
	cd latex && pdflatex figure2.tex

# Auxiliary targets
latex-clean : | latex
	cd latex && rm *.aux *.log *.pdf	
	
data : 
	mkdir -p data 
	
plots :
	mkdir -p plots

clean:
	rm -rf data && rm -rf plots

# cartoons
cartoon-panels : cartoons/plots/contiguous_lang2.pdf cartoons/plots/negative-selection.pdf

cartoons/plots/contiguous_lang2.pdf cartoons/plots/negative-selection.pdf : 
	cd cartoons && $(MAKE)

# simulations
simulation-panels : plots/Xhosa-example-r3-n500.pdf plots/Xhosa-example-r3-n0.pdf \
	plots/pf-xh-r3.pdf plots/top-xh-r3.pdf
	


plots/Xhosa-example-r3-n500.pdf : 

plots/Xhosa-example-r3-n0.pdf :

plots/pf-xh-r3.pdf :

plots/top-xh-r3.pdf : 