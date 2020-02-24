.SECONDARY:
.DELETE_ON_ERROR:

all : plots/contiguous_lang2.pdf plots/negative-selection.pdf

plots : 
	mkdir -p $@
	
plots/contiguous_lang2.pdf : latex/contiguous_lang2.tex latex-plots | data plots
	cd data && latexmk -pdf ../$< && cp contiguous_lang2.pdf ../plots/

plots/negative-selection.pdf : negative-selection.pdf
	cp $< $@

latex-plots : plots/tcr.pdf plots/contiguous_affinity4.pdf plots/contiguous_affinity3.pdf plots/contiguous_affinity2.pdf

#plots/tcr.pdf : latex/tcr.tex | data
#	cd data && latexmk -pdf ../$< && cp $*.pdf ../plots/


plots/%.pdf : latex/%.tex | data plots
	cd data && latexmk -pdf ../$< && cp $*.pdf ../plots/


# control targets
data :
	mkdir -p data

plots : 
	mkdir -p plots

clean:
	rm -rf data/ plots/