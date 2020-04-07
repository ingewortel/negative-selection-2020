.SECONDARY:
.DELETE_ON_ERROR:

all : plots/contiguous-pep2.pdf plots/contiguous_langvspeps3.pdf


plots/contiguous-pep2.pdf : latex/contiguous-pep2.tex latex-plots | data plots
	cd data && latexmk -pdf ../$< && cp contiguous-pep2.pdf ../plots/



latex-plots : plots/tcr-peps.pdf plots/tcr.pdf plots/contiguous_affinity4-peps.pdf plots/contiguous_affinity2-peps.pdf

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