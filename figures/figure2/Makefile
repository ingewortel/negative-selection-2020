all : figure2.pdf

clean : 
	rm -f figure2.pdf && $(MAKE) latex-clean

latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	

figure2.pdf : latex/figure2.pdf
	cp $< $@

latex/figure2.pdf : latex/figure2.tex ../shared-scripts/latex/figure-preamble.tex figure-cartoons
	cd latex && pdflatex figure2.tex


figure-cartoons : 
	cd cartoons && $(MAKE) plots/contiguous-pep2.pdf plots/contiguous_langvspeps3.pdf
