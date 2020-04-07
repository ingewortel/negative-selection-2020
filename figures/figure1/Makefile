all : figure1.pdf


latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	

figure1.pdf : latex/figure1.pdf
	cp $< $@

latex/figure1.pdf : latex/figure1.tex ../shared-scripts/latex/figure-preamble.tex cartoons/tolerance-discrimination.pdf
	cd latex && pdflatex figure1.tex

clean :
	rm -f figure1.pdf && $(MAKE) latex-clean
