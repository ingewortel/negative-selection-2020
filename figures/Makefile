all : all-figures

figs=1 2 3 4 5 6



all-figures : $(foreach f, $(figs), figure$(f)/figure$(f).pdf)

figure1/figure1.pdf :
	cd figure1 && $(MAKE)

figure2/figure2.pdf : figure1/figure1.pdf
	cd figure2 && $(MAKE)

figure3/figure3.pdf : figure2/figure2.pdf
	cd figure3 && $(MAKE)

figure4/figure4.pdf : figure3/figure3.pdf
	cd figure4 && $(MAKE)

figure5/figure5.pdf : figure4/figure4.pdf
	cd figure5 && $(MAKE)

figure6/figure6.pdf : figure5/figure5.pdf
	cd figure6 && $(MAKE)
