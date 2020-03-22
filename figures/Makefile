.SECONDARY:
.DELETE_ON_ERROR:

all : all-figures

figs=1 2 3 4 5 6 7 S1 S2 S3 S4 S5 S6


all-figures : $(foreach f, $(figs), figure$(f)/figure$(f).pdf)

figure1/figure1.pdf :
	cd figure1 && $(MAKE)

figureS2/figureS1.pdf : figure1/figure1.pdf
	cd figureS1 && $(MAKE)

figure2/figure2.pdf : figure1/figure1.pdf
	cd figure2 && $(MAKE)

figure3/figure3.pdf : figure2/figure2.pdf
	cd figure3 && $(MAKE)

figure4/figure4.pdf : figure3/figure3.pdf
	cd figure4 && $(MAKE)

figure5/figure5.pdf : figure4/figure4.pdf
	cd figure5 && $(MAKE)
	
figureS2/figureS2.pdf : figure5/figure5.pdf
	cd figureS2 && $(MAKE)

figureS3/figureS3.pdf : figure5/figure5.pdf
	cd figureS3 && $(MAKE)

figure6/figure6.pdf : figure5/figure5.pdf
	cd figure6 && $(MAKE)

figureS4/figureS4.pdf : figure6/figure6.pdf
	cd figureS4 && $(MAKE)

figureS5/figureS5.pdf : figure6/figure6.pdf
	cd figureS5 && $(MAKE)

figureS6/figureS6.pdf : figureS3/figureS3.pdf
	cd figureS6 && $(MAKE)

figure7/figure7.pdf : figure6/figure6.pdf figureS6/figureS6.pdf
	cd figure7 && $(MAKE)


# Setup

aux/setup : aux/csetup

aux/csetup : aux/Rsetup
	@echo Compiling C++ code in shared-scripts/negselmodel/src... &&\
	cd shared-scripts/negselmodel/src && $(MAKE) && cd ../../.. && touch $@

aux/Rsetup : shared-scripts/Rsetup.R aux/Rpackages.txt
	Rscript shared-scripts/Rsetup.R aux/Rpackages.txt && touch $@

aux/Rpackages.txt :  aux/basics | aux
	@echo Searching for used R packages... && \
	grep -r "library(" shared-scripts | sed "s/library(/ /g" | awk '{print $$2}' | sed 's/[),\,]//g' | sort | uniq > $@  && \
	grep -r "require(" shared-scripts | grep -v Rsetup.R |  sed 's/require(/ /g' | awk '{print $$2}' | sed 's/[),\,]//g' | sort | uniq >> $@

aux/basics : | aux
	bash shared-scripts/check-setup.sh && touch $@

aux : 
	mkdir -p $@