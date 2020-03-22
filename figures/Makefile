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