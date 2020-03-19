.SECONDARY:
.DELETE_ON_ERROR:


MODELDIR=../shared-scripts/negselmodel/src
lang_selfsizes=0 1 5 10 50 75 100 150 200 250 300 400 500 600 700 800
NSIM=30
NTEST=50
trainset=../input-data/languages/en_t.txt
languages=xh en en2 me la pd ta hi
rvec=1 2 3 4 5 6

nsim_control=100 #for panelF, computing the english-vs-english control

all : figure4.pdf

# ==================== FIGURE 4 ================================================
# This code combines all panels into a single pdf for figure 4. 
figure4.pdf : latex/figure4.pdf
	@cp $< $@

latex/figure4.pdf : latex/figure4.tex ../shared-scripts/latex/figure-preamble.tex simulation-panels
	@cd latex && pdflatex figure4.tex

# These are the simulation plots to generate
panels=panelA panelB-xh panelB-self panelC-r1 panelC-r3 panelC-r6 panelD-xh panelD-me panelE panelF
simulation-panels :$(foreach p, $(panels), plots/F4$(p).pdf)



# ==================== PANEL A ================================================
# For the TCR survival plot, the workflow is slightly different than for figure3.
# Here, we need to see how the number of TCRs per string drop as the training set
# increases, which means we need the *same* test strings for simulations of
# different training set sizes. We therefore select those test sets first (one test
# set for each replicate simulation number, but the same test set is used across
# different values for ntrain). We then select training strings from all strings not
# in the test set, and generate the repertoire on those. Frequencies can then be
# analyzed in the same way as before, but using the fixed test sets chosen in step 1.

clean-panelA : 
	rm -f data/test-xh data/test-self data/test-all data/train-lang* data/repertoires/panelA &&\
	rm -f data/frequencies-panelA*
	rm -rf data/fixtest-repertoires data/testsets data/fixtest-pfout

# Step 1: Generate test sets of NTEST strings each, for self/xhosa and each simulation
# that is needed for the output panels.
data/test-xh : ../input-data/languages/xh-unseen.txt | data/testsets
	@echo Making test set for xh && \
	for sim in $$(seq 1 $(NSIM)) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk -v id="xh" '{print id, $$1}' > data/testsets/xh-sim$$sim.txt ;\
	done && touch $@

data/test-self : $(trainset) | data/testsets
	@echo Making test set for self && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk '{print "self", $$1}' > data/testsets/self-sim$$sim.txt ;\
	done && touch $@

data/test-all : data/test-xh data/test-self
	@echo merging test sets && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat data/testsets/self-sim$$sim.txt data/testsets/xh-sim$$sim.txt > data/testsets/all-sim$$sim.txt ;\
	done && touch $@


# Step 2: Generate repertoires using strings *not* in those predefined test sets.
# this is done in the same way as before, but now we make the training sets first
# (so they are not generated automatically anymore).
# First create the training sets for t=1, then just use the same ones for other t values:
data/train-lang-r1 : data/test-all | data/fixtest-repertoires/languages/trainsets
	@for sim in $$(seq 1 $(NSIM) ) ; do \
		bash ../shared-scripts/repertoires/foreign_unseen.sh data/testsets/all-sim$$sim.txt $(trainset) data/fixtest-repertoires/languages/trainsets/tmp-sim$$sim.txt && \
		for n in $(lang_selfsizes) ; do \
			bash ../shared-scripts/repertoires/train-epitopes.sh data/fixtest-repertoires/languages/trainsets/tmp-sim$$sim.txt $$n > data/fixtest-repertoires/languages/trainsets/traineps-r1-R-n$$n-sim$$sim.txt ;\
		done ;\
	done && touch $@

data/train-lang : data/train-lang-r1
	@for r in $$(seq 2 6) ; do \
		for f in $$( ls data/fixtest-repertoires/languages/trainsets/traineps-r1-R-*.txt); do \
			cp $$f $$( echo $$f | sed "s/r1/r$$r/g" ); \
		done ;\
	done && touch $@

data/fixtest-repertoires/panelA-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh data/train-lang | data/fixtest-repertoires/languages/contiguous data/fixtest-repertoires/languages/trainsets
	@bash $< $(trainset) "$(lang_selfsizes)" "3" $(MODELDIR) data/fixtest-repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

data/repertoires/panelA :  data/fixtest-repertoires/panelA-repertoires.mk data/train-lang
	@echo ".......Generating repertoires for Figure4A..." && $(MAKE) -f $< && touch $@


# Step 3: Analyze precursor frequencies for the fixed self and xh test sets.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings from step1. It uses the corresponding negatively 
# selected repertoires from the previous step.
data/frequencies-panelA.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh  data/repertoires/panelA
	@bash $< data/fixtest-repertoires/languages "$(lang_selfsizes)" "3" $(MODELDIR) \
		$(NSIM) "R" -n 6 -m "contiguous" -l "-lang" > $@

# Use the makefile above to generate the actual files with precursor frequencies.
# They will be stored in data/fixtest-pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/fixtest-pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/fixtest-pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
data/frequencies-panelA : data/frequencies-panelA.mk ../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection-lang data/repertoires/panelA
	$(MAKE) -f $< && touch $@


# Step 4: TCR survival rates are computed from the precursor frequencies computed in the 
# previous step.
# For a given value of t, merge the data of self and xh. 
data/tsurv-r%.txt : data/frequencies-panelA
	for n in $(lang_selfsizes) ; do \
		cat data/fixtest-pfout/contiguous/R/output-r$*-lang-all-n$$n-r$*.txt | awk -v n="$$n" -v r="$*" '{print $$0, n, r}' >> $@ ;\
	done && touch $@


# Step 5: plot the result.
plots/F4panelA.pdf : ../shared-scripts/plotting/plot-tsurv2.R data/tsurv-r3.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/tsurv-r3.txt $(NSIM) $(NTEST) "lang" $@


# ==================== PANEL B ================================================

# Step 1: generate the graph data (nodes and edges)
# This does not actually involve simulations, so no dependencies except folders.
data/gself-self.dot : ../shared-scripts/analysis/concgraph.py | data
	python3 $< ../input-data/languages/en-unseen.txt ../input-data/languages/en2-unseen.txt 4 10 "big" > $@
data/gself-foreign.dot : ../shared-scripts/analysis/concgraph.py | data
	python3 $< ../input-data/languages/en-unseen.txt ../input-data/languages/xh-unseen.txt 4 10 "big" > $@


# Step 2: Build graph svg graphic from .dot file using graphviz program.
# some options for graphviz:
GRAPHVIZ_OPTS=-Nlabel="" -Nstyle=filled -Npenwidth=0.0 \
	   -Nfixedsize=true -Nwidth=0.2 -Epenwidth=4\
	   -Gsize=14,14 -Grotate=90 -Goutputorder=edgesfirst -Gsplines=line
	   
plots/F4panelB-self.svg : data/gself-self.dot | plots
	neato $(GRAPHVIZ_OPTS) -Tsvg -o $@ $<

plots/F4panelB-xh.svg : data/gself-foreign.dot | plots
	neato $(GRAPHVIZ_OPTS) -Tsvg -o $@ $<
	
# Step 3: convert svg graph to pdf
plots/%.pdf : plots/%.svg
	rsvg-convert $< -f pdf -o $@

# ==================== PANEL C ================================================

# Step 0: choose strings for the example graph and put into simple txt file.
data/example-self.txt :
	@printf "%s\n%s\n%s\n" "any_th" "ady_th" "say_th" > $@
data/example-foreign.txt :
	@printf "%s\n%s\n%s\n" "uba_ub" "uba_zo" "uba_uh" > $@
	
# Step 1: generate the graph data (nodes and edges)
# This does not actually involve simulations, so no dependencies except folders.
data/ex-r%.dot : ../shared-scripts/analysis/concgraph.py data/example-self.txt data/example-foreign.txt | data
	python3 $< data/example-self.txt data/example-foreign.txt $* 1 "small" > $@

# Step 2: Build graph svg graphic from .dot file using graphviz program.
# some options for graphviz:
GRAPHVIZ_OPTS_SMALL=-Nlabel="" -Nshape=box -Npenwidth=0.0 -Nfontcolor=lightgray \
	   -Nfontname=courier -Nmargin=0.05 -Nwidth=0 -Nheight=0\
	   -Epenwidth=2\
	   -Goverlap=false -Gsize=4,4 -Goutputorder=edgesfirst -Gsplines=line

plots/F4panelC-r1.svg : data/ex-r1.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<

plots/F4panelC-r3.svg : data/ex-r3.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<

plots/F4panelC-r6.svg : data/ex-r6.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<
	
# Step 3: convert svg graph to pdf.
# This is done automatically using the same rule as in panel B. 

# ==================== PANEL D ================================================

# Step 1: Make files of 6mers to compute concordance on. This is just the input 
# strings used for all the figures, so this does not involve any simulation of 
# negative selection or any sampling of "training strings".
# In this panel, we only need xh and me.
data/self-me-6mers.txt : $(trainset) ../input-data/languages/me-unseen.txt
	@cat $< | awk '{print "self", $$1}' > $@;\
	cat ../input-data/languages/me-unseen.txt | awk '{print "me", $$1}' >> $@

data/self-xh-6mers.txt : $(trainset) ../input-data/languages/xh-unseen.txt
	@cat $< | awk '{print "self", $$1}' > $@ ;\
	cat ../input-data/languages/xh-unseen.txt | awk '{print "xh", $$1}' >> $@

# Step 2: Compute xh/me concordances for different values of the threshold t 
# (here called r for historical reasons).
data/concordances-xh-r%.txt : ../shared-scripts/analysis/compute-concordance.R data/self-xh-6mers.txt ../shared-scripts/analysis/concordance-functions2.R
	@echo "....Computing all xh concordances at t = $*" && Rscript $< data/self-xh-6mers.txt $* 1 $@

data/concordances-me-r%.txt : ../shared-scripts/analysis/compute-concordance.R data/self-me-6mers.txt ../shared-scripts/analysis/concordance-functions2.R
	@echo "....Computing all me concordances at t = $*" && Rscript $< data/self-me-6mers.txt $* 1 $@

# Step 3: Combine the files of different values of t into a single file, using the
# files data/concordances-xh-r[t-value].txt generated in step 2. 
data/concordances-xh.txt : $(foreach r, $(rvec), data/concordances-xh-r$(r).txt)
	@for r in $(rvec); do \
		conc=$$(cat data/concordances-xh-r$$r.txt | grep xh | awk '{print $$2}');\
		echo $$r $$conc >> $@ ;\
	done
data/concordances-me.txt : $(foreach r, $(rvec), data/concordances-me-r$(r).txt)
	@for r in $(rvec) ; do\
		conc=$$(cat data/concordances-me-r$$r.txt | grep me | awk '{print $$2}');\
		echo $$r $$conc >> $@;\
	done

# Step 4: Plot the result.	
plots/F4panelD-xh.pdf : ../shared-scripts/plotting/plot-concordance-rvalue.R data/concordances-xh.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/concordances-xh.txt $@
plots/F4panelD-me.pdf : ../shared-scripts/plotting/plot-concordance-rvalue.R data/concordances-me.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/concordances-me.txt $@


# ==================== PANEL E ====================================================

# Step 1: Concordances for Xhosa at different values for t we already have from panel D,
# in the file data/concordances-xh.txt.

# Step 2: Now we also need the outcome from negative selection simulations with an
# english vs xhosa testset, using ntrain=800 strings. 
# First, we need the repertoires, which we make in the same way as in figure 3 (see
# the makefile there for details).
data/repertoires/panelE-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets
	@bash $< $(trainset) "800" "$(rvec)" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

data/repertoires/panelE :  data/repertoires/panelE-repertoires.mk
	@echo ".......Generating repertoires for Figure4E..." && $(MAKE) -f $< && touch $@

# Step 3: Compute the precursor frequencies for self and xh on these repertoires.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
data/frequencies-panelE-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	data/repertoires/panelE
	@rm -f data/pfout/contiguous/R/output-r*-lang-self-n800-*.txt;\
	bash $< data/repertoires/languages "800" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

data/frequencies-panelE-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt data/repertoires/panelE 
	@rm -f data/pfout/contiguous/R/output-r*-lang-xh-n800-*.txt;\
	bash $< data/repertoires/languages "800" "$(rvec)" \
		../input-data/languages/xh-unseen.txt "xh" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

# Use the makefiles below to generate the actual files with precursor frequencies.
# They will be stored in data/pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
data/frequencies-panelE-% : data/frequencies-panelE-%.mk data/repertoires/panelE \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

data/frequencies-panelE : data/frequencies-panelE-self data/frequencies-panelE-xh
	touch $@

# Step 4: Analyze the outputs to get discrimination scores.
data/panelE-precursors-calc.txt : data/analysis-panelE.mk ../shared-scripts/analysis/analyze-precursors3.R data/frequencies-panelE
	@echo "Analyzing repertoires for Figure 4E" && $(MAKE) -f $< $@

data/analysis-panelE.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh data/frequencies-panelE
	@bash $< data/pfout/contiguous/R "800" "$(rvec)" "xh" "$(NSIM)" "data/panelE-precursors" -l "-lang" > $@


# Step 5: Plot
plots/F4panelE.pdf : ../shared-scripts/plotting/top-vs-concordance-xh.R data/concordances-xh.txt data/panelE-precursors-calc.txt ../shared-scripts/plotting/mytheme.R ../shared-scripts/plotting/repositionlabel.R | plots
	Rscript $< data/concordances-xh.txt data/panelE-precursors-calc.txt 800 $@


# ==================== PANEL F ====================================================

# Step 0: Make a single file with 6mers from all languages to compute concordance on. 
# This is just the input strings used for all the figures, so this does not involve any 
# simulation of negative selection or any sampling of "training strings".
data/all-6mers.txt : $(trainset) $(foreach l,$(languages),../input-data/languages/$(l)-unseen.txt)
	@cat $< | awk '{print "self", $$1}' > $@ ;\
	for l in $(languages); do \
		cat ../input-data/languages/$$l-unseen.txt | awk -v l=$$l '{print l, $$1}' >> $@ ;\
	done

# Step 1: Compute concordances for all these languages at t = 3
data/concordances-r3.txt : ../shared-scripts/analysis/compute-concordance.R data/all-6mers.txt ../shared-scripts/analysis/concordance-functions2.R
	@echo "....Computing all language concordances at t = 3" && Rscript $< data/all-6mers.txt 3 1 $@

# Step 2: Now we also need the outcome from negative selection simulations with an
# english vs xhosa testset, using ntrain=800 strings. 
# First, we need the repertoires, but since making a repertoire only involves selecting
# against self-strings, we can use the same repertoires as for panel E. 
data/repertoires/panelF : data/repertoires/panelE
	@touch $@

# Step 3: Compute the precursor frequencies for self and other languages on these repertoires.
# We don't need to do self and xh, since these have been analyzed already for panelE.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
data/frequencies-panelF-%.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/%-unseen.txt data/repertoires/panelF 
	@bash $< data/repertoires/languages "800" "3" \
		../input-data/languages/$*-unseen.txt "$*" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

# Use the makefiles below to generate the actual files with precursor frequencies.
# They will be stored in data/pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
data/frequencies-panelF-% : data/frequencies-panelF-%.mk data/repertoires/panelF \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

other_languages=en en2 me la pd ta hi
data/frequencies-panelF : $(foreach l, $(other_languages), data/frequencies-panelF-$(l)) data/frequencies-panelE
	@touch $@

# Step 4: Analyze the outputs to get discrimination scores.
data/panelF-precursors-calc.txt : data/analysis-panelF.mk ../shared-scripts/analysis/analyze-precursors3.R data/frequencies-panelF
	@echo "Analyzing repertoires for Figure 4F" && $(MAKE) -f $< $@

data/analysis-panelF.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh data/frequencies-panelF
	@bash $< data/pfout/contiguous/R "800" "3" "$(languages)" "$(NSIM)" "data/panelF-precursors" -l "-lang" > $@


# Step 5: Repeat but with extra simulations for the "control" (english versus more english).
# We do extra simulations to see that this point really ends up at 0.5.
data/control-data.txt : ../shared-scripts/analysis/lang-analyze-control.R ../shared-scripts/analysis/concordance-functions2.R data/control
	Rscript $< "data/lang-concordance/out-sim" $(nsim_control) 3 $@

data/control : $(shell for i in $$(seq 1 $(nsim_control)); do echo data/lang-concordance/out-sim$$i.txt; done)
	touch $@
data/lang-concordance/out-sim%.txt : ../shared-scripts/analysis/lang-concordance-control.sh \
	../input-data/languages/samples/english-train.txt ../shared-scripts/get6mers/chunkify.py \
	../shared-scripts/repertoires/repertoires.sh ../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection-lang $(MODELDIR)/contiguous-fa-lang \
	$(MODELDIR)/makerep-contiguous-fa-lang
	@rm -f data/lang-concordance/repertoires/rep-sim$*.fst ;\
	bash $< ../input-data/languages/samples/english-train.txt $* $(MODELDIR)

# Step 6: Plot
plots/F4panelF.pdf : ../shared-scripts/plotting/top-vs-concordance.R data/concordances-r3.txt data/panelF-precursors-calc.txt data/control-data.txt ../shared-scripts/plotting/mytheme.R ../shared-scripts/plotting/repositionlabel.R | plots
	Rscript $< data/concordances-r3.txt data/panelF-precursors-calc.txt data/control-data.txt 800 3 $@


# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	mkdir -p data 
	
	
data/fixtest-repertoires : 
	mkdir -p $@

data/fixtest-repertoires/languages/contiguous :
	mkdir -p $@
	
data/fixtest-repertoires/languages/trainsets :
	mkdir -p $@
	
data/repertoires : 
	mkdir -p $@

data/repertoires/languages/contiguous :
	mkdir -p $@
	
data/repertoires/languages/trainsets :
	mkdir -p $@
	
data/testsets :
	mkdir -p $@
	
plots :
	mkdir -p plots

clean: latex-clean
	rm -rf data && rm -rf plots
