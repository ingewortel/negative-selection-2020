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
	@cd latex && pdflatex figure4.tex > texlog.log

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
	

# Step 1: Generate test sets of NTEST strings each, for self/xhosa and each simulation
# that is needed for the output panels.
progress/panelA-test-xh : ../input-data/languages/xh-unseen.txt | data/testsets progress
	@echo Making test set for xh && \
	for sim in $$(seq 1 $(NSIM)) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk -v id="xh" '{print id, $$1}' > data/testsets/xh-sim$$sim.txt ;\
	done && touch $@

progress/panelA-test-self : $(trainset) | data/testsets progress
	@echo Making test set for self && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat $< | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n $(NTEST) | awk '{print "self", $$1}' > data/testsets/self-sim$$sim.txt ;\
	done && touch $@

progress/panelA-test-all : progress/panelA-test-xh progress/panelA-test-self | progress
	@echo merging test sets && \
	for sim in $$(seq 1 $(NSIM) ) ; do \
		cat data/testsets/self-sim$$sim.txt data/testsets/xh-sim$$sim.txt > data/testsets/all-sim$$sim.txt ;\
	done && touch $@


# Step 2: Generate repertoires using strings *not* in those predefined test sets.
# this is done in the same way as before, but now we make the training sets first
# (so they are not generated automatically anymore).
# First create the training sets for t=1, then just use the same ones for other t values:
progress/panelA-train-r1 : progress/panelA-test-all | data/fixtest-repertoires/languages/trainsets progress
	@for sim in $$(seq 1 $(NSIM) ) ; do \
		bash ../shared-scripts/repertoires/foreign_unseen.sh data/testsets/all-sim$$sim.txt $(trainset) data/fixtest-repertoires/languages/trainsets/tmp-sim$$sim.txt && \
		for n in $(lang_selfsizes) ; do \
			bash ../shared-scripts/repertoires/train-epitopes.sh data/fixtest-repertoires/languages/trainsets/tmp-sim$$sim.txt $$n > data/fixtest-repertoires/languages/trainsets/traineps-r1-R-n$$n-sim$$sim.txt ;\
		done ;\
	done && touch $@

progress/panelA-train : progress/panelA-train-r1
	@for r in $$(seq 2 6) ; do \
		for f in $$( ls data/fixtest-repertoires/languages/trainsets/traineps-r1-R-*.txt); do \
			cp $$f $$( echo $$f | sed "s/r1/r$$r/g" ); \
		done ;\
	done && touch $@

data/mkfiles/panelA-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh progress/panelA-train | data/fixtest-repertoires/languages/contiguous data/fixtest-repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "$(lang_selfsizes)" "3" $(MODELDIR) data/fixtest-repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

progress/panelA-repertoires :  data/mkfiles/panelA-repertoires.mk progress/panelA-train
	@echo ".......Generating repertoires for Figure4A..." && $(MAKE) -f $< && touch $@


# Step 3: Analyze precursor frequencies for the fixed self and xh test sets.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings from step1. It uses the corresponding negatively 
# selected repertoires from the previous step.
data/mkfiles/frequencies-panelA.mk : ../shared-scripts/analysis/loop-precursor-counts-fixtest-makeout.sh progress/panelA-repertoires | data/mkfiles
	@bash $< data/fixtest-repertoires/languages "$(lang_selfsizes)" "3" $(MODELDIR) \
		$(NSIM) "R" -n 6 -m "contiguous" -l "-lang" > $@

# Use the makefile above to generate the actual files with precursor frequencies.
# They will be stored in data/fixtest-pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/fixtest-pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/fixtest-pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
progress/panelA-frequencies : data/mkfiles/frequencies-panelA.mk ../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/analysis/precursor-counts-persim-fixtest.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection-lang progress/panelA-repertoires
	$(MAKE) -f $< && touch $@

# Step 4: TCR survival rates are computed from the precursor frequencies computed in the 
# previous step.
# For a given value of t, merge the data of self and xh. 
data/tsurv-r%.txt : progress/panelA-frequencies
	for n in $(lang_selfsizes) ; do \
		cat data/fixtest-pfout/contiguous/R/output-r$*-lang-all-n$$n-r$*.txt | awk -v n="$$n" -v r="$*" '{print $$0, n, r}' >> $@ ;\
	done && touch $@


# Step 5: plot the result.
plots/F4panelA.pdf : ../shared-scripts/plotting/plot-tsurv2.R data/tsurv-r3.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/tsurv-r3.txt $(NSIM) $(NTEST) "lang" $@

clean-panelA :
	rm -f figure4.pdf latex/figure4.pdf progress/panelA* data/mkfiles/*panelA* data/testsets/xh-sim*.txt \
		data/testsets/self-sim*.txt data/testsets/all-sim*.txt data/tsurv-r3.txt plots/F4panelA.pdf && \
	rm -rf data/fixtest-repertoires

# ==================== PANEL B ================================================

# Step 1: generate the graph data (nodes and edges)
# This does not actually involve simulations, so no dependencies except folders.
data/graphs/gself-self.dot : ../shared-scripts/analysis/concgraph.py | data/graphs
	python3 $< ../input-data/languages/en-unseen.txt ../input-data/languages/en2-unseen.txt 4 10 "big" > $@
	
data/graphs/gself-foreign.dot : ../shared-scripts/analysis/concgraph.py | data/graphs
	python3 $< ../input-data/languages/en-unseen.txt ../input-data/languages/xh-unseen.txt 4 10 "big" > $@


# Step 2: Build graph svg graphic from .dot file using graphviz program.
# some options for graphviz:
GRAPHVIZ_OPTS=-Nlabel="" -Nstyle=filled -Npenwidth=0.0 \
	   -Nfixedsize=true -Nwidth=0.2 -Epenwidth=4\
	   -Gsize=14,14 -Grotate=90 -Goutputorder=edgesfirst -Gsplines=line
	   
plots/F4panelB-self.svg : data/graphs/gself-self.dot | plots
	neato $(GRAPHVIZ_OPTS) -Tsvg -o $@ $<

plots/F4panelB-xh.svg : data/graphs/gself-foreign.dot | plots
	neato $(GRAPHVIZ_OPTS) -Tsvg -o $@ $<
	
# Step 3: convert svg graph to pdf
plots/%.pdf : plots/%.svg
	rsvg-convert $< -f pdf -o $@


clean-panelB : 
	rm -f data/graphs/gself-*.dot plots/F4panelB* figure4.pdf latex/figure4.pdf

# ==================== PANEL C ================================================

# Step 0: choose strings for the example graph and put into simple txt file.
data/graphs/example-self.txt : | data/graphs
	@printf "%s\n%s\n%s\n" "any_th" "ady_th" "say_th" > $@
data/graphs/example-foreign.txt : | data/graphs
	@printf "%s\n%s\n%s\n" "uba_ub" "uba_zo" "uba_uh" > $@
	
# Step 1: generate the graph data (nodes and edges)
# This does not actually involve simulations, so no dependencies except folders.
data/graphs/ex-r%.dot : ../shared-scripts/analysis/concgraph.py data/graphs/example-self.txt data/graphs/example-foreign.txt | data/graphs
	python3 $< data/graphs/example-self.txt data/graphs/example-foreign.txt $* 1 "small" > $@

# Step 2: Build graph svg graphic from .dot file using graphviz program.
# some options for graphviz:
GRAPHVIZ_OPTS_SMALL=-Nlabel="" -Nshape=box -Npenwidth=0.0 -Nfontcolor=lightgray \
	   -Nfontname=courier -Nmargin=0.05 -Nwidth=0 -Nheight=0\
	   -Epenwidth=2\
	   -Goverlap=false -Gsize=4,4 -Goutputorder=edgesfirst -Gsplines=line

plots/F4panelC-r1.svg : data/graphs/ex-r1.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<

plots/F4panelC-r3.svg : data/graphs/ex-r3.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<

plots/F4panelC-r6.svg : data/graphs/ex-r6.dot | plots
	neato $(GRAPHVIZ_OPTS_SMALL) -Tsvg -o $@ $<
	
# Step 3: convert svg graph to pdf.
# This is done automatically using the same rule as in panel B. 


clean-panelC : 
	rm -f data/graphs/example-*.txt data/graphs/ex-r*.dot plots/F4panelC* figure4.pdf latex/figure4.pdf

# ==================== PANEL D ================================================

# Step 1: Make files of 6mers to compute concordance on. This is just the input 
# strings used for all the figures, so this does not involve any simulation of 
# negative selection or any sampling of "training strings".
# In this panel, we only need xh and me.
data/concordances/self-me-6mers.txt : $(trainset) ../input-data/languages/me-unseen.txt | data/concordances
	@cat $< | awk '{print "self", $$1}' > $@;\
	cat ../input-data/languages/me-unseen.txt | awk '{print "me", $$1}' >> $@

data/concordances/self-xh-6mers.txt : $(trainset) ../input-data/languages/xh-unseen.txt | data/concordances
	@cat $< | awk '{print "self", $$1}' > $@ ;\
	cat ../input-data/languages/xh-unseen.txt | awk '{print "xh", $$1}' >> $@

# Step 2: Compute xh/me concordances for different values of the threshold t 
# (here called r for historical reasons).
data/concordances/xh-r%.txt : ../shared-scripts/analysis/compute-concordance.R \
	data/concordances/self-xh-6mers.txt ../shared-scripts/analysis/concordance-functions2.R | data/concordances
	@echo "....Computing all xh concordances at t = $*" && Rscript $< data/concordances/self-xh-6mers.txt $* 1 $@

data/concordances/me-r%.txt : ../shared-scripts/analysis/compute-concordance.R \
	data/concordances/self-me-6mers.txt ../shared-scripts/analysis/concordance-functions2.R | data/concordances
	@echo "....Computing all me concordances at t = $*" && Rscript $< data/concordances/self-me-6mers.txt $* 1 $@

# Step 3: Combine the files of different values of t into a single file, using the
# files data/concordances-xh-r[t-value].txt generated in step 2. 
data/concordances/xh.txt : $(foreach r, $(rvec), data/concordances/xh-r$(r).txt) | data/concordances
	@for r in $(rvec); do \
		conc=$$(cat data/concordances/xh-r$$r.txt | grep xh | awk '{print $$2}');\
		echo $$r $$conc >> $@ ;\
	done
data/concordances/me.txt : $(foreach r, $(rvec), data/concordances/me-r$(r).txt) | data/concordances
	@for r in $(rvec) ; do\
		conc=$$(cat data/concordances/me-r$$r.txt | grep me | awk '{print $$2}');\
		echo $$r $$conc >> $@;\
	done

# Step 4: Plot the result.	
plots/F4panelD-xh.pdf : ../shared-scripts/plotting/plot-concordance-rvalue.R data/concordances/xh.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/concordances/xh.txt $@
plots/F4panelD-me.pdf : ../shared-scripts/plotting/plot-concordance-rvalue.R data/concordances/me.txt ../shared-scripts/plotting/mytheme.R | plots
	Rscript $< data/concordances/me.txt $@


clean-panelD : 
	rm -f figure4.pdf latex/figure4.pdf plots/F4panelD* data/concordances/me*.txt \
		data/concordances/xh*.txt data/concordances/self-me-6mers.txt data/concordances/self-xh-6mers.txt

# ==================== PANEL E ====================================================

# Step 1: Concordances for Xhosa at different values for t we already have from panel D,
# in the file data/concordances-xh.txt.

# Step 2: Now we also need the outcome from negative selection simulations with an
# english vs xhosa testset, using ntrain=800 strings. 
# First, we need the repertoires, which we make in the same way as in figure 3 (see
# the makefile there for details).
data/mkfiles/panelE-repertoires.mk : ../shared-scripts/repertoires/loop-repertoires-makeout.sh | data/repertoires/languages/contiguous data/repertoires/languages/trainsets data/mkfiles
	@bash $< $(trainset) "800" "$(rvec)" $(MODELDIR) data/repertoires/languages/contiguous -T "R" -S "1 $(NSIM)" -n 6 -m "contiguous" -l "-lang" > $@

progress/panelE-repertoires :  data/mkfiles/panelE-repertoires.mk | progress
	@echo ".......Generating repertoires for Figure4E/F..." && $(MAKE) -f $< && touch $@

# Step 3: Compute the precursor frequencies for self and xh on these repertoires.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
data/mkfiles/frequencies-panelE-self.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh $(trainset) \
	progress/panelE-repertoires | data/mkfiles
	@bash $< data/repertoires/languages "800" "$(rvec)" \
		$(trainset) "self" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

data/mkfiles/frequencies-panelE-xh.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/xh-unseen.txt progress/panelE-repertoires | data/mkfiles
	@bash $< data/repertoires/languages "800" "$(rvec)" \
		../input-data/languages/xh-unseen.txt "xh" $(MODELDIR) $(NTEST) $(NSIM) "R" \
		-n 6 -m "contiguous" -l "-lang" -u > $@

# Use the makefiles below to generate the actual files with precursor frequencies.
# They will be stored in data/pfout/contiguous/R/. In the first step, a separate file
# will be generated for each individual simulation: 
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-sim[sim].txt
# These will then be concatenated to a single file with all simulations for every
# combination of t, language, and ntrain:
#		data/pfout/contiguous/R/output-r[t value]-lang-[language]-n[ntrain]-r[t value].txt
progress/panelE-frequencies-% : data/mkfiles/frequencies-panelE-%.mk progress/panelE-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo ".......Analyzing precursor frequencies for $* ..." && $(MAKE) -f $< && touch $@

progress/panelE-frequencies : progress/panelE-frequencies-self progress/panelE-frequencies-xh | progress
	touch $@

# Step 4: Analyze the outputs to get discrimination scores.
data/panelE-precursors-calc.txt : data/mkfiles/analysis-panelE.mk ../shared-scripts/analysis/analyze-precursors3.R progress/panelE-frequencies
	@echo "Analyzing repertoires for Figure 4E" && $(MAKE) -f $< $@

data/mkfiles/analysis-panelE.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh progress/panelE-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "800" "$(rvec)" "xh" "$(NSIM)" "data/panelE-precursors" -l "-lang" > $@


# Step 5: Plot
plots/F4panelE.pdf : ../shared-scripts/plotting/top-vs-concordance-xh.R data/concordances/xh.txt data/panelE-precursors-calc.txt ../shared-scripts/plotting/mytheme.R ../shared-scripts/plotting/repositionlabel.R | plots
	Rscript $< data/concordances/xh.txt data/panelE-precursors-calc.txt 800 $@

clean-panelE : 
	rm -f figure4.pdf latex/figure4.pdf plots/F4panelE* data/mkfiles/*panelE* \
		progress/panelF-repertoires \
		progress/*panelE* data/panelE* data/pfout/contiguous/R/output-r*-lang-self-n800-r*.txt \
		data/pfout/contiguous/R/output-r*-lang-xh-n800-r*.txt \
		data/repertoires/languages/contiguous/repertoire-r*-R-n800-sim*.fst.gz



# ==================== PANEL F ====================================================

# Step 0: Make a single file with 6mers from all languages to compute concordance on. 
# This is just the input strings used for all the figures, so this does not involve any 
# simulation of negative selection or any sampling of "training strings".
data/concordances/all-6mers.txt : $(trainset) $(foreach l,$(languages),../input-data/languages/$(l)-unseen.txt) | data/concordances
	@cat $< | awk '{print "self", $$1}' > $@ ;\
	for l in $(languages); do \
		cat ../input-data/languages/$$l-unseen.txt | awk -v l=$$l '{print l, $$1}' >> $@ ;\
	done

# Step 1: Compute concordances for all these languages at t = 3
data/concordances/r3.txt : ../shared-scripts/analysis/compute-concordance.R data/concordances/all-6mers.txt ../shared-scripts/analysis/concordance-functions2.R
	@echo "...Computing all language concordances at t = 3 for Figure 4F" && Rscript $< data/concordances/all-6mers.txt 3 1 $@

# Step 2: Now we also need the outcome from negative selection simulations with an
# english vs xhosa testset, using ntrain=800 strings. 
# First, we need the repertoires, but since making a repertoire only involves selecting
# against self-strings, we can use the same repertoires as for panel E. 
progress/panelF-repertoires : progress/panelE-repertoires
	@echo "...Using repertoires from Figure 4E also for Figure 4F" && touch $@

# Step 3: Compute the precursor frequencies for self and other languages on these repertoires.
# We don't need to do self and xh, since these have been analyzed already for panelE.
# This script generates makefiles to count precursor frequencies for all simulations,
# for all test strings of the given language. It uses the corresponding negatively 
# selected repertoires from the previous step and compares them to a sample of unseen
# test strings of size NTEST (which is set on top of this page). 
data/mkfiles/frequencies-panelF-%.mk : ../shared-scripts/analysis/loop-precursor-counts-makeout.sh \
	../input-data/languages/%-unseen.txt progress/panelF-repertoires | data/mkfiles
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
progress/panelF-frequencies-% : data/mkfiles/frequencies-panelF-%.mk progress/panelF-repertoires \
	../shared-scripts/analysis/precursor-counts-persim.sh \
	../shared-scripts/analysis/precursor-frequencies-c.sh \
	../shared-scripts/negselmodel/src/contiguous-negative-selection-lang \
	../shared-scripts/negselmodel/src/countpaths
	@echo "...Analyzing precursor frequencies for $* in Figure 4F..." && $(MAKE) -f $< && touch $@

other_languages=en en2 me la pd ta hi
progress/panelF-frequencies : $(foreach l, $(other_languages), progress/panelF-frequencies-$(l)) progress/panelE-frequencies
	@touch $@

# Step 4: Analyze the outputs to get discrimination scores.
data/panelF-precursors-calc.txt : data/mkfiles/analysis-panelF.mk ../shared-scripts/analysis/analyze-precursors3.R progress/panelF-frequencies
	@echo "...Analyzing repertoires for Figure 4F" && $(MAKE) -f $< $@

data/mkfiles/analysis-panelF.mk : ../shared-scripts/analysis/analyze-precursors-makeout.sh progress/panelF-frequencies | data/mkfiles
	@bash $< data/pfout/contiguous/R "800" "3" "$(languages)" "$(NSIM)" "data/panelF-precursors" -l "-lang" > $@


# Step 5: Repeat but with extra simulations for the "control" (english versus more english).
# We do extra simulations to see that this point really ends up at 0.5.
data/lang-concordance/out-sim%.txt : ../shared-scripts/analysis/lang-concordance-control.sh \
	../input-data/languages/samples/english-train.txt ../shared-scripts/get6mers/chunkify.py \
	../shared-scripts/repertoires/repertoires.sh ../shared-scripts/analysis/precursor-frequencies-c.sh \
	$(MODELDIR)/countpaths $(MODELDIR)/contiguous-negative-selection-lang $(MODELDIR)/contiguous-fa-lang \
	$(MODELDIR)/makerep-contiguous-fa-lang progress/panelF-message
	@rm -f data/lang-concordance/repertoires/rep-sim$*.fst ;\
	bash $< ../input-data/languages/samples/english-train.txt $* $(MODELDIR)

progress/panelF-message : | progress
	@echo "...Running control simulations Figure 4F" && touch $@

progress/panelF-controlsims : $(shell for i in $$(seq 1 $(nsim_control)); do echo data/lang-concordance/out-sim$$i.txt; done) | data/concordances progress
	touch $@

data/concordances/control-data.txt : ../shared-scripts/analysis/lang-analyze-control.R \
	../shared-scripts/analysis/concordance-functions2.R progress/panelF-controlsims | data/concordances
	@echo "...Analyzing control concordances Figure 4F" && \
	Rscript $< "data/lang-concordance/out-sim" $(nsim_control) 3 $@

	
# Step 6: Plot
plots/F4panelF.pdf : ../shared-scripts/plotting/top-vs-concordance.R data/concordances/r3.txt data/panelF-precursors-calc.txt data/concordances/control-data.txt ../shared-scripts/plotting/mytheme.R ../shared-scripts/plotting/repositionlabel.R | plots
	@echo "Plotting Figure 4F." && Rscript $< data/concordances/r3.txt data/panelF-precursors-calc.txt data/concordances/control-data.txt 800 3 $@


clean-panelF : 
	rm -f figure4.pdf latex/figure4.pdf plots/F4panelF* data/mkfiles/*panelF* \
		progress/panelE-repertoires \
		progress/*panelF* data/panelF* data/concordances/control-data.txt data/concordances/r3.txt \
		data/pfout/contiguous/R/output-r3-lang-*-n800-r3.txt data/concordances/all-6mers.txt \
		data/repertoires/languages/contiguous/repertoire-r*-R-n800-sim*.fst.gz &&\
	rm -rf data/lang-concordance


# ==================== FOLDER STRUCTURE ================================================
# This code automatically generates the required folders.

# Auxiliary targets
latex-clean : | latex
	@cd latex && rm -f *.aux *.log *.pdf	
	
data : 
	@mkdir -p data 
	
data/graphs :
	@mkdir -p $@

data/concordances :
	@mkdir -p $@
	
data/fixtest-repertoires : 
	@mkdir -p $@

data/fixtest-repertoires/languages/contiguous :
	@mkdir -p $@
	
data/fixtest-repertoires/languages/trainsets :
	@mkdir -p $@
	
data/repertoires : 
	@mkdir -p $@

data/repertoires/languages/contiguous :
	@mkdir -p $@
	
data/repertoires/languages/trainsets :
	@mkdir -p $@
	
data/testsets :
	@mkdir -p $@
	
data/mkfiles :
	@mkdir -p $@
	
	
progress :
	@mkdir -p $@
	
plots :
	@mkdir -p plots

clean: latex-clean
	@rm -rf data && rm -rf plots && rm -rf progress
