# Input data

This folder contains the input data, the 6-mer letter/AA sequences used for the analyses
in this manuscript. See the description below and the methods section of the manuscript
for details.

## Languages
The `languages/` folder contains the 6-mer sequences used in figures 2 and 3 of the manuscript:

- en_t.txt : 6-mer sequences extracted from the book Moby Dick, which were used to sample non-overlapping training- and test sets. Sequences labelled as "self" in Figures 2/3 were derived from this set.
- en2-unseen.txt : 6-mer sequences extracted from the English translation of the bible, of which none also occur in en_t.txt. These sequences were used to generate the datapoint "English (different book)" in Figure 3F.
- me-unseen.txt : 6-mer sequences extracted from a Medieval English translation of the bible. None overlap with the strings in en_t.txt. 
- pd-unseen.txt : 6-mer sequences extracted from a Plautdietsch translation of the bible. None overlap with the strings in en_t.txt.
- la-unseen.txt : 6-mer sequences extracted from a Latin translation of the bible. None overlap with the strings in en_t.txt.
- hi-unseen.txt : 6-mer sequences extracted from a Hiligaynon translation of the bible. None overlap with the strings in en_t.txt
- ta-unseen.txt : 6-mer sequences extracted from a Tagalog translation of the bible. None overlap with the strings in en_t.txt.
- xh-unseen.txt : 6-mer sequences extracted from a Xhosa translation of the bible. None overlap with the strings in en_t.txt.

## Peptides
The `peptides/` folder contains the 6-mer sequences extracted from human/pathogenic HLA-A2 binders, which are used in Figures 4-6 of the manuscript.
See below for a brief description, and see the methods section of the manuscript for details.

- self-6mers.txt : all unique 6-mers extracted from the 263,216 HLA-A2 binders downloaded from the human proteome using NetMHCPan (see methods section).
- self1-6mers.txt : A subset of 262,000 6-mers sampled randomly from self-6mers.txt. This was used to sample the training sets for the self-other self graph in Figure 6B.
- self2-6mers.txt : The remaining 1,216 6-mers from self-6mers.txt, which do not occur in self1-6mers.txt. This was used to sample the test sets for the self-other self graph in Figure 6B.
- hiv-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the HIV proteome.
- ebola-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the ebola proteome.
- hcmv-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the HCMV proteome.
- hepb-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Hepatitis B proteome.
- hepc-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Hepatitis C proteome.
- lis-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Listeria proteome.
- mal-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Malaria proteome.
- vac-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Vaccinia proteome.
- zika-6mers.txt : All unique 6-mers extracted from the HLA-A2 binders downloaded from the Zika proteome.

The file `peptides/pathogeninfo.txt` contains an overview of the pathogens, the abbreviations used, and the ID numbers of the proteomes used to predict HLA-A2 binders.
