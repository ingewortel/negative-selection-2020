# 2020-negative-selection

This repository contains the simulation code for the manuscript "Is T Cell Negative Selection a Learning Algorithm?".
It currently contains the code for the negative selection model itself (in the `model/` folder), as well as the 
sequences used in the paper (in the `data/` folder). See below for a description on how to use this code, as well as
the README.md files in each of those folders. 

All analyses are described in the methods section of the manuscript, so the content of the `model/` and `data/` folders
should together suffice to reproduce all figures in the manuscript. For convenience, we will later also upload a complete folder
with all (simulation, analysis and plotting) code that can be used to generate the figures automatically.


## How to run the code

To use the model code, first go to the `model/` folder and compile the code by typing

```
make
```

(If you do not have Make, you can also compile the code manually using the flags as shown in the `Makefile`).

Using this code, we can now run a negative selection simulation in several steps.

### Step 1 : Generate an "automaton" description of a TCR repertoire containing all possible TCRs

We first generate a TCR repertoire that contains TCRs for all possible sequences. For the languages,
this is all possible 6-mer strings that can be generated using the letters a-z and the underscore
(which we use to replace any kind of interpunction). For the peptides, this is all possible 6-mer
peptides that can be generated from the 20 amino acids.

For the language repertoire, we run:

```
makerep-contiguous-fa-lang 6 3 | fstcompile --acceptor > complete-repertoire.fst
```

Here, the '6' specifies that we look at sequences of 6 letters/amino acids. The '3' is
the value of the t parameter. The process works the same for a TCR repertoire recognizing peptides,
except that we then use `makerep-contiguous-fa` instead of `makerep-contiguous-fa-lang`.



### Step 2 : Select and compress the training set

Next, we sample our training peptides. For example, to reproduce Figure 2C of the paper, we make
an empty training set for n = 0, or a training set of 500 english strings for n = 500. Example
files are included in the `example/` folder.

We then compress these for use by our program:

```
cat ../example/trainset-n0.txt | contiguous-fa-lang 6 3 | fstcompile --acceptor > trainset-n0-compressed.fst
cat ../example/trainset-n500.txt | contiguous-fa-lang 6 3 | fstcompile --acceptor > trainset-n500-compressed.fst
```

This produces a 'postively selected repertoire', so a compressed description of all TCRs recognizing one of the
sequences in the training set.

The process works the same for a TCR repertoire recognizing peptides,
except that we then use `contiguous-fa` instead of `contiguous-fa-lang`.


### Step 3 : Negatively select the repertoire

We now use the complete repertoire and the positively selected repertoire to arrive at a negatively selected repertoire.
In essence, this means that we remove all positively selected TCRs (which all recognize one of the TCRs in the trainset)
from the complete repertoire:

```
fstdifference complete-repertoire.fst trainset-n0-compressed.fst | fstminimize > repertoire-n0.fst
fstdifference complete-repertoire.fst trainset-n500-compressed.fst | fstminimize > repertoire-n500.fst
```
The process works exactly the same for a TCR repertoire recognizing peptides instead of languages.


### Step 4 : Count remaining TCRs in the selected repertoire

We can also count the number of TCRs left in the repertoire, which shows that there are fewer after selection on a non-empty training set:

```
cat repertoire-n0.fst | fstprint | countpaths
# yields 387420489

cat repertoire-n500.fst | fstprint | countpaths
# yields 359727693

```

### Step 5 : Count reacting TCRs for a test set

We can compute the number of reacting TCRs in the post-selection repertoire for a test set of "unseen" sequences that were *not*
part of the training set used for negative selection in step 2.

Examples of such test sets for Figure 2C are included in the `example/` folder.

We run:

```
contiguous-negative-selection-lang ../example/testset-english-n0.txt 6 3 < repertoire-n0.fst > test-english-n0.txt
contiguous-negative-selection-lang ../example/testset-english-n500.txt 6 3 < repertoire-n500.fst > test-english-n500.txt
contiguous-negative-selection-lang ../example/testset-xhosa-n0.txt 6 3 < repertoire-n0.fst > test-xhosa-n0.txt
contiguous-negative-selection-lang ../example/testset-xhosa-n500.txt 6 3 < repertoire-n500.fst > test-xhosa-n500.txt
```

Dividing these numbers by the total repertoire size obtained in step 4 and multiplying by 1 million yields Figure 2C.


