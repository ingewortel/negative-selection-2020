# Figure by figure code

This repository contains all figure-by-figure code for the manuscript 
"Is T Cell Negative Selection a Learning Algorithm?".


## About the model
The used model is a string-based model of thymic selection in the immune system.
It implements the techniques described in:

Johannes Textor, Katharina Dannenberg, Maciej Liskiewicz:
__A Generic Finite Automata Based Approach to Implementing Lymphocyte Repertoire Models.__
In _Proceedings of the 2014 conference on genetic and evolutionary computation (GECCO'14)_, pp. 129-137. ACM, 2014. http://dx.doi.org/10.1145/2576768.2598331

The code generates deterministic automata (DFAs) that recognize certain sets of strings of fixed length. The alphabet 
is defined in "proteins.hpp" and it is normally taken to be the 20-letter amino acid alphabet (of course, this can be 
adapted, and we also use "latinalphabet.hpp" in te current manuscript).

For example, the file "contiguous-fa.cpp" implements the so-called r-contiguous matching rule. 

## Dependencies and installation

#### 1 Install Make and other command line essentials

We here use GNU Make to generate the figures automatically using all the code.
If you get an error that you don't have make, you may still have to install it.

On Mac OS X, install it by getting the Xcode command line tools. In the terminal, type:

```
xcode-select --install
```

This will also get you other essentials that you will need.

On Linux, you can run:

```
sudo apt-get install build-essential
```

which again will install make along with other essentials (such as C++ compilers).


#### 2 Install Openfst

You will need a C++ compiler and the OpenFST library binaries installed (http://www.openfst.org). On Mac OS X, you can 
install these with homebrew using

```
brew install openfst
```

OpenFST is also part of the libfst-dev package that you can install using the APT package
manager on Linux systems:

```
sudo apt-get install libfst-dev
```

#### 3 Install Graphviz

The graph visualizations in figures 4 and 5 require Graphviz and rsvg-convert. 
On Mac OS X, install this using homebrew:

```
brew install graphviz
brew install librsvg
```

On Linux:

```
sudo apt-get install graphviz
sudo apt-get install librsvg2-bin
```


#### 4 Install R and some packages

Most of the analysis and plotting was done in R, which you can get at 
[https://cloud.r-project.org/](https://cloud.r-project.org/).

You will also need to install several R packages for all the code to run. 


Once you have R (and Xcode CLT, step 1) installed, open R from the terminal using:
```
R
```

Then type the following to install the required packages:

```
install.packages( "ggplot2" )
install.packages( "dplyr" )
install.packages( "ggbeeswarm" )
install.packages( "Matrix" )
install.packages( "grid" )
install.packages( "RColorBrewer" )
```
If R asks you to select a CRAN mirror, you can just choose 1 or another of the listed
numbers (preferably a location nearby).


#### Python and packages

You will also need Python 3. On Linux:

```
sudo apt-get install python3.6
```

On Mac OS X via homebrew:

```
brew install python
```

This should also automatically install pip, the python package manager.
We can then use pip to install the required networkx/numpy packages:

```
python3 -m pip install networkx
python3 -m pip install numpy

```

#### Compiling the model code

To use the model code, first go to the `shared-scripts/negselmodel/src` folder and compile the code by typing

```
make
```

(If you do not have Make, you can also compile the code manually using the flags as shown in the `Makefile`).

This should work, but if you run into problems because the compiler cannot find the
OpenFST installation, consider setting the `-L` and `-I` flags in the `FSTFLAGS` variable
in the Makefile. 



## How to build the figures

To build a figure, go inside a folder (eg `figure1/`) and type `make`:

```
cd figure1
make
```

This will automatically build the figure. Note that simulations will be run for this and
this can take very long, so you may want to do it inside a `screen`.

If you have multiple cores, you may try using them to run simulations in parallel; e.g.
allow make to use 4 cores by saying:

```
make -j 4
```

This should work, but if it doesn't, try `make clean` and `make` without the parallel
option.