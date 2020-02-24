all : contiguous-fa-lang contiguous-fa countpaths mjc-fa pattern-fa contiguous-negative-selection-lang pattern-exact-fa makerep-contiguous-fa makerep-contiguous-fa-lang makerep-pattern-exact-fa contiguous-negative-selection #compact-stringset-representation

CXXFLAGS=-O3 -std=c++0x
FSTLIB=-lfst
FSTFLAGS=-L /usr/local/Cellar/openfst/1.6.3/lib/ -I /usr/local/Cellar/openfst/1.6.3/include/ 

%-fa : %-fa.cpp dfa.hpp idfa.hpp refinable-partition.hpp aminoacids.hpp miyazawa-jernigan-matrix.hpp
	g++ $(CXXFLAGS) -o $*-fa $*-fa.cpp

printlang: printlang.cpp
	g++ -O3 -o $@ $<

countpaths: countpaths.cpp
	g++ -O3 -o $@ $<

contiguous-negative-selection-lang : contiguous-negative-selection-lang.cpp latinalphabet.hpp
	g++ $(FSTFLAGS) $(CXXFLAGS) -o $@ $< -lfst

compact-stringset-representation : compact-stringset-representation.cpp aminoacids.hpp
	g++ $(FSTFLAGS) $(CXXFLAGS) -o $@ $< $(FSTLIB) 

contiguous-negative-selection : contiguous-negative-selection.cpp aminoacids.hpp
	g++ $(FSTFLAGS) $(CXXFLAGS) -o $@ $< -lfst

clean : 
	rm -f contiguous-fa-lang contiguous-fa countpaths mjc-fa pattern-fa contiguous-negative-selection-lang pattern-exact-fa makerep-contiguous-fa makerep-contiguous-fa-lang makerep-pattern-exact-fa contiguous-negative-selection
