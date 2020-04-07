/*
  This code is a minor variation of the DFA minimization algorithm originally downloaded 
  from http://www.cs.tut.fi/~ava/DFA_minimizer.cc
  
  The license allows use and adaptation for scientific purposes; it is reproduced below.

  Copyright Antti Valmari 2012, this commented version 2013.
  This program is from

    Antti Valmari: Fast brief practical DFA minimization,
    Information Processing Letters 112 (2012) 213-217

  You may use and adapt the program for scientific purposes at your own risk,
  but you must give credit to the original author and source. Please negotiate
  with me about uses for other purposes.

  If you do not have access to the above-mentioned publication, please see
  A. Valmari, P. Lehtinen: Efficient minimization of DFAs with partial
  transition functions, Symposium on Theoretical Aspects of Computer Science,
  2008, pp. 645-656, http://drops.dagstuhl.de/volltexte/2008/1328/
  That publication explains part of the background. However, this program is
  much further optimized.

  This program inputs a deterministic finite automaton whose transition
  function is not necessarily full, and outputs the minimal automaton
  accepting the same language. The program also contains the removal of
  irrelevant parts of the DFA.

  This program runs in O(n + m log m) time, where n is the number of states
  and m is the number of defined transitions. If the transitions are given in
  the input such that all transitions with the same label are given together,
  then the transitions with another label, and so on, then the lines
  "#include <algorithm>" and "std::sort( C.E, C.E+mm, cmp );" can be removed,
  improving the running time to O(n + m log n). These should be compared to
  the running time of Hopcroft's algorithm, which is O(nk log n), where k is
  the size of the alphabet.

  This program is also fast in practice, I believe.
*/

#ifndef REFINABLE_PARTITION_HPP
#define REFINABLE_PARTITION_HPP

/* Refinable partition */
class RefinablePartition{
public:
   int z, *E, *L, *S, *F, *P;

	RefinablePartition( int n ){
		z = bool( n );  E = new int[n];
		L = new int[n]; S = new int[n];
		F = new int[n]; P = new int[n];
		for( int i = 0; i < n; ++i ){
			E[i] = L[i] = i; S[i] = 0; 
		}
		if( z ){
			F[0] = 0; P[0] = n;
		}
	}
   
	~RefinablePartition(){
	   delete[] E;
	   delete[] L;
	   delete[] S;
	   delete[] F;
	   delete[] P;
	}

	void mark( int *M, int * W, int &w, int e ){
		int s = S[e], i = L[e], j = F[s]+M[s];
		E[i] = E[j]; L[E[i]] = i;
		E[j] = e; L[e] = j;
		if( !M[s]++ ){ W[w++] = s; }
	}

	void split( int *M, int * W, int &w  ){
		while( w ){
			int s = W[--w], j = F[s]+M[s];
			if( j == P[s] ){M[s] = 0; continue;}
			if( M[s] <= P[s]-j ){
				F[z] = F[s]; P[z] = F[s] = j;
			}
			else{
				P[z] = P[s]; F[z] = P[s] = j;
			}
			for( int i = F[z]; i < P[z]; ++i ){
				S[E[i]] = z;
			}
			M[s] = M[z++] = 0;
		}
	}
};

#endif