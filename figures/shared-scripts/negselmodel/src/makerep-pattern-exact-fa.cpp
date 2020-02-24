#include <iostream>
#include <vector>
#include <cstdlib> 

#include "dfa.hpp"
#include "aminoacids.hpp"
#include "pmbec-subst.hpp"

using namespace std;

inline int c2i( int c ){
	return aa_goedel[c-'A'];
}

DFA<N_AMINOACIDS+1> branch( int n, int r ){
	DFA<N_AMINOACIDS+1> R;
	if( n == 0 ){ return R; };

	map<long,long> Vcode;
	int i, l;
	int acc = n*n+1;
	int row = 0;
	R.addState( 0, Vcode );
	for( l = 0 ; l <= r && l <= n ; l ++ ){
		for( i = 0 ; i < n-l-1 ; i ++ ){
			for( int k = 0 ; k < N_AMINOACIDS ; k ++ ){
				R.addEdge( row+i, row+i+1, c2i(aminoacids[k]), Vcode );
			}
			if( l < r ){
				R.addEdge( row+i, row+i+n-l, N_AMINOACIDS, Vcode );
			}
		}
		if(  l < n && l == r ){
			for( int k = 0 ; k < N_AMINOACIDS ; k ++ ){
				R.addEdge( row+n-l-1, acc, c2i(aminoacids[k]), Vcode );
			}
		}
		if(  l < n && l == r-1 ){
			R.addEdge( row+n-l-1, acc, N_AMINOACIDS, Vcode );
		}
		row += n-l;
	}

	R.setAccepting( acc, Vcode );

	return R;
}

int main( int argc, char * argv[] )
{
	string pep, affinity;
	if( argc < 3 ){
		cerr << "Usage : " << argv[0] << " [n] [r] " << endl;
		exit(1);
	}
	int n = atoi( argv[1] ), r = atoi( argv[2] ); 
	DFA<N_AMINOACIDS+1> A = branch( n, r );
	A = A.minimize();
	A.print();
}
