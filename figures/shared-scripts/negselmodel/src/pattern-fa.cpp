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

void add_match_edge( DFA<N_AMINOACIDS+1> & R , int s1, int s2, int c, 
	map<long,long> & Vcode ){
	for( int k = 0 ; k < N_AMINOACIDS ; k ++ ){
		if( (aminoacids[k] != c) && 
			pmbec_subst[aminoacids[k]-'A'][c-'A'] ){
				R.addEdge( s1, s2, c2i(aminoacids[k]), Vcode );
		}
	}
}

DFA<N_AMINOACIDS+1> branch_for_peptide( string pep, int n, int r ){
	DFA<N_AMINOACIDS+1> R;
	if( n == 0 ){ return R; };
	if( r > n ){ r = n; };
	map<long,long> Vcode;
	int i, l;
	int acc = (2*n-r)*(r+1)/2;
	// cerr <<"A " << acc << endl;
	int row = 0;
	R.addState( 0, Vcode );
	for( l = 0 ; l <= r && l <= n ; l ++ ){
		for( i = 0 ; i < n-l-1 ; i ++ ){
			// uncomment line below for fuzzy match
			// add_match_edge( R, row+i, row+i+1, pep[i+l], Vcode );
			R.addEdge( row+i, row+i+1, c2i(pep[i+l]), Vcode );
			//cerr << row+i << "\t" << row+i+1 << "\t" << pep[i+l] << endl;
			if( l < r ){
				R.addEdge( row+i, row+i+n-l, N_AMINOACIDS, Vcode );
				//cerr << row+i << "\t" << row+i+n-l << "\t#\n";
			}
		}
		if( l < n ){
			// uncomment line below for fuzzy match
			// add_match_edge( R, row+i, acc, pep[n-1], Vcode );
			R.addEdge( row+i, acc, c2i(pep[n-1]), Vcode );
			//cerr << row+i << "\t" << acc << "\t" << pep[n-1] << "\n";
		}
		if( l < r ){
			R.addEdge( row+i, acc, N_AMINOACIDS, Vcode );
			//cerr << row+i << "\t" << acc << "\t#\n";
		}
		row += i+1;
	}
	R.setAccepting( acc, Vcode );
	return R;
}

vector< DFA<N_AMINOACIDS+1> > fold( vector< DFA<N_AMINOACIDS+1> > A ){
	vector< DFA<N_AMINOACIDS+1> > R;
	int n =  A.size();
	for( int i = 0 ; i < n ; i += 2 ){
		if( i+1 < n ){
			R.push_back( A[i].join(A[i+1]).minimize() );
		} else {
			R.push_back( A[i] );
		}
	}
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
	long current_node = 2;
	string line;
	DFA<N_AMINOACIDS+1> A;
	vector< DFA<N_AMINOACIDS+1> > As;
	while( cin ){
		getline( cin, line );
		if( line.length() >= n ){
			pep = line.substr( 0, n );
			As.push_back( branch_for_peptide( pep, n, r ) ); 
		}
	}
	
	while( As.size() > 1 ){
		As = fold( As );
	}
	if( As.size() == 1 ){
		As[0].print();
	}

	//As[0].print();
}
