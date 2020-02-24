#include <iostream>
#include <vector>
#include <cstdlib>

#include "dfa.hpp"

using namespace std;

char aminoacids[27];
inline int c2i( int c ){
	if( c == '_' ){
		return 0;
	}
	return c-'a'+1;
}
const int N_AMINOACIDS = 27;

DFA<N_AMINOACIDS> branch_for_peptide( string pep, int n, int r ){
	DFA<N_AMINOACIDS> A;
	map<long,long> Vcode;
	for( int i = 0 ; i <= r ; i ++ ){
		for( int j = i ; j < n ; j ++ ){
			for( int k = 0 ; k < N_AMINOACIDS ; k ++ ){
				if( i == r ){
					A.addEdge( (n+1)*i+j, (n+1)*i+j+1, c2i(aminoacids[k]), Vcode );
				} else {
					if( aminoacids[k] == pep[j] ){
						if( r-i <= n-j ){
							A.addEdge( (n+1)*i+j, (n+1)*(i+1)+j+1,
								c2i(aminoacids[k]), Vcode );
						}
					} else {
						// mismatch crash
						if( r < n-j ){
							A.addEdge( (n+1)*i+j, j+1, c2i(aminoacids[k]), Vcode );
						}
					}
				}
			}
		}
	}
	A.setAccepting( (n+1)*r+n, Vcode );
	return A;
}

vector< DFA<N_AMINOACIDS> > fold( vector< DFA<N_AMINOACIDS> > A ){
	vector< DFA<N_AMINOACIDS> > R;
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
	aminoacids[0] = '_';
	for( int k = 0 ; k < 26 ; k ++ ){
		aminoacids[k+1] = 'a'+k;
	}

	string pep, affinity;
	if( argc < 3 ){
		cerr << "Usage : " << argv[0] << " [n] [r] " << endl;
		exit(1);
	}
	int n = atoi( argv[1] ), r = atoi( argv[2] );
	long current_node = 2;
	string line;
	DFA<N_AMINOACIDS> A;
	vector< DFA<N_AMINOACIDS> > As;
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
	//A.printDebug();
}
