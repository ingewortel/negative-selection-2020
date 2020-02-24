#include <fst/fstlib.h>

#include "dfa.hpp"
#include "aminoacids.hpp"


using namespace fst;

typedef StdArc::StateId StateId;

map< StateId, long > n;



long dfs( StdVectorFst & f, StateId start ){
	if( start == f.Start() ){
		n.clear();
	}
	if( start < 0 ){
		return 0;
	}
	if( f.Final( start ) == 0 ){
		return 1;
	}
	if( n.count(start)==0){
		long r = 0;

		for (ArcIterator<StdFst> aiter(f, start); !aiter.Done(); aiter.Next()){
		    	const StdArc &arc = aiter.Value();
			r += dfs( f, arc.nextstate );
		}

		n[start] = r;
	}
	return n[start];
}


inline int c2i( int c ){
	return aa_goedel[c-'A'];
}

StdVectorFst branch2fst( DFA<N_AMINOACIDS> A ){
	// output edges out of starting state first, for fstcompile
	StdVectorFst f;
	int i;
	for( i = 0 ; i < A.nr_states ; i ++ ){
		f.AddState();
	}
	f.SetStart( 0 );
	for( int c = 0 ; c < N_AMINOACIDS ; c ++ ){
		for( unordered_map< long, long >::const_iterator ei = 
			A.EE[c].begin() ; ei != A.EE[c].end() ; ++ei ){
			f.AddArc( ei -> first, StdArc( c+1, c+1, 0, ei -> second ) );
		}
	}
	f.SetFinal( A.finalstate, 0 );
	return f;
}

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


int main( int argc, char ** argv ){

	string pep, line;
	if( argc < 4 ){
		cerr << "Usage : " << argv[0] << " [nonself] [n] [r] with FST as standard input" << endl;
		exit(1);
	}
	ifstream nonself;

	StdVectorFst *A = StdVectorFst::Read("");

	int n = atoi( argv[2] ), r = atoi( argv[3] ), i, j;

	nonself.open( argv[1] );

	//std::cout << dfs( *A, A->Start() ) << std::endl;

	DFA<N_AMINOACIDS> pepA;
	StdVectorFst pepAFst;
	StdVectorFst result;

	while( !nonself.eof() ){
		getline( nonself, line );
		if( line.length() >= n ){
			pep = line.substr( 0, n );
			pepA = branch_for_peptide( pep, n, r );
			//branch2fst( pepA ).Write("p.fst");
			pepAFst = branch2fst( pepA );
			//std::cout << dfs( pepAFst, pepAFst.Start() ) << std::endl;
			Intersect( pepAFst, *A, &result );
			//pepAFst.Write("p.fst");
			//result.SetFinal( result.NumStates()-1, 0 );
			//result.Write("r.fst");
			//std::cout << result.NumStates() << std::endl;
			//std::cout << result.Final( 6 ) << std::endl;
			std::cout << dfs( result, result.Start() ) << std::endl;
		}
	}


}
