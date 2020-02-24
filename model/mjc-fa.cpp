#include <iostream>
#include <vector>
#include <cstdlib>
#include <iterator>
#include <math.h>

#include "dfa.hpp"
#include "aminoacids.hpp"
#include "miyazawa-jernigan-matrix.hpp"


vector<int> pep2intseq( const string& pep ){
	vector<int> r(pep.length());
	for( int i = 0 ; i < pep.length() ; i ++ ){
		r[i]=aa_goedel[pep[i]-'A'];
	}
	return r;
}

string intseq2pep( const vector<int>& pep ){
	string r(pep.size(),' ');
	for( int i = 0 ; i < pep.size() ; i ++ ){
		r[i]=aminoacids[pep[i]];
	}
	return r;
}

vector<double> min_energy( const vector<int> & pep ){
	int n = pep.size();
	vector<double> r(n);
	for( int i = 0 ; i < n ; i ++ ){
		r[i] = mj_minima[pep[i]];
	}
	for( int i = 0 ; i < n-1 ; i ++ ){
		r[n-i-2] += r[n-i-1];
	}
	return r;
}

vector<double> max_energy( const vector<int> & pep ){
	int n = pep.size();
	vector<double> r(n);
	for( int i = 0 ; i < n ; i ++ ){
		r[i] = mj_maxima[pep[i]];
	}
	for( int i = 0 ; i < n-1 ; i ++ ){
		r[n-i-2] += r[n-i-1];
	}
	return r;
}

template<typename T>
long add_strings_near( IDFA<T> &A, DFA<N_AMINOACIDS> &Aacc, long qacc,
	const unordered_set<long> &acc,
	const vector<int>& pep,  vector<int> &pre,
	double thr_pos, const vector<double>& emin, const vector<double>& emax,
	int i ){
	if( thr_pos < emin[i] ){
		//cerr << "breaking " << thr_pos << " " << emin[i] << endl;
		return 0;
	}
	if( i == pep.size() ){
		if( thr_pos >= 0 && acc.find( qacc ) == acc.end() ){
 			//cerr << "adding single string" << endl;
			A.addString( pre );
			return 1;
		} else {
			return 0;
		}
	}
	if( acc.find( qacc ) != acc.end() ){
		/*if( i < 4 ){
			cerr << "branch off at " << i << endl;
		}*/
		return 0;
	}
	if( thr_pos >= emax[i] ){
		if( i < 4 ){
			//DFA<N_AMINOACIDS>(A).printDebug();
			//cerr << intseq2pep(pre) << " to insert complete subtree at " << i << endl;
		}
		A.addTree( pre, pep.size()-i, N_AMINOACIDS );
		return pow( N_AMINOACIDS, pep.size()-i );
	}
	long r = 0;
	for( int j = 0 ; j < N_AMINOACIDS ; j ++ ){
		pre[i] = j;
		double eij = mj_matrix[pep[i]][j];
		r += add_strings_near( A, Aacc, Aacc.nextState( qacc, j ), acc, pep, pre,
			thr_pos-eij, emin, emax, i+1 );
	}
	return r;
}

long trie_recurse( DFA<N_AMINOACIDS> &A, IDFA<int> &Apep,
		IDFAState<int> *pep,
		vector<int> &pre,
		const double thr_pos, int i ){
	//cerr << i << " " << pre.size() << endl;
	if( i == pre.size() ){
		IDFA<int> A2;
		vector<int> pre2( pre );
		long snew = add_strings_near( A2, A, 0, A.pureAccepting(),
			pre, pre2, thr_pos, max_energy(pre), min_energy(pre), 0 );
		if( snew > 0 ){
			/*cout << "1st guy: " << endl;
			A.printDebug();
			cout << "2nd guy: " << endl;
			DFA<N_AMINOACIDS>(A2).printDebug();*/
			A = A.join( DFA<N_AMINOACIDS>(A2) ).minimize();
		}
		//cerr << snew << endl;
		return snew;
	}
	long r = 0;
	float max_of_minima = 0.0;

	vector<int> dominated(N_AMINOACIDS,0);
	for( int j = 0 ; j < N_AMINOACIDS ; j ++ ){
		if( IDFAState<int> * tgt = pep->next(aa_by_min_energy[j]) ){
			pre[i] = aa_by_min_energy[j];
			if( pep -> fanOut() > 1 ){
				if( dominated[pre[i]] ){
					continue;
				}
				for( int k = 0 ; k < N_AMINOACIDS ; k ++ ){
					if( pre[i] != k && mj_domination[pre[i]][k] && ( tgt == pep->next(k) ) ){
						/*cerr << aminoacids[pre[i]] << " (" << pre[i] << ") dom. " <<
							aminoacids[k] << " (" << k << ")" << endl;*/
						dominated[k]=1;
					}
				}
			}
			r += trie_recurse( A, Apep, tgt, pre, thr_pos, i+1 );
		}
	}

	return r;
}

int main( int argc, char * argv[] )
{
	if( argc < 3 ){
		cerr << "Usage : " << argv[0] << " [n] [thr_pos] " << endl;
		exit(1);
	}
	int n = atoi( argv[1] );
	double thr_pos = atof( argv[2] )*n;

	string line;

	/** new method with trie/trie interaction */
	IDFA<int> Apep;
	DFA<N_AMINOACIDS> TCR;
	while( cin ){
		getline( cin, line );
		if( line.length() >= n ){
			Apep.addString( pep2intseq( line.substr( 0, n ) ) );
		}
	}
	vector<int> pre(n);
	trie_recurse( TCR, Apep, Apep.initialState(), pre, thr_pos, 0 );
	TCR.print();
}
