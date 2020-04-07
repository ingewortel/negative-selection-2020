#ifndef DFA_HPP
#define DFA_HPP

#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <algorithm>

#include <assert.h>     /* assert */

#include "idfa.hpp"
#include "refinable-partition.hpp"

// "Vcode" is a hash map used to go from non-dense state encodings
// to a dense state encoding.
// e.g. in a non-dense encoding, the dfa 1 -(a)-> 10 could be specified.
//      in a dense encoding, this must be 1 -(a)-> 2.

// If an application does not wish to take care itself of dense encoding,
// it can initialize an empty map<long,long> and pass this every time 
// a DFA manipulation routine such as addState is called.

// This automaton class assumes that all stored strings are of the same length.
// This means we only need one final state.

// The automaton does nothing by itself to actually check that only equal-length
// strings are stored.

using namespace std;

template <unsigned int S> class DFA {
	public:
		DFA(){ finalstate = -1; nr_states = 0; };

		// generate a DFA recognizing only the given string.
		DFA( const vector<int> & s ){
			nr_states = s.size() + 1;
			for( int i = 0 ; i < s.size() ; i ++ ){
				EE[s[i]][(long)i]=(long)i+1L;
			}
			finalstate = s.size();
		}
		
		// generate a DFA from given IDFA
		DFA( const IDFA<int>& other ) : nr_states(2), finalstate(1){
			vector<IDFAState<int> *> P;
			unordered_map<IDFAState<int> *,long> state_nr;
			P.push_back( other.initialState() );
			state_nr[other.initialState()]=0;
			state_nr[other.finalState()]=1;
			while( !P.empty() ){
				IDFAState<int> * s = P.back();
				P.pop_back();
				for( vector<IDFATransition<int> >::iterator it = s->getTrans()-> begin();
					it != s->getTrans()-> end(); ++it ){
					IDFATransition<int> ti = *it;
					IDFAState<int> *t = ti.getTarget();
					if( state_nr.find(t) == state_nr.end() ){
						state_nr[t] = nr_states++;
						P.push_back( t );
					}
					EE[ti.getLabel()][state_nr[s]]=state_nr[t];
				}
			}
		}

		//void addEdge( long u, long v, int c );
		//void setAccepting( long a );
		pair<long,int> commonPrefix( const vector<int>& s );
		pair<long,int> commonSuffix( const vector<int>& s );
		void addAllStrings( long q, int l );
		void addSuffix( long q, const vector<int>& s, int i );
		void addString( const vector<int>& s );
		long addState( long s, map<long,long>& Vcode );
		void addEdge( long u, long v, int c, map<long,long>& Vcode );
		void setAccepting( long a, map<long,long>& Vcode );
		
		void print( char * alphabet = NULL );
		void printFi( char * alphabet = NULL );
		void printDebug();
		void printComplexity();
		inline long nextState( long u, int c ){
			unordered_map<long,long>::iterator ci = EE[c].find(u);
			if( ci != EE[c].end() ){
				return ci -> second;
			} else {
				return -1;
			}
		};
		inline long nrTransitions(){
			long r = 0;
			for( int c = 0 ; c < S ; c ++ ){
				r += EE[c].size();
			}
			return r;
		};
		inline bool isAccepting( long s ){
			return finalstate == s; //A.find(s) != A.end();
		};
		DFA join( DFA other );
		DFA minimize();
		static DFA minimize( IDFA<int> &other );
		
		inline DFA clone(){
			DFA r;
			r.finalstate = this->finalstate;
			r.nr_states = this->nr_states;
			for( int i = 0 ; i < S ; i ++ ){
				r.EE[i].insert(this->EE[i].begin(),this->EE[i].end());
			}
			return r;
		}
		
		long finalstate;

		long nr_states; // 0 - initial  1 - accepting are always in
		
		unordered_set<long> pureAccepting();

		unordered_map< long, long > EE[S];
		
	private:

		/* Helper functions for minimization algorithm */
		static void make_adjacent( int K[], int *A, int *F, const int nn, const int mm );
		static void reach( RefinablePartition& B, int *L, int & rr, int q );
		static void rem_unreachable( RefinablePartition& B, int *A, int *F, 
			int & rr, const int nn, int & mm, int T[], int H[], int *L );
		static DFA runMinimize( int mm, int nn, int *T, int *H, int *L, 
			int finalstate, bool input_sorted = true );

		int pureAcceptingRec( long state, unordered_set<long>& R, unordered_set<long>& V );
	
		inline long t2( long i, long j, long Npad ){
			return (i*(Npad)+j);
		}

		inline long ti( long t, long Npad ){
			return (int)(t/Npad);
		}

		inline long tj( long t, long Npad ){
			return t%Npad;
		}
};

template <unsigned int S>
DFA<S> DFA<S>::join( DFA<S> other ){
	//set<long>& V1 = this->V, V2 = other.V,
	if( this->nr_states == 0 || this->finalstate == -1){
		return other.clone(); 
	}
	if( other.nr_states == 0 || other.finalstate == -1 ){
		return this->clone();
	}
	unordered_set<long> PA1 = this->pureAccepting(), PA2 = other.pureAccepting();
	//map< long, map<int,long> >& E1 = this->E, E2 = other.E;

	DFA<S> R;
	map<long,long> Vcode;
	
	long N1=this->nr_states; long N2=other.nr_states;
	long Npad = N2+1;
	long Nacc = t2(N1,N2,Npad);
	
	set< long > Li, Lj;
	set< long > * current_level = &Li,
		* next_level = &Lj;
	current_level -> insert( 0 );
	while( current_level -> size() > 0 ){
		//cout << "level " << current_level -> size() << endl;
		next_level -> clear();
		for( set<long>::iterator it=current_level->begin();
			it != current_level->end() ; ++ it ){
			long i = ti(*it,Npad),j=tj(*it,Npad);
			//cout << "node " << i << " " << j << endl;
			if( i != N1 && j != N2 ){
				for( int c = 0 ; c < S ; c ++ ){
					long inext=this->nextState(i,c),
						jnext=other.nextState(j,c);
					if( inext >= 0 ){
						if( jnext >= 0 ){
 							if ( PA1.count( inext ) ){
 								jnext = N2;
 							} else if( PA2.count( jnext ) ){
 								inext = N1;
 							}
 						} else {
							jnext = N2;
						}
					} else if( jnext >= 0 ){
						inext = N1;
					}
					if( inext > 0 ){
						if( this->isAccepting(inext) || 
							other.isAccepting(jnext) ){
							R.addEdge( *it, Nacc, c, Vcode );
						} else {
							long jt = t2( inext, jnext, Npad );
							R.addEdge( *it, jt, c, Vcode );
							next_level -> insert( jt );
						}
					}
				}
			} else if( j==N2 ){
				bool inserted = false;
				for( int c = 0 ; c < S ; c ++ ){
					long v = this -> nextState( i, c );
					if( v >= 0 ){
						if( this -> isAccepting(v) ){
							R.addEdge( *it, Nacc, c, Vcode );
						} else {
							long jt = t2( v, N2, Npad );
							R.addEdge( *it, jt, c, Vcode );
							next_level -> insert( jt );
						}
					}
				}
			} else if( i==N1 ){
				for( int c = 0 ; c < S ; c ++ ){
					long v = other.nextState( j, c );
					if( v >= 0 ){
						if( other.isAccepting(v) ){
							R.addEdge( *it, Nacc, c, Vcode );
						} else {
							long jt = t2( N1, v, Npad );
							R.addEdge( *it, jt, c, Vcode );
							next_level -> insert( jt );
						}
					}
				}
			}
			Vcode.erase( *it );
		} // loop within level
		if( current_level == &Li ){
			next_level = &Li;
			current_level = &Lj;
		} else {
			next_level = &Lj;
			current_level = &Li;
		}
	} // loop across levels
	R.setAccepting( Nacc, Vcode );
	return R;
}

template <unsigned int S>
void DFA<S>::addSuffix( long st, const vector<int>& s, int start_i ){
	assert( s.size()-start_i > 0 );
	long u = st, v = st;
	int i;
	for (i = start_i ; i < s.size()-1 ; i ++ ){
		v = nr_states ++;
		EE[s[i]][u]=v;
		u = v;
	}
	EE[s[i]][u]=finalstate;
}

template <unsigned int S>
pair<long,int> DFA<S>::commonPrefix( const vector<int>& s )
{
	int i = 0;
	long q = 0; // initial state
	while( EE[s[i]].find(q) != EE[s[i]].end() ){
		q = EE[s[i]][q]; i ++;
	}
	return pair<long,int>(q,i);
}//common_prefix

template <unsigned int S>
void DFA<S>::addString( const vector<int>& s ){
	// create initial state if it didn't already exist
	assert( s.size() > 0 );
	if( nr_states == 0 ){ 
		nr_states = 2; 
		finalstate = 1;
	}
	pair<long,int> pre = commonPrefix( s );
	addSuffix( pre.first, s, pre.second );
}

template <unsigned int S>
void DFA<S>::addAllStrings( long st, int l ){
	long u = st, v = st;
	for (int i = 0 ; i < l ; i ++ ){
		v = nr_states ++;
		for( int j = 0 ; j < S ; j ++ ){
			EE[j][u]=v;
		}
		u = v;
	}
}

template <unsigned int S>
void DFA<S>::setAccepting( long a, map<long,long>& Vcode ){
	finalstate = addState( a, Vcode );
}

template <unsigned int S>
long DFA<S>::addState( long s, map<long,long>& Vcode ){
	map<long,long>::iterator it = Vcode.find(s);
	if(it==Vcode.end()){
		long v = nr_states++;
		Vcode[s]=v;
		return v;
	} else {
		return it->second;
	}
}

template <unsigned int S>
void DFA<S>::addEdge( long u, long v, int c, map<long,long>& Vcode ){
	long uc = addState( u, Vcode );
	long vc = addState( v, Vcode );
	EE[c][uc]=vc;
}

template <unsigned int S>
void DFA<S>::printDebug(){

	unordered_set<long> acc = pureAccepting();
	unordered_map<long,long> printed;
	
	for( int c = 0 ; c < S ; c ++ ){
		for( unordered_map< long, long >::const_iterator ei = 
			EE[c].begin() ; ei != EE[c].end() ; ++ei ){
			if( acc.find( ei->first ) != acc.end() ){
				if( printed.find( ei->first ) != printed.end() ){
					if( printed[ei->first] != ei -> second ){
						cerr << "DFA not minimal!" << endl;
						exit(1);
					}
				} else {
					std::cout << ei -> first << " " << ei -> second << " *\n";
					printed[ei->first]=ei->second;
				}
			}
			else {
				std::cout << ei -> first << " " << ei -> second << 
					" " << 
					(c+1)
				<< "\n";
			}
		}
	}

	std::cout << finalstate << "\n";
}

template <unsigned int S>
unordered_set<long> DFA<S>::pureAccepting(){
	unordered_set<long> R, V;
	pureAcceptingRec( 0, R, V );
	return R;
}

template <unsigned int S>
int DFA<S>::pureAcceptingRec( long state, unordered_set<long>& R, unordered_set<long>& V ){
	V.insert( state );
	if( isAccepting(state) ){
		R.insert( state );
		return 1;
	} else {
		int pureChildren = 0;
		for( int c = 0 ; c < S ; c ++ ){
			long v = nextState( state, c );
			if( v >= 0 ){
				if( V.find( v ) == V.end() ){
					pureChildren += pureAcceptingRec( v, R, V );
				} else {
					pureChildren += R.count(v );
				}
			}
		}; 
		if( pureChildren == S ){
			R.insert( state );
			return 1;
		}
		return 0;
	}
}

template <unsigned int S>
void DFA<S>::print( char * alphabet ){
	// output edges out of starting state first, for fstcompile
	for( int c = 0 ; c < S ; c ++ ){
		if( EE[c].find(0) != EE[c].end() ){
			std::cout << 0 << " " << EE[c][0] << 
				" " << 
				(alphabet != NULL ? alphabet[c] : c+1)
			<< "\n";
		}
	}

	// now output other edges
	for( int c = 0 ; c < S ; c ++ ){
		for( unordered_map< long, long >::const_iterator ei = 
			EE[c].begin() ; ei != EE[c].end() ; ++ei ){
			if( ei -> first != 0 ){
				std::cout << ei -> first << " " << ei -> second << 
					" " << 
					(alphabet != NULL ? alphabet[c] : c+1)
				<< "\n";
			}
		}
	}

	std::cout << finalstate << "\n";
}

template <unsigned int S>
void DFA<S>::printFi( char * alphabet ){

	cout << nr_states << " " << nrTransitions() << " 0 1 " << endl;

	for( int c = 0 ; c < S ; c ++ ){
		for( unordered_map< long, long >::const_iterator ei = 
			EE[c].begin() ; ei != EE[c].end() ; ++ei ){
			std::cout << ei -> first << 
				" " << 
				(alphabet != NULL ? alphabet[c] : c+1) << 
				" " << ei -> second <<
			 "\n";
		}
	}
	
	std::cout << finalstate << "\n";
}

/* Adjacent transitions */
template <unsigned int S>
void DFA<S>::make_adjacent( int K[], int *A, int *F, const int nn, const int mm ){
  int q, t;
  for( q = 0; q <= nn; ++q ){ F[q] = 0; }
  for( t = 0; t < mm; ++t ){ ++F[K[t]]; }
  for( q = 0; q < nn; ++q )F[q+1] += F[q];
  for( t = mm; t--; ){ A[--F[K[t]]] = t; }
}

/* Removal of irrelevant parts */
template <unsigned int S>
void DFA<S>::reach( RefinablePartition& B, int *L, int & rr, int q ){
  int i = B.L[q];
  if( i >= rr ){
    B.E[i] = B.E[rr]; B.L[B.E[i]] = i;
    B.E[rr] = q; B.L[q] = rr++; }
}

template <unsigned int S>
void DFA<S>::rem_unreachable( RefinablePartition& B, int *A, int *F, 
	int & rr, const int nn, int & mm, int T[], int H[], int *L ){
  make_adjacent( T, A, F, nn, mm ); int i, j;
  for( i = 0; i < rr; ++i ){
    for( j = F[B.E[i]];
         j < F[B.E[i] + 1]; ++j ){
      reach( B, L, rr, H[A[j]]); } }
  j = 0;
  for( int t = 0; t < mm; ++t ){
    if( B.L[T[t]] < rr ){
      H[j] = H[t]; L[j] = L[t];
      T[j] = T[t]; ++j; } }
  mm = j; B.P[0] = rr; rr = 0;
}

struct cmpclass {
	int *L;
	bool operator() ( int i,int j ) { return (L[i]<L[j]);}
} cmp;

template <unsigned int S>
DFA<S> DFA<S>::runMinimize( int mm, int nn, int *T, int *H, int *L, int finalstate, bool input_sorted ){
	int
		q0=0,    // initial state
		rr = 0,   // number of reached states
		ff = 1;    // number of final states
	int *A, *F, *W, *M; // temporary worksets
	
	//cerr << "enter run minimize " << mm << " " << nn << endl;
	
	RefinablePartition B( nn );  // blocks (consist of states)
	A = new int[ mm ]; F = new int[ nn+1 ];

	/* Remove states that cannot be reached
	   from the initial state, and from which
       final states cannot be reached */
	reach( B, L, rr, q0 ); rem_unreachable( B, A, F, rr, nn, mm, T, H, L );
	if( B.L[finalstate] < B.P[0] ){ reach( B, L, rr, finalstate ); }
	ff = rr; rem_unreachable( B, A, F, rr, nn, mm, H, T, L );

	/* Make initial partition */
	W = new int[ mm+1 ]; 
	M = new int[ mm+1];
	int w = 0;
	M[0] = ff;
	if( ff ){ W[w++] = 0; B.split( M, W, w ); }

	/* Make transition partition */
	RefinablePartition C( mm ); // cords (consist of transitions)
	if( mm ){
		if( !input_sorted ){
			//cerr << "sorting input" << endl;
			cmp.L=L;
			std::sort( C.E, C.E+mm, cmp );
		}
		C.z = M[0] = 0; int a = L[C.E[0]];
		for( int i = 0; i < mm; ++i ){
			int t = C.E[i];
			if( L[t] != a ){
				a = L[t]; C.P[C.z++] = i;
				C.F[C.z] = i; M[C.z] = 0;
			}
			C.S[t] = C.z; C.L[t] = i;
		}
		C.P[C.z++] = mm;
	}

	//cerr << "make adjacent" << endl;
	
	/* Split blocks and cords */
	make_adjacent( H, A, F, nn, mm );
	int b = 1, c = 0, i, j;
	while( c < C.z ){
		for( i = C.F[c]; i < C.P[c]; ++i ){
			B.mark( M, W, w, T[C.E[i]] ); 
		}
		//cerr << "splitting at " << c << endl;
		B.split( M, W, w ); ++c;
		while( b < B.z ){
			for( i = B.F[b]; i < B.P[b]; ++i ){
				for(
					j = F[B.E[i]];
					j < F[B.E[i]+1]; ++j
				){
					C.mark( M, W, w, A[j] ); 
				}
			}
			C.split( M, W, w ); ++b; 
		}
	}

	//cerr << "bringing dfa back" << endl;
	
	DFA r;
	r.nr_states = B.z;

	for( int b = 0; b < B.z; ++b ){
		if( B.F[b] == 0 ){
			r.finalstate=b;
			break;
		}
	}
    
	if( r.finalstate == 0 ){
		r.finalstate = B.S[q0];
	}
    
	/* Store the result */
	for( int t = 0 ; t < mm ; t ++ ){
		int c = L[t], tt = T[t], ht = H[t];
		if( B.L[tt] == B.F[B.S[tt]] ){
			int vt = B.S[tt];
			int vh = B.S[ht];
			// make sure initial state gets # 0
			if( vt == B.S[q0] ){
				vt = 0;
			} else if( vt == 0 ){
				vt = B.S[q0];
			}
			if( vh == B.S[q0] ){
				vh = 0;
			} else if( vh == 0 ){
				vh = B.S[q0];
			}
			r.EE[c][vt]=vh;
		}
	}
	
	/* Free temporary worksets */
	delete [] M;
	delete [] W;
	delete [] A;
	delete [] F;
	
	return r;
}

template <unsigned int S>
DFA<S> DFA<S>::minimize(){
	int 
		nn = nr_states,    // number of states
		mm = 0;    // number of transitions

	int
		*T,    // tails of transitions
		*L,    // labels of transitions
		*H;    // heads of transitions

	/* Read sizes and reserve most memory */
	for( int c = 0 ; c < S ; c ++ ){
		mm+=EE[c].size();
	}
	T = new int[ mm ]; L = new int[ mm ];
	H = new int[ mm ]; 
	/* Read transitions */
	int t = 0;
	for( int c = 0 ; c < S ; c ++ ){
		for( unordered_map<long,long>::iterator it = EE[c].begin() ;
				it != EE[c].end(); ++it ){
			H[t]=it->second;
			L[t]=c;
			T[t]=it->first;
			t++;
		}
	}
	
	DFA r = runMinimize( mm, nn, T, H, L, finalstate, true );
	
	delete [] H;
	delete [] T;
	delete [] L;

	return r;
};

template <unsigned int S>
DFA<S> DFA<S>::minimize( IDFA<int> &other ){
	int 
		nn = 2,    // number of states
		mm = 0;    // number of transitions

	int
		*T,    // tails of transitions
		*L,    // labels of transitions
		*H;    // heads of transitions

	//cerr << "enter minimize method" << endl;
	
	vector<IDFAState<int> *> P;
	unordered_map<IDFAState<int> *,long> state_nr;
	P.push_back( other.initialState() );
	state_nr[other.initialState()]=0;
	state_nr[other.finalState()]=1;
	while( !P.empty() ){
		IDFAState<int> * s = P.back();
		P.pop_back();
		mm += s -> fanOut();
		for( vector<IDFATransition<int> >::iterator it = s->getTrans()-> begin();
			it != s->getTrans()-> end(); ++it ){
			IDFATransition<int> ti = *it;
			IDFAState<int> *t = ti.getTarget();
			if( state_nr.find(t) == state_nr.end() ){
				state_nr[t] = nn++;
				P.push_back( t );
			}
		}
	}
	
	//cerr << "found " << mm << " transitions" << endl;
	T = new int[ mm ]; L = new int[ mm ];
	H = new int[ mm ];

	/* Read transitions */
	int t=0;
	P.push_back( other.initialState() );
	while( !P.empty() ){
		IDFAState<int> * s = P.back();
		P.pop_back();
		for( vector<IDFATransition<int> >::iterator it = s->getTrans()-> begin();
			it != s->getTrans()->end(); ++it ){
			T[t] = -state_nr[s];
			IDFATransition<int> ti = *it;
			L[t] = ti.getLabel();
			IDFAState<int> * tgt = ti.getTarget();
			unordered_map<IDFAState<int> *,long>::iterator sit = state_nr.find(tgt);
			if( sit->second > 0 ){
				P.push_back( tgt );
				sit->second = -sit->second;
			}
			H[t] = -sit->second;
			//cerr << T[t] << " " << H[t] << " " << L[t] << endl;
			t++;
		}
	}
	
	//cerr << "completed reading transitions" << endl;

	DFA r = runMinimize( mm, nn, T, H, L, 1, false );

	delete [] H;
	delete [] T;
	delete [] L;

	return r;
};



#endif
