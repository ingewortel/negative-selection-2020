#ifndef IDFA_HPP
#define IDFA_HPP

#include <vector>
#include <unordered_map>

#include "assert.h"

using namespace std;

template<typename T> class IDFATransition;

template<typename T> class IDFAState{
protected:
  vector<IDFATransition<T> > out_trans; // outgoing transitions
  T payload;
  // Create new state with no transitions
public:
	IDFAState(void) {
	};
	IDFAState(int n) {
	};
	void setNext(const int label, IDFAState *target);
	// Traverse a transition
	IDFAState* next(const int label);
	/*inline int getNumber(){ return nr; }
	inline int setNumber(int n){ nr=n; }*/
	inline T getPayload(){ return payload; }
	inline void setPayload( T p ){ payload = p; }
	vector<IDFATransition<T> > *getTrans() { return &out_trans; }
	inline int fanOut(void) const { return out_trans.size(); }
};

template<typename T>
class IDFATransition {
	int label;
	IDFAState<T> *target;
public:
	IDFATransition(const int l, IDFAState<T> *t) : label(l), target(t){
	}
	inline int getLabel(void) const { return label; }
	inline IDFAState<T> *getTarget(void) const { return target; }
	inline void setTarget(IDFAState<T> *t) { target = t; }
	bool operator==(const IDFATransition &t) const {
		return label == t.label && target == t.target;
	}
	bool operator!=(const IDFATransition &t) const {
		return label != t.label || target != t.target;
	}
	bool operator<(const IDFATransition &t) const {
		return (label == t.label ? target < t.target : label < t.label);
	}
};

template<typename T>
class IDFA {
protected:
	IDFAState<T> * initial_state;
	IDFAState<T> * finalstate;
	unordered_multimap<int,IDFAState<T> *> Register;
	vector<IDFAState<T> *> tree_states;
public:
	IDFA() {
		initial_state = new IDFAState<T>();
		finalstate = new IDFAState<T>();
	}
	bool deleteChildStates( IDFAState<T> *s, map<IDFAState<T> *,bool> &visited ){
		visited[s] = true;
		vector<IDFATransition<T> > *fo = s->getTrans();
		bool final_deleted = false;
		for( typename vector<IDFATransition<T> >::iterator it = fo->begin(); 
			it != fo->end(); ++it ){
			IDFAState<T> * qn = (*it).getTarget();
			if( !visited[qn] ){
				final_deleted = deleteChildStates( qn, visited ) || final_deleted;
				visited[qn] = true;
			}
		}
		if( s == finalstate ){
			final_deleted = true;
		}
		delete s;
		return final_deleted;
	}
	~IDFA() {
		map<IDFAState<T> *,bool> v;
		if( !deleteChildStates( initial_state, v ) ){
			delete finalstate;
		}
	}
	void addSuffix( IDFAState<T> *st, const vector<int>& s, int start_i,
		int end_i, IDFAState<T> * target ){
		/*if( target != finalstate ){
			cerr << s.size() << " " << start_i << " < :: " << endl;
		}*/
		assert( s.size()-start_i > 0 );
		IDFAState<T> *u = st, *v;
		int i;
		for (i = start_i ; i < end_i ; i ++ ){
			v = new IDFAState<T>();
			u -> setNext( s[i], v );
			u = v;
		}
		u -> setNext( s[i], target );
	}
	void addTree( const vector<int>&s, int suffix_length, int nr_transitions ){
		assert( s.size() > 0 );
		IDFAState<T> * t = tree_states.size() == 0 ? finalstate : tree_states.back();
		while( suffix_length > tree_states.size() ){
			IDFAState<T> * n = new IDFAState<T>();
			tree_states.push_back( n );
			for( int i = 0 ; i < nr_transitions ; i ++ ){
				n->setNext( i, t );
			}
			t = n;
		}
		t = tree_states[suffix_length-1];
		if( s.size()==suffix_length ){
			tree_states.pop_back();
			delete t;
			if( tree_states.size() == 0 ){
				t = finalstate;
			} else {
				t = tree_states.back();
			}
			for( int i = 0 ; i < nr_transitions ; i ++ ){
				initial_state->setNext( i, t );
			}
		} else {
			//cerr << s[0] << " : " << s.size() << " " << s.size()-suffix_length-1 << endl;
			addString( s, s.size()-suffix_length-1, t );
		}
	}
	void addString( const vector<int>&s, int end_i = -1, 
			IDFAState<T> * target = NULL ){
		assert( s.size() > 0 );
		if( end_i == -1 ){
			end_i = s.size()-1;
		}
		if( target == NULL ){
			target = finalstate;
		}
		pair<IDFAState<T>*,int> pre = commonPrefix( s );
		if( pre.first == finalstate ){
			// automaton already contains the string
			return;
		}
		if (pre.first->fanOut()) {
			IDFATransition<T> &lt = pre.first->getTrans()->back();
			IDFAState<T>* t = lt.getTarget();
			if( t != finalstate ){
				lt.setTarget(replOrReg(t));
			}
		}
		addSuffix( pre.first, s, pre.second, end_i, target );
	}
	bool containsString( const vector<int>&s ){
		pair<IDFAState<T>*,int> pre = commonPrefix( s );
		
	}
	pair<IDFAState<T> *,int> commonPrefix(const vector<int> &l){
		int hops = 0; IDFAState<T> *s = initial_state;
		for(IDFAState<T> *n = s; (n = s->next(l[hops])) != NULL; s = n){
			hops++;
		}
		return pair<IDFAState<T> *,int>(s, hops);
	}
	IDFAState<T> *initialState() const{
		return initial_state;
	}
	IDFAState<T> *finalState() const{
		return finalstate;
	}
	int hash(IDFAState<T> *s){
		int sum = s == finalstate;
		int i = 0;
		for( typename vector<IDFATransition<T> >::iterator it = s->getTrans()->begin();
			it != s->getTrans()->end() ; ++it ){
			IDFATransition<T> ti = *it;
			sum += ti.getLabel()*7 + ((long(ti.getTarget())>>2)*101) * (11 + 2*i++);
		}
		return sum;
	}//hash
	
	bool statesEquivalent( IDFAState<T> * p1, IDFAState<T> * p2 ){
		if( p1 == p2 ){
			return true;
		}
		if( p1->fanOut() != p2->fanOut() || p1 == finalstate || p2 == finalstate ||
			find(tree_states.begin(),tree_states.end(),p1) != tree_states.end() ||
			find(tree_states.begin(),tree_states.end(),p2) != tree_states.end() ) {
			return false;
		}
		return equal(p1->getTrans()->begin(),
			p1->getTrans()->end(),
			p2->getTrans()->begin());
	}
	
	IDFAState<T> *registerGetOrPut( IDFAState<T> * s ){
		int hs = hash(s);
		//cerr << s << " " << "hv : " << hs << endl;
		pair <typename unordered_multimap<int,IDFAState<T>* >::iterator,
			typename unordered_multimap<int,IDFAState<T>* >::iterator> ret;
		ret = Register.equal_range(hs);
		int t = 0;
		for( typename unordered_multimap<int,IDFAState<T>* >::iterator
			it=ret.first; it!=ret.second; ++it ){
			if( statesEquivalent( s, it->second ) ){
				return it -> second;
			}
			t++;
		}
		// cerr << s << " " << hs << " " << t << endl;
		Register.insert( make_pair( hs, s ) );
		return s;
		//cerr << "done" << endl;
	}

	IDFAState<T> * replOrReg(IDFAState<T> *s){
		if( s == finalstate ){
			return s;
		}
		//cerr << "checking " << s << endl;
		if (s->fanOut()) {
			IDFATransition<T> &last_trans = s->getTrans()->back();
			last_trans.setTarget(replOrReg(last_trans.getTarget()));
		}
		IDFAState<T> *s1 = registerGetOrPut(s);
		if (s1 != s) {
			delete s;
			s = s1;
		}
		return s;
	}
};

template <typename T>
void IDFAState<T>::setNext(const int label, IDFAState<T> *target){
	typename vector<IDFATransition<T> >::iterator p;
	for (p = out_trans.begin(); p != out_trans.end() && p->getLabel() < label; p++);
	assert( p == out_trans.end() || p -> getLabel() != label );
	out_trans.insert(p, IDFATransition<T>(label, target));
};

template <typename T>
IDFAState<T>* IDFAState<T>::next(const int label){
	for (typename vector<IDFATransition<T> >::iterator p = out_trans.begin();
		 p != out_trans.end(); p++) {
		if (p->getLabel() == label) {
			return p->getTarget();
		}
	}
	return NULL;
};

#endif