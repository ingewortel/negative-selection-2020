/** Counts the number of strings accepted by the input DFA.

    Usage: countpaths [initial]

    where the optional argument gives the number of the initial state.
    By default, the initial state is 0. If a negative argument is given,
    the initial state is taken to be the starting state of the first
    input edge.
 */

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <cstdlib>
#include <random>
#include "aminoacids.hpp"

using namespace std;

map< int, long > n;
map< int, vector<int> > edges;
map< int, vector<int> > labels;
map< int, vector<double> > weights;

set<int> accepting;

mt19937 generator (time(NULL) + getpid());
uniform_real_distribution<double> dis(0.0, 1.0);


long dfs( int start ){
	if( accepting.count(start)>0 ){
		return 1;
	}
	if( edges.count( start ) == 0 ){
		return 0;
	}
	if( n.count(start)==0){
		for( int i = 0 ; i < edges.at(start).size() ; i ++ ){
			weights[start].push_back( dfs( edges[start][i] ) );
			n[start] += weights[start][i];
		}
		for( int i = 0 ; i < edges.at(start).size() ; i ++ ){
			weights[start][i] /= n[start];
		}
	}
	return n[start];
}

void sample( int start ){
	if( accepting.count(start) > 0 || edges.count(start)==0 ){
		cout << endl;
		return;
	}
	double r = dis( generator );
	for( int i = 0 ; i < edges.at(start).size() ; i ++ ){
		r -= weights[start][i];
		if( r < 0 ){
			cout << aminoacids[labels[start][i]-1];
			sample( edges[start][i] );
			return;
		}
	}
}


int main( int argc, char * argv[] )
{
	string line; int tok;
	vector<int> tokens;
	int N;

	int initial=0;
	if( argc >= 2 ){
		N = atoi(argv[1]);
		if( N < 0 ){
			N = 1;
		}
	} else {
		N = 1;
	}
	if( argc >= 3 ){
		initial = atoi(argv[2]);
	}
	while( cin ){
		getline( cin, line );
		tokens.clear();
		stringstream ss(line);
		while( ss >> tok ){
			tokens.push_back( tok );
		}
		if( tokens.size() >= 1 ){
			if( initial < 0 ){
				initial = tokens[0];
			}
		}
		if( tokens.size() == 1 ){
			accepting.insert( tokens[0] );
		}
		if( tokens.size() >= 2 ){
			edges[tokens[0]].push_back( tokens[1] );
			labels[tokens[0]].push_back( tokens[2] );
		}
	}
	long npaths = 0;
	npaths += dfs( initial );
	while( N-->0 ) sample( initial );
}
