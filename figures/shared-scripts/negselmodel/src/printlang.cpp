/** Prints all strings in the language accepted by the input DFA.

	Assumes that the input DFA contains no self-loops (in other words, it accepts only
	finite strings).

    Usage: printlang [initial]

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

#include "aminoacids.hpp"

using namespace std;

map< int, vector<int> > edges;
map< int, vector<int> > labels;
vector< int > prefix;
set<int> accepting;


void dfs( int start ){
	if( accepting.count( start )>0 ){
		for( int i = 0 ; i < prefix.size() ; i ++ ){
			cout << aminoacids[prefix[i]-1];
		}
		cout << endl;
	}
	if( edges.count( start ) > 0 ){
		for( int i = 0 ; i < edges.at(start).size() ; i ++ ){
			prefix.push_back( labels[start][i] );
			dfs( edges[start][i] );
			prefix.pop_back();
		}
	}
}

int main( int argc, char * argv[] )
{
	string line; int tok;
	vector<int> tokens;
	int initial=0;
	if( argc == 2 ){
		initial = atoi(argv[1]);
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
			//cout << "acc " << tokens[0] << endl;
			accepting.insert( tokens[0] );
		}
		if( tokens.size() >= 3 ){
			//cout << "edg " << tokens[0] << " " << tokens[1] << endl;
			edges[tokens[0]].push_back( tokens[1] );
			labels[tokens[0]].push_back( tokens[2] );
		}
	}
	dfs( initial );
}
