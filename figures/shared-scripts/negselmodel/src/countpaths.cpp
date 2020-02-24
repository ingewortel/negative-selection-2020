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

using namespace std;

map< int, long > n;
map< int, vector<int> > edges;
set<int> accepting;

long dfs( int start ){
	if( accepting.count(start)>0 ){
		return 1;
	}
	if( edges.count( start ) == 0 ){
		return 0;
	}
	if( n.count(start)==0){
		for( vector<int>::const_iterator c = edges.at(start).begin();
			c != edges.at(start).end(); ++c ){
			n[start] += dfs( *c );
		}
	}
	return n[start];
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
		if( tokens.size() >= 2 ){
			//cout << "edg " << tokens[0] << " " << tokens[1] << endl;
			edges[tokens[0]].push_back( tokens[1] );
		}
	}
	long npaths = 0;
	npaths += dfs( initial );
	cout << npaths << endl;
}
