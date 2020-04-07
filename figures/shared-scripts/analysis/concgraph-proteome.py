# import fst 
import itertools
import networkx as nx
from sys import argv

# SOME SETTINGS
# --------------------------------------------------------------------------

l = 6				# string motiflength
selffile = argv[1]		# file with self strings
nonselffile = argv[2]		# file with nonself strings
r = int( argv[3] ) 		# contiguous stretch length
labels = argv[4]		# "on" or "off"

all_foreign_nodes = True

# layout options
foreignline =  " [color=dodgerblue2]"
discordantline = " [color=black]"
selfline = " [color=gray]"
foreignlayout = " [shape=square width=1 fixedsize=true color=dodgerblue2 ]"
selflayout = " [shape=circle color=gray]"


# READ FILES
# --------------------------------------------------------------------------

txt = []
cls = []
nonself = []

# function of file f and class indicator c
# after running, all strings in f will have been added to txt,
# and for every string the corresponding c will have been added to cls.

def rlines( f, c ):
    global txt, cls
    r = []
    with open(f) as ff:
        for l in ff:
            r.append( len(txt) )
            txt.append(l.strip())
            cls.append(c)
    return r		# r will be a sequence of indices along the length of txt.


			#print "reading data"

# Read in the self and nonself files, assign them id 0 (self) and 1 (nonself)
rlines(selffile, 0)
nonself = rlines(nonselffile, 1)	# will contain indices of nonself strings in txt and cls.


# BUILD GRAPH
# --------------------------------------------------------------------------
store = {}
edges = set()


G = nx.Graph() # create an empty graph
			#print "building dictionaries"

# this loop builds the list 'store', with an element for
# every rmer containing the node numbers of the strings it occurs in.

for j,x in enumerate(txt):
    for i in range(0,l-r+1):	# the starting position of rmers in strings of length l
        h = str(i)+x[i:i+r]	# current starting position + rmer
        if h not in store:	# store keeps track of the rmers already seen
            store[h] = []
        store[h].append( j )	# add current node to the list element of its rmers.

		#print "creating edges"

		#print( ''.join( store.keys() ) )
# now use the list 'store' to make edges between every pair of nodes that share an rmer.
for i in store.keys():
    for pair in itertools.combinations( store[i], 2 ):
        G.add_edge(*pair)


nedges = set()
keep_nodes = set()

		#print "starting to make trees" 

# Add a node for every nonself peptide
# Also select the self peptides that are neighbors,
# Then select the subgraph of G that only contains these
# nonself peptides + their (self/nonself) neighbors.
for i in nonself:
    G.add_node( i )
    T = nx.bfs_tree( G, i )
    keep_nodes.add( i )
    for j in T.successors( i ):
        keep_nodes.add( j )
        #keep_nodes.update( T.successors( j ) )
G = G.subgraph( keep_nodes )


# WRITE DOT FILE
# --------------------------------------------------------------------------
print ( "graph G { " )

# function to determine edge layout based on the nodes
# it connects (concordant/discordant)
def ecol( e ):
    if cls[e[1]] == 1 and cls[e[0]] == 1:
        return foreignline 		#" [color=lightpink]"
    if cls[e[1]] + cls[e[0]] == 1:
        return discordantline		#" [color=black penwidth=15]"
    return selfline			#""

# Print every edge + its layout as determined by ecol()
for e in G.edges():
    print ( "v"+str(e[0])+" -- v"+str(e[1])+ ecol( e ) )

# Print every nonself node + layout
for n in nonself:
    print ( "v"+str(n)+ foreignlayout 	)	#" [shape=square color=lightpink width=5]"

# Give every node a label with the corresponding string (only for the small graph)
if labels == "on":
    for i,n in enumerate(G.nodes()):
        print ( "v"+str(n)+" [label="+txt[n]+"]" )

print ( " } " )
