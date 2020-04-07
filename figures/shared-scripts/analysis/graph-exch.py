# import fst 
import itertools
import numpy as np
import networkx as nx
from sys import argv

# SOME SETTINGS
# --------------------------------------------------------------------------

l = 6				# string motiflength
selffile = argv[1]		# file with self strings
nonselffile = argv[2]		# file with nonself strings
r = int( argv[3] ) 		# contiguous stretch length

all_foreign_nodes = True


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
sself = rlines(selffile, 0)
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
    G.add_node(j)
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


#nedges = set()
#keep_nodes = set()

		#print "starting to make trees" 

# Add a node for every nonself peptide
# Also select the self peptides that are neighbors,
# Then select the subgraph of G that only contains these
# nonself peptides + their (self/nonself) neighbors.
#for i in nonself:
#    G.add_node( i )
#    T = nx.bfs_tree( G, i )
#    keep_nodes.add( i )
#    for j in T.successors( i ):
#        keep_nodes.add( j )
#        #keep_nodes.update( T.successors( j ) )
#for i in sself:
#	G.add_node( i )

#G = G.subgraph( keep_nodes )


# WRITE DEGREE FILE
# --------------------------------------------------------------------------
# output only for the self nodes
for n in G.nodes():
	#if cls[n] == 0:
#for n in sself:
	neighbors = G.neighbors(n)
	num_neighbors = len( list( neighbors ) )
	num_foreign = 0
	for nn in G.neighbors(n):
		#print(str(cls[nn]))
		if cls[nn] == 1:
			num_foreign = num_foreign + 1
	num_self = num_neighbors - num_foreign	
	print ( txt[n]+" "+str(cls[n])+" "+str(num_self)+" "+str(num_foreign)+" "+str(num_neighbors) )