# import fst
import itertools
import networkx as nx
from sys import argv
from random import randint,shuffle,seed
from collections import Counter

# SOME SETTINGS
# --------------------------------------------------------------------------

seed(124)

l = 6				# string motiflength
selffile = argv[1]		# file with self strings
nonselffile = argv[2]		# file with nonself strings
r = int(argv[3]) 		# contiguous stretch length
each_cc = int(argv[4])		# plot every ...th connected component
graphtype = argv[5]		# type of graph: "small" to indicate the small graph, 
				# anything else for the big graph

all_foreign_nodes = True

# layout according to type of graph (big/small)
if graphtype == "small":
    foreignlayout = " [fontcolor=dodgerblue3]"
    selflayout = " [fontcolor=dimgray]"
else:
    foreignlayout = " [shape=square color=dodgerblue3]"
    selflayout = " [shape=circle color=gray]"

foreignline = " [color=dodgerblue3]"
discordantline = " [color=red2]"
selfline = " [color=gray]"



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
    return r	# r will be a sequence of indices along the length of txt.


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

# this loop adds a node for every string x in txt
# it will build the list 'store', with an element for
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


# SELECT COMPONENTS FOR VISUALIZATION
# --------------------------------------------------------------------------

# Go over all connected components in G. Store their number in tk_i,
# and the number of nodes in tk_l.
tk_i = []
tk_l = []

ccon=nx.connected_components( G )

for i,C in enumerate(ccon):
    tk_l.append( len(C) )
    tk_i.append( i )

# Sort connected components by their number of nodes, then select every
# [each_cc]th component for visualization.
tk_i=[x for _,x in sorted( zip( tk_l, tk_i ) )]
tk_i=set([tk_i[x] for x in range(0,len(tk_i),each_cc)])


			#print tk_i

# based on the numbers of the connected components, construct a new graph
# with only those connected components
tk_c = set()

for i,C in enumerate(nx.connected_components( G )):
    if i in tk_i:
			#        print C
        tk_c.update(C)

			#print tk_c

G = G.subgraph( tk_c )

edges = G.edges()


# WRITE DOT FILE
# --------------------------------------------------------------------------
print( "graph G { " )

# function to determine edge layout based on the nodes
# it connects (concordant/discordant)
def ecol( e ):
    if cls[e[1]] == 1 and cls[e[0]] == 1:
            return foreignline
    if cls[e[1]] + cls[e[0]] == 1:
            return discordantline
    return selfline

# Print every edge + its layout as determined by ecol()
for e in edges:
    print ( "v"+str(e[0])+" -- v"+str(e[1])+ ecol( e ) )


# Print every node + layout according to class
for n in G.nodes():
    if G.degree(n) >= 0:
        if cls[n] == 1:
            print ( "v"+str(n)+ foreignlayout ) #" [fontcolor=lightpink]" # [shape=square color=lightpink]"
        else:
            print ( "v"+str(n)+ selflayout )

# Give every node a label with the corresponding string (only for the small graph)
if graphtype == "small":
    for i,n in enumerate(G.nodes()):
        print( "v"+str(n)+" [label="+txt[n]+"]" )

print( " } ")


