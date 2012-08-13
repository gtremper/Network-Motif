import networkx as nx
from random import *
import time


#input: a directed graph and the number of swaps you want to make
#output: a copy of the swapped graph.		
def randomize_graph(graph,numpasses):
    G = graph.copy()
    for i in xrange(numpasses):
        success = False
        while not success:
            edges = G.edges()
            edgeSet = set(edges)
            edge1 = choice(edges)
            a,b = edge1
            shuffle(edges)
            for edge2 in edges:
                c,d = edge2
                if (a,d) not in edgeSet and (c,b) not in edgeSet:
                    success = True
                    break
        G.add_edge(a,d)
        G.add_edge(c,b)
        G.remove_edge(a,b)
        G.remove_edge(c,d)
    return G
	
