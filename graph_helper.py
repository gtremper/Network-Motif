import networkx as nx
import random
import numpy as np
from itertools import repeat

def randomize_graph(G, numpasses):
	"Perfoms numpasses edge swaps in place on G"
	for i in xrange(numpasses):
		success = False
		while not success:
			edges = G.edges()
			edgeSet = set(edges)
			edge1 = random.choice(edges)
			a,b = edge1
			random.shuffle(edges)
			for edge2 in edges:
				c,d = edge2
				if (a, d) not in edgeSet and (c, b) not in edgeSet:
					success = True
					break
		G.add_edge(a, d)
		G.add_edge(c, b)
		G.remove_edge(a, b)
		G.remove_edge(c, d)

def randomize_graph_count(G, numpasses):
	"Performs numpasses edgeswaps and returns diff from original"
	original = nx.to_numpy_matrix(G)
	diff = [0]
	for i in xrange(numpasses):
		success = False
		while not success:
			edges = G.edges()
			edgeSet = set(edges)
			edge1 = random.choice(edges)
			a,b = edge1
			random.shuffle(edges)
			for edge2 in edges:
				c,d = edge2
				if (a, d) not in edgeSet and (c, b) not in edgeSet:
					success = True
					break
		G.add_edge(a, d)
		G.add_edge(c, b)
		G.remove_edge(a, b)
		G.remove_edge(c, d)
		newGraph = nx.to_numpy_matrix(G)
		diff.append( np.sum(abs(original-newGraph))/2 )
	return np.array(diff)
	
"""	
def randomize_graph(G, numpasses):
	for i in xrange(numpasses):
		edges = set(G.edges())
		edge1,edge2 = random.sample(edges,2)
		a,b = edge1
		c,d = edge2
		counter = 0
		while (a, d) in edges or (c, b) in edges:
			edge1,edge2 = random.sample(edges,2)
			a,b = edge1
			c,d = edge2
			counter += 1
			if counter > 100:
				break
		G.add_edge(a, d)
		G.add_edge(c, b)
		G.remove_edge(a, b)
		G.remove_edge(c, d)
	return G
"""