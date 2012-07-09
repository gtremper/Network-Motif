import networkx as nx
import random
import numpy as np
from itertools import repeat

def randomize_graph(G, numpasses):
	for i in xrange(numpasses):
		edges = set(G.edges())
		edge1,edge2 = random.sample(edges,2)
		a,b = edge1
		c,d = edge2
		while (a, d) in edges or (c, b) in edges:
			edge1,edge2 = random.sample(edges,2)
			a,b = edge1
			c,d = edge2
		G.add_edge(a, d)
		G.add_edge(c, b)
		G.remove_edge(a, b)
		G.remove_edge(c, d)
	return G

def randomize_graph2(G, numpasses):
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
	return G