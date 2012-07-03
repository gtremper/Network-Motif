import networkx as nx
import random

def randomize_graph(G, numpasses):
	for i in xrange(numpasses*G.number_of_edges()):
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

