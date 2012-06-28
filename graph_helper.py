import networkx as nx
import random

def randomize_graph(G, numpasses):
    for i in range(numpasses*len(G.edges())):
        edge1 = random.choice(G.edges())
        edge2 = random.choice(G.edges())
        a,b = edge1
        c,d = edge2
        if (a, d) in G.edges() or (c, b) in G.edges():
            continue
        G.add_edge(a, d)
        G.add_edge(c, b)
        G.remove_edge(a, b)
        G.remove_edge(c, d)
    return G

