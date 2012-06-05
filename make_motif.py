'''
Created on Jun 3, 2012

@author: jason
'''
"""
generates a graph contain motifs from a text file specified in the following format:

motif [number of motifs to generate]
[list of edges in motif]

sample:

motif 10
1 2
2 3
1 4

"""
import networkx as nx
import sys
import os
import matplotlib.pyplot as plt
from random import choice

def addRandomEdge(G):
  while True:
    A = choice(G.nodes())
    B = choice(G.nodes())
    
    if A != B and (A, B) not in G.edges():
      G.add_edge(A, B)
      return

f = open(sys.argv[1], 'r')

listofgraphs = []
edgelist = None
number = None

for line in f:
  if 'motif' in line:
    if edgelist:
      G = nx.read_edgelist(edgelist, create_using=nx.DiGraph())
      listofgraphs.append((G, int(number)))
    number = line.strip().split()[1]
    edgelist = []
  else:
    edgelist.append(line.rstrip())
f.close()

G = nx.read_edgelist(edgelist, create_using=nx.DiGraph())
listofgraphs.append((G, int(number)))                 

F = None
for (G, number) in listofgraphs:
  for i in range(number):
    if F == None:
      F = G.copy()
    else:
      F = nx.disjoint_union(F,G)

while nx.number_weakly_connected_components(F) > 1:
  addRandomEdge(F)

F = nx.convert_node_labels_to_integers(F,1)

f2 = open('output.txt', 'wb')
f2.write(str(len(F.nodes())) + '\n')
nx.write_edgelist(F, f2, data=False)
f2.close()
os.system("./Kavosh -i output.txt -r 1000 -s 3")

#nx.draw(F)
#plt.show()