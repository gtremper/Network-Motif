import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
 
#Out put a text file representation of graph
def outputGraph(graph,name="OUTPUT.txt"):
	G = nx.convert_node_labels_to_integers(graph,1)
	f = open("result/"+name,'w')
	f.write(str(len(G)) + '\n')
	for e in G.edges():
		f.write(str(e[0])+' '+str(e[1])+'\n')
		if not G.is_directed():
			f.write(str(e[1])+' '+str(e[0])+'\n')
	f.close()

#Find motifs of motifSize in G
def findMotifs(G,motifSize):
	outputGraph(G)
	os.system("./Kavosh -i result/OUTPUT.txt -r 1000 -s "+str(motifSize))

	
def convertIDToGraph(id,motifSize):
	binary = bin(id);
	adj = np.zeros(motifSize*motifSize)
	for x in xrange(motifSize*motifSize):
		if binary[-x] == 'b':
			break
		adj[-x] = int(binary[-x])
	adj.shape = (motifSize,motifSize)
	M = np.mat(adj)
	graph = nx.from_numpy_matrix(M,create_using=nx.DiGraph())
	nx.draw(graph)
	plt.show()

if __name__ == '__main__':
	#G = nx.gn_graph(100)
	#findMotifs(G,4)
	convertIDToGraph(642,4)
	
