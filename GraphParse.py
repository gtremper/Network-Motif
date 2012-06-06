import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
 
#Out put a text file representation of graph
def outputGraph(graph,name="OUTPUT.txt"):
	G = nx.convert_node_labels_to_integers(graph,1)
	f = open("result/"+name,'wb')
	f.write(str(len(G)) + '\n')
	nx.write_edgelist(G,f,data=False)
	f.close()

#Find motifs of motifSize in G
def findMotifs(G,motifSize):
	outputGraph(G)
	os.system("./Kavosh -i networks/ecoli -s "+str(motifSize))
	data = np.loadtxt("result/MotifCount.txt")
	motifList = []
	for d in data:
		motifList.append( (int(d[0]),int(d[1])) )
	sort = sorted(motifList,key=lambda derp:-derp[1])
	for d in sort:
		print d
	graph = convertIDToGraph(sort[0][0],motifSize)
	nx.draw_circular(graph)
	plt.show()

# convert Graph ID to a networkx graph	
def convertIDToGraph(id,motifSize):
	binary = bin(id);
	adj = np.zeros(motifSize*motifSize)
	for x in xrange(motifSize*motifSize):
		if binary[-x] == 'b':
			break
		adj[-x] = int(binary[-x])
	adj.shape = (motifSize,motifSize)
	graph = nx.to_networkx_graph(adj,create_using=nx.DiGraph())
	return graph
	

if __name__ == '__main__':
	G = nx.gn_graph(100)
	#nx.draw(G)
	#plt.show()
	findMotifs(G,4)
	convertIDToGraph(642,4)
	
