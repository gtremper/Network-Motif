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
	os.system("./Kavosh -i networks/ecoli -r 100 -s "+str(motifSize))
	f = open("result/ZS.txt")
	data = np.loadtxt(f)
	f.close()
	motifList = []
	for d in data:
		motifList.append( (int(d[0]),d[1]) )
	sort = sorted(motifList,key=lambda derp:-derp[1])
	print sort
	graph = convertIDToGraph(sort[0][0],motifSize)
	nx.draw(graph)
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
	
