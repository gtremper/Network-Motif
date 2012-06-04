import networkx as nx
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
	
	

if __name__ == '__main__':
	G = nx.gn_graph(100)
	findMotifs(G,4)