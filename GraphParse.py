import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cPickle as pickle
 
#Out put a text file representation of graph
def outputGraph(graph,name="OUTPUT.txt"):
	G = nx.convert_node_labels_to_integers(graph,1)
	f = open("result/"+name,'wb')
	f.write(str(len(G)) + '\n')
	nx.write_edgelist(G,f,data=False)
	f.close()

#Find motifs of motifSize in G
def findSingleMotif(G,motifSize):
	outputGraph(G)
	os.system("./Kavosh -i result/OUTPUT.txt -s "+str(motifSize))
	data = np.loadtxt("result/MotifCount.txt")
	return data	


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
	nx.draw_circular(graph)
	plt.show()
	
def findMotifs(graphs,motifSize,threshold):
	motifs = {}
	numstring ="/"+str(len(graphs))
	counter = 0
	for G in graphs:
		sys.stdout.write("\r")
		sys.stdout.write("Motif Finding Progress: "+str(counter)+numstring)
		sys.stdout.flush()
		counter+=1
		graph = nx.DiGraph(G>=threshold)
		outputGraph(graph)
		#Jenky way to use c++ motif finder in python
		os.system("./Kavosh "+str(motifSize))
		data = np.loadtxt("result/MotifCount.txt")
		for d in data:
			iD,size,percent = d
			if iD in motifs:
				motifs[iD].append(size)
			else:
				motifs[iD] = [size]
	
	print '\nMotifs Done!'
	
	for key,value in motifs.iteritems():
		numZero = len(graphs)-len(value)
		value.extend([0 for derp in xrange(numZero)])
		array = np.array(value)
		stats = (array.sum(),array.mean(),array.std(),numZero)
		motifs[key] = stats
	
	return motifs

		
		
	

if __name__ == '__main__':
	f = open("aznorbert_corrsd.pkl","rb")
	data = pickle.load(f)
	f.close()
	
	motifs = findMotifs(data[('NL','lcorr')], 3, 0.98)
	
	for key,value in motifs.iteritems():
		print "ID: "+str(key)+" Stats: "+str(value)

	
	
	
