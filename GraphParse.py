import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cPickle as pickle
import scipy.stats as stats
 
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
	
def findMotifs(graphs,motifSize,degree):
	motifs = {}
	subgraphs = []
	numstring ="/"+str(len(graphs))
	counter = 1
	for G in graphs:
		sortedWeights = np.sort(G,axis=None)
		threshold = sortedWeights[-len(G)*degree-1]
		
		sys.stdout.write("\r")
		sys.stdout.write("Motif Finding Progress: "+str(counter)+numstring)
		sys.stdout.write(" Threshold: "+str(threshold))
		sys.stdout.flush()
		counter+=1
		
		graph = nx.DiGraph(G>threshold)
		outputGraph(graph)
		#Jenky way to use c++ motif finder in python
		os.system("./Kavosh "+str(motifSize))
		data = np.loadtxt("result/MotifCount.txt",ndmin=2)
		sub = int(data[0][1]/data[0][2]+0.5)
		subgraphs.append(sub)
		for iD,total,percent in data:
			if iD in motifs:
				motifs[iD].append(percent)
			else:
				motifs[iD] = [percent]
	
	print '\nMotifs Done!'
	
	for key,value in motifs.iteritems():
		numZero = len(graphs)-len(value)
		value.extend([0 for derp in xrange(numZero)])
		motifs[key] = np.array(value)
		
	return (motifs,np.array(subgraphs),)


def aveDegree(G):
	return G.size()/float(len(G))
		
	

if __name__ == '__main__':
	f = open("aznorbert_corrsd.pkl","rb")
	data = pickle.load(f)
	f.close()	
	
	motifSize = 3
	
	motifsNL = findMotifs(data[('NL','corr')], motifSize, 15)[0]
	motifsAD = findMotifs(data[('MCI','corr')], motifSize, 15)[0]
	
	allMotifs = list(set(motifsNL.keys()) & set(motifsAD.keys()))
	
	for key in allMotifs:
		tup = stats.ks_2samp(motifsNL[key],motifsAD[key])
		print str(key)+": "+str(tup)
		
		
	print "Ave, std per key NL then AD"
	for key in allMotifs:
		nl = motifsNL[key]
		ad = motifsAD[key]
		tup = stats.ks_2samp(nl,ad)
		print str(key)
		print str(nl.mean())+"  "+str(ad.mean())
		print str(nl.std())+"  "+str(ad.std())
		print str(len(nl)-np.count_nonzero(nl)) + " "+str(len(ad)-np.count_nonzero(ad))
		if tup[1]<0.01:
			convertIDToGraph(int(key),motifSize)
			
	
	
	#motifList = sorted(motifList,key=lambda x: -x[1][1])
	
	#for i in xrange(10):
	#	print str(motifList[i][0])+" "+str(motifList[i][1])
	#	convertIDToGraph(int(motifList[i][0]),motifSize)
		
	
	#for key,value in motifList:
		#print str(key)+"  "+str(value)
		#print 
	
	
	
	
	
	
"""	
def plotThresholds():
	f = open("aznorbert_corrsd.pkl","rb")
	DATA = pickle.load(f)
	f.close()
	
	for name,graphs in DATA.iteritems():
		thresholds = np.linspace(.9999,0.9,11)
		totalAveDegree = []
		print str(name)
		for thresh in thresholds:
			minDeg = 999999
			array = np.zeros(len(graphs))
			counter = 0
			for graph in graphs:
				ave = aveDegree(nx.DiGraph(graph>thresh))
				array[counter]=ave
				counter+=1
			totalAveDegree.append(array.mean())
		plt.plot(thresholds,totalAveDegree)
		plt.savefig(str(name))
"""
	

	
	
	
