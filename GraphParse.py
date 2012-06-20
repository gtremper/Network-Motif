import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cPickle as pickle
import scipy.stats as stats
from collections import defaultdict

USECACHE = False

def outputGraph(graph,name="OUTPUT.txt"):
	"""Output a textfile edgelist of graph"""
	G = nx.convert_node_labels_to_integers(graph,1)
	with open("result/"+name,'wb') as f:
		f.write(str(len(G)) + '\n')
		nx.write_edgelist(G,f,data=False)

def findSingleMotif(G,motifSize):
	"""Find motifs of motifSize in G"""
	outputGraph(G)
	os.system("./Kavosh -i result/OUTPUT.txt -s "+str(motifSize))
	data = np.loadtxt("result/MotifCount.txt")
	return data	

def convertIDToGraph(id,motifSize,save=False):
	"""Plot graph with id and motifSize"""
	binary = bin(id);
	adj = np.zeros(motifSize*motifSize)
	for x in xrange(motifSize*motifSize):
		if binary[-x] == 'b':
			break
		adj[-x] = int(binary[-x])
	adj.shape = (motifSize,motifSize)
	graph = nx.to_networkx_graph(adj,create_using=nx.DiGraph())
	nx.draw_circular(graph)
	if save:
		plt.savefig("result/id-"+str(id)+"size-"+str(motifSize))
	else:
		plt.show()

def findMotifs(data,key,motifSize=3,degree=10,usetotal=False):
	"""Main finding motifs routine"""
	#Check cache
	filename = str(key)+'s'+str(int(motifSize))+'d'+str(int(degree))+str(usetotal)+".pkl"
	if os.path.exists('cache/'+filename) and USECACHE:
		print "in cache"
		with open('cache/'+filename,"rb") as f:
			return pickle.load(f)
	
	if key == "rand":
		graphs = []
		for i in xrange(100):	
			x = np.random.rand(88,88)
			x -= np.diag(np.diag(x))
			graphs.append(x)
	else:
		graphs = data[key]
	
	motifs = defaultdict(list)
	numstring ="/"+str(len(graphs))
	rejected = 0
	for index,G in enumerate(graphs):
		#Cull bad graphs
		if np.count_nonzero(G)<len(G)*degree:
			rejected += 1
			continue
			
		#calculate threshold
		sortedWeights = np.sort(G,axis=None)
		threshold = sortedWeights[-len(G)*degree-1]
		#Print progress
		sys.stdout.write("\rMotif Finding Progress: "+str(index)+numstring)
		sys.stdout.write(" Threshold: "+str(threshold))
		sys.stdout.flush()
		
		#Output graph to txt file
		graph = nx.DiGraph(G>threshold)
		graph = nx.convert_node_labels_to_integers(graph,1)
		with open('result/OUTPUT.txt','wb') as f:
			f.write(str(len(graph)) + '\n')
			nx.write_edgelist(graph,f,data=False)
		#Jenky way to use c++ motif finder in python
		os.system("./Kavosh "+str(motifSize))
		data = np.loadtxt("result/MotifCount.txt",ndmin=2)
		
		for iD,total,percent in data:
			motifs[int(iD)].append(total if usetotal else percent)
		
	print '\nMotifs Done! Graphs Rejected: '+str(rejected)
	
	#add zeros to graphs that didn't contain motifs
	for key,value in motifs.iteritems():
		numZero = len(graphs)-len(value)-rejected
		value.extend([0 for derp in xrange(numZero)])
		motifs[int(key)] = np.array(value)
	
	motifs = dict(motifs)
	#add motifs to cache
	if USECACHE:
		with open('cache/'+filename,'wb') as f:
			pickle.dump(motifs,f)
		
	return motifs

def motifOrder(data,key,orderSize=3,motifSize=3,degree=10):
	"""Sorts graphs into bins based on motif frequency ordering"""	
	if key == "rand":
		graphs = []
		for i in xrange(100):	
			x = np.random.rand(88,88)
			x -= np.diag(np.diag(x))
			graphs.append(x)
	else:
		graphs = data[key]
	
	pattern = defaultdict(float)
	for G in graphs:
		#calculate threshold
		sortedWeights = np.sort(G,axis=None)
		threshold = sortedWeights[-len(G)*degree-1]
		#Output graph to txt file
		graph = nx.DiGraph(G>threshold)
		graph = nx.convert_node_labels_to_integers(graph,1)
		with open('result/OUTPUT.txt','wb') as f:
			f.write(str(len(graph)) + '\n')
			nx.write_edgelist(graph,f,data=False)
		#Jenky way to use c++ motif finder in python
		os.system("./Kavosh "+str(motifSize))
		data = np.loadtxt("result/MotifCount.txt",ndmin=2)
		
		order = []
		for iD,total,percent in data:
			order.append((iD,total))
		keys = sorted(order,key=lambda x:-x[1])
		keys = [int(k[0]) for k in keys]
		pat = tuple(keys[:orderSize])
		pattern[pat] += 1/float(len(graphs))
		
		
	total = sorted(pattern.items(), key = lambda x: -x[1])
	
	for key,value in total:
		print str(key)+": " + str(value)
		

		

def motifStats(data, motifSize=3, degree=10, usetotal=False):
	"""Outputs text file with stats on the motifs in data"""
			
	filename = "result/t_test_Deg-{0}_Size-{1}.txt".format(degree,motifSize)
	with open(filename,'w') as f:
		f.write("Student's T test comparing both MCI and AD to NL.\n\n")
		for corr in ('corr',):
			title = "P-values for "+corr+" data set compared to normal patients\n"
			f.write(title)
	
			
			#motifsNL=findMotifs(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCI=findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
			motifsAD=findMotifs(data,('AD',corr), motifSize, degree, usetotal)
			
			mats = []
			for i in xrange(108):	
				x = np.random.rand(88,88)
				x -= np.diag(np.diag(x))
				mats.append(x)
			rand = {}
			rand['derp'] = mats
			motifsNL=findMotifs(rand,'derp', motifSize, degree, usetotal)
		
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys()) )
			allMotifs.sort()
			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}\n".format('MOTIF ID','MCI','AD','NORM Var','MCI Var', 'AD Var'))
			for key in allMotifs:
				tMCI, probMCI = stats.ttest_ind(motifsMCI[key], motifsNL[key])
				tAD, probAD = stats.ttest_ind(motifsAD[key], motifsNL[key])
				normVar = motifsNL[key].var()
				mciVar = motifsMCI[key].var()
				adVar = motifsAD[key].var()
				if probMCI<0.01 or probAD<0.01:
					line = "**{0:>8}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}\n"
				elif probMCI<0.1 or probAD<0.1:
					line = "*{0:>9}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}\n"
				else:
					line = "{0:>10}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}\n"
				f.write(line.format(str(int(key)), probMCI, probAD,normVar,mciVar,adVar))
			f.write("\n\n")

					
def plotMotifGraphs(data,motifSize=3,degree=10,numofmotifs=10,usetotal=False):
	"""Draws graph compairing average motif count between samples in the data"""
	for corr in ('corr',):
		#mats = []
		#for i in xrange(108):	
		#	x = np.random.rand(88,88)
		#	x -= np.diag(np.diag(x))
		#	mats.append(x)
		#rand = {}
		#rand['derp'] = mats
		#nl=findMotifs(rand,'derp', motifSize, degree, usetotal)
		
		nl=findMotifs(data,('NL',corr), motifSize, degree,usetotal)
		mci=findMotifs(data,('MCI',corr), motifSize, degree,usetotal)
		ad=findMotifs(data,('AD',corr), motifSize, degree,usetotal)
		
		motifs = nl.items()
		motifs = sorted(motifs,key=lambda x:-x[1].mean())
		keys = [int(x[0]) for x in motifs[:numofmotifs]]
		
		meansNL = []
		meansMCI = []
		meansAD = []
		stdNL = []
		stdMCI = []
		stdAD = []
		for key in keys:
			meansNL.append(nl[key].mean() if key in nl else 0.)
			stdNL.append(nl[key].std() if key in nl else 0.)
			meansMCI.append(mci[key].mean() if key in mci else 0.)
			stdMCI.append(mci[key].std() if key in mci else 0.)
			meansAD.append(ad[key].mean() if key in mci else 0.)
			stdAD.append(ad[key].std() if key in mci else 0.)
		
		ind = np.arange(numofmotifs)
		width = 0.2 

		NLplt = plt.bar(ind, meansNL, width, color='b', yerr=stdNL, ecolor='y')
		MCIplt = plt.bar(ind+width, meansMCI, width, color='y', yerr=stdMCI, ecolor='b')
		ADplt = plt.bar(ind+width+width, meansAD, width, color='g', yerr=stdAD, ecolor='r')

		plt.ylabel('Average number of motifs')
		plt.xlabel('Motif ID')
		plt.title('Motif size '+str(motifSize) +' distribution for '+corr+" with average degree "+str(degree))
		plt.xticks(ind+width+width/2., keys)
		plt.ylim(ymin=0.0)
		plt.legend( (NLplt[0], MCIplt[0], ADplt[0]), ('NL', 'MCI', 'AD') )
		plt.grid(True)
		header = 'result/TotalMotifDis-' if usetotal else 'result/PercentMotifDis-'
		plt.savefig(header+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		plt.clf()
	

if __name__ == '__main__':
	with open("aznorbert_corrsd.pkl","rb") as f:
		data = pickle.load(f)
	
	#motifStats(data,3,10)
	plotMotifGraphs(data,3,10,10)
	
	
	
	
	
	
	
	"""
	print 'Normal'
	motifOrder(data,('NL','corr'),3,3,10)
	print 'MCI'
	motifOrder(data,('MCI','corr'),3,3,10)
	print 'AD'
	motifOrder(data,('AD','corr'),3,3,10)
	print 'Random'
	motifOrder(data,"rand",3,3,10)
	"""
	
	#for i in xrange(108):
	#	x = np.random.rand(88,88)
	#	x -= np.diag(x)
		
	

	
	
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
	

	
	
	
