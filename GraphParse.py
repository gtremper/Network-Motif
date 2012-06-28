import math
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cPickle as pickle
import scipy.stats as stats
import graph_helper as gh
import random
from collections import defaultdict
from math import floor

USECACHE = True

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
		x+=1
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
	plt.clf()

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
		#if np.count_nonzero(G)<len(G)*degree:
		#	rejected += 1
		#	continue
			
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

def motifOrder(data,key,num,orderSize=3,motifSize=3,degree=10):
	"""Sorts graphs into bins based on motif frequency ordering"""	
	if key == "rand":
		graphs = []
		for i in xrange(num):	
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
		
def genRandomGraphs(graphs, degree, num):
	mats = []
	for i in xrange(num):	
		x = np.random.rand(88,88)
		x -= np.diag(np.diag(x))
		mats.append(x)
	return mats

def findMotifs2(data,key,motifSize=3,degree=10,usetotal=False):
	"""Main finding motifs routine"""
	#Check cache
	filename = str(key)+'s'+str(int(motifSize))+'d'+str(int(degree))+str(usetotal)+"random.pkl"
	
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
		graph = gh.randomize_graph(graph, 3)
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

def motifStats2(data, motifSize=3, degree=10, usetotal=False):
	"""Outputs text file with stats on the motifs in data"""
			
	filename = "result/t_test_Deg-{0}_Size-{1}.txt".format(degree,motifSize)
	with open(filename,'w') as f:
		f.write("Student's T test comparing both AD/NL/MCI to Random.\n\n")
		for corr in ['corr']:
			title = "P-values for "+corr+" data set compared to random generated graphs\n"
			f.write(title)

			motifsNL=findMotifs(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCI=findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
			motifsAD=findMotifs(data,('AD',corr), motifSize, degree, usetotal)
			motifsNLR=findMotifs2(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCIR=findMotifs2(data,('MCI',corr), motifSize, degree, usetotal)
			motifsADR=findMotifs2(data,('AD',corr), motifSize, degree, usetotal)
		
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys())
							& set(motifsNLR.keys())
							& set(motifsADR.keys())
							& set(motifsMCIR.keys()) )
			allMotifs.sort()
			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}{6:>15}{7:>15}{8:>15}{9:>15}\n".format(
				'MOTIF ID','NL', 'MCI','AD', 'NLR Mean','MCIR Mean','ADR Mean', 'NLR Std','MCIR Std', 'ADR Std'))
			
			motifStats = []
			for key in allMotifs:
				tMCI, probMCI = stats.ttest_ind(motifsMCI[key], motifsMCIR[key])
				tAD, probAD = stats.ttest_ind(motifsAD[key], motifsADR[key])
				tNL, probNL = stats.ttest_ind(motifsNL[key], motifsNLR[key])
				motifStats.append((key,probNL,probMCI,probAD))
			
			motifStats.sort(key=lambda x: min(x))
				
			for key, probNL, probMCI, probAD in motifStats:
				normRMean = motifsNLR[key].mean()
				mciRMean = motifsMCIR[key].mean()
				adRMean = motifsADR[key].mean()
				normRVar = motifsNLR[key].std()
				mciRVar = motifsMCIR[key].std()
				adRVar = motifsADR[key].std()
				if probMCI<0.01 or probAD<0.01 or probNL<0.01:
					star = "**"
				elif probMCI<0.1 or probAD<0.1 or probNL<0.01:
					star = "*"
				else:
					star = ""
				line = star+"{0:>"+str(10-len(star))+"}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}{6:>15.3}{7:>15.3}{8:>15.3}{9:>15.3}\n"
				f.write(line.format(str(int(key)), probNL, probMCI, probAD,normRMean,mciRMean,adRMean,normRVar,mciRVar,adRVar))
			f.write("\n\n")	

def motifStats3(data, motifSize=3, degree=10, usetotal=False):
	"""Outputs text file with stats on the motifs in data"""
			
	filename = "result/t_test_Deg-{0}_Size-{1}.txt".format(degree,motifSize)
	with open(filename,'w') as f:
		f.write("Student's T test comparing both AD/NL/MCI to Random.\n\n")
		for corr in ['corr']:
			title = "P-values for "+corr+" data set compared to random generated graphs\n"
			f.write(title)
			
			data[('MCIR', corr)] = genRandomGraphs(data[('MCI', corr)], degree, 119) 
			data[('ADR', corr)] = genRandomGraphs(data[('AD', corr)], degree, 39) 
			data[('NLR', corr)] = genRandomGraphs(data[('NL', corr)], degree, 108)
			 
			motifsNL=findMotifs(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCI=findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
			motifsAD=findMotifs(data,('AD',corr), motifSize, degree, usetotal)
			motifsNLR=findMotifs(data,('NLR',corr), motifSize, degree, usetotal)
			motifsMCIR=findMotifs(data,('MCIR',corr), motifSize, degree, usetotal)
			motifsADR=findMotifs(data,('ADR',corr), motifSize, degree, usetotal)
		
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys())
							& set(motifsNLR.keys())
							& set(motifsADR.keys())
							& set(motifsMCIR.keys()) )
			allMotifs.sort()
			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}{6:>15}{7:>15}{8:>15}{9:>15}\n".format(
				'MOTIF ID','NL', 'MCI','AD', 'NLR Mean','MCIR Mean','ADR Mean', 'NLR Std','MCIR Std', 'ADR Std'))
			
			motifStats = []
			for key in allMotifs:
				tMCI, probMCI = stats.ttest_ind(motifsMCI[key], motifsMCIR[key])
				tAD, probAD = stats.ttest_ind(motifsAD[key], motifsADR[key])
				tNL, probNL = stats.ttest_ind(motifsNL[key], motifsNLR[key])
				motifStats.append((key,probNL,probMCI,probAD))
			
			motifStats.sort(key=lambda x: min(x))
				
			for key, probNL, probMCI, probAD in motifStats:
				normRMean = motifsNLR[key].mean()
				mciRMean = motifsMCIR[key].mean()
				adRMean = motifsADR[key].mean()
				normRVar = motifsNLR[key].std()
				mciRVar = motifsMCIR[key].std()
				adRVar = motifsADR[key].std()
				if probMCI<0.01 or probAD<0.01 or probNL<0.01:
					star = "**"
				elif probMCI<0.1 or probAD<0.1 or probNL<0.01:
					star = "*"
				else:
					star = ""
				line = star+"{0:>"+str(10-len(star))+"}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}{6:>15.3}{7:>15.3}{8:>15.3}{9:>15.3}\n"
				f.write(line.format(str(int(key)), probNL, probMCI, probAD,normRMean,mciRMean,adRMean,normRVar,mciRVar,adRVar))
			f.write("\n\n")	
			
def motifStats3(data, motifSize=3, degree=10, usetotal=False):
	"""Outputs text file with stats on the motifs in data"""
			
	filename = "result/t_test_Deg-{0}_Size-{1}.txt".format(degree,motifSize)
	with open(filename,'w') as f:
		f.write("Student's T test comparing both AD/NL/MCI to Random.\n\n")
		for corr in ['corr']:
			title = "P-values for "+corr+" data set compared to random generated graphs\n"
			f.write(title)
			
			data[('MCIR', corr)] = genRandomGraphs(data[('MCI', corr)], degree, 119) 
			data[('ADR', corr)] = genRandomGraphs(data[('AD', corr)], degree, 39) 
			data[('NLR', corr)] = genRandomGraphs(data[('NL', corr)], degree, 108)
			 
			motifsNL=findMotifs(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCI=findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
			motifsAD=findMotifs(data,('AD',corr), motifSize, degree, usetotal)
			motifsNLR=findMotifs(data,('NLR',corr), motifSize, degree, usetotal)
			motifsMCIR=findMotifs(data,('MCIR',corr), motifSize, degree, usetotal)
			motifsADR=findMotifs(data,('ADR',corr), motifSize, degree, usetotal)
		
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys())
							& set(motifsNLR.keys())
							& set(motifsADR.keys())
							& set(motifsMCIR.keys()) )
			allMotifs.sort()
			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}{6:>15}{7:>15}{8:>15}{9:>15}\n".format(
				'MOTIF ID','NL', 'MCI','AD', 'NLR Mean','MCIR Mean','ADR Mean', 'NLR Std','MCIR Std', 'ADR Std'))
			
			motifStats = []
			for key in allMotifs:
				tMCI, probMCI = stats.ttest_ind(motifsMCI[key], motifsMCIR[key])
				tAD, probAD = stats.ttest_ind(motifsAD[key], motifsADR[key])
				tNL, probNL = stats.ttest_ind(motifsNL[key], motifsNLR[key])
				motifStats.append((key,probNL,probMCI,probAD))
			
			motifStats.sort(key=lambda x: min(x))
				
			for key, probNL, probMCI, probAD in motifStats:
				normRMean = motifsNLR[key].mean()
				mciRMean = motifsMCIR[key].mean()
				adRMean = motifsADR[key].mean()
				normRVar = motifsNLR[key].std()
				mciRVar = motifsMCIR[key].std()
				adRVar = motifsADR[key].std()
				if probMCI<0.01 or probAD<0.01 or probNL<0.01:
					star = "**"
				elif probMCI<0.1 or probAD<0.1 or probNL<0.01:
					star = "*"
				else:
					star = ""
				line = star+"{0:>"+str(10-len(star))+"}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}{6:>15.3}{7:>15.3}{8:>15.3}{9:>15.3}\n"
				f.write(line.format(str(int(key)), probNL, probMCI, probAD,normRMean,mciRMean,adRMean,normRVar,mciRVar,adRVar))
			f.write("\n\n")	
			
def choose_and_remove( items ):
    # pick an item index
    if items:
        index = random.randrange( len(items) )
        return items.pop(index)
    # nothing left!
    return None
			
def createfakeGroups(data):
	newdata ={}
	for corr in  ('corr','lcorr','lacorr'):
		newdata[('NL', corr)] = []
		newdata[('AD', corr)] = []
		newdata[('MCI', corr)] = []
		
		adlen = len(data[('AD', corr)]) 
		nllen = len(data[('NL', corr)])
		mcilen = len(data[('MCI', corr)])
		total = adlen + nllen + mcilen
		
		myList = [nllen, adlen, mcilen]
		for i in range(3):
			n = int(float(myList[i])/float(total) * nllen)
			for j in range(n):
				if i == 0:
					g = choose_and_remove(data[('NL', corr)])
					newdata[('NL', corr)].append(g)
				if i == 1:
					g = choose_and_remove(data[('AD', corr)])
					newdata[('NL', corr)].append(g)
				if i == 2:
					g = choose_and_remove(data[('MCI', corr)])
					newdata[('NL', corr)].append(g)
		
		for i in range(3):
			n = int(float(myList[i])/float(total) * adlen)
			for j in range(n):
				if i == 0:
					g = choose_and_remove(data[('NL', corr)])
					newdata[('AD', corr)].append(g)
				if i == 1:
					g = choose_and_remove(data[('AD', corr)])
					newdata[('AD', corr)].append(g)
				if i == 2:
					g = choose_and_remove(data[('MCI', corr)])
					newdata[('AD', corr)].append(g)
					
		for i in range(3):
			n = int(float(myList[i])/float(total) * mcilen)
			for j in range(n):
				if i == 0:
					g = choose_and_remove(data[('NL', corr)])
					newdata[('MCI', corr)].append(g)
				if i == 1:
					g = choose_and_remove(data[('AD', corr)])
					newdata[('MCI', corr)].append(g)
				if i == 2:
					g = choose_and_remove(data[('MCI', corr)])
					newdata[('MCI', corr)].append(g)

	return newdata
		
def motifStats(data, motifSize=3, degree=10, usetotal=False):
	"""Outputs pdf file with stats on the motifs in data"""
			
	filename = "result/t_test_Deg-{0}_Size-{1}.txt".format(degree,motifSize)
	with open(filename,'w') as f:
		f.write("Student's T test comparing both MCI and AD to NL.\n\n")
		for corr in ('corr','lcorr','lacorr'):
			title = "P-values for "+corr+" data set compared to normal patients\n"
			f.write(title)
	
			
			motifsNL=findMotifs(data,('NL',corr), motifSize, degree, usetotal)
			motifsMCI=findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
			motifsAD=findMotifs(data,('AD',corr), motifSize, degree, usetotal)
			
			#mats = []
			#for i in xrange(108):	
			#	x = np.random.rand(88,88)
			#	x -= np.diag(np.diag(x))
			#	mats.append(x)
			#rand = {}
			#rand['derp'] = mats
			#motifsNL=findMotifs(rand,'derp', motifSize, degree, usetotal)
		
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys()) )

			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}{6:>15}{7:>15}{8:>15}\n".format(
				'MOTIF ID','MCI','AD','Norm Mean','MCI Mean','AD Mean','NORM Std','MCI Std', 'AD Std'))
			
			motifStats = []
			for key in allMotifs:
				tMCI, probMCI = stats.ttest_ind(motifsMCI[key], motifsNL[key])
				tAD, probAD = stats.ttest_ind(motifsAD[key], motifsNL[key])
				motifStats.append((key,probMCI,probAD))
			
			motifStats.sort(key=lambda x: min(x))
				
			for key, probMCI, probAD in motifStats:
				normMean = motifsNL[key].mean()
				mciMean = motifsMCI[key].mean()
				adMean = motifsAD[key].mean()
				normVar = motifsNL[key].std()
				mciVar = motifsMCI[key].std()
				adVar = motifsAD[key].std()
				if probMCI<0.01 or probAD<0.01:
					star = "**"
				elif probMCI<0.1 or probAD<0.1:
					star = "*"
				else:
					star = ""
				line = star+"{0:>"+str(10-len(star))+"}{1:>15.3}{2:>15.3}{3:>15.3}{4:>15.3}{5:>15.3}{6:>15.3}{7:>15.3}{8:>15.3}\n"
				f.write(line.format(str(int(key)), probMCI, probAD,normMean,mciMean,adMean,normVar,mciVar,adVar))
			f.write("\n\n")

def genRandMats(num):
	mats = []
	for i in xrange(num):
		x = np.random.rand(88,88)
		x -= np.diag(np.diag(x))
		mats.append(x)
	return mats
					
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
		keys = [int(key) for key,value in motifs[:numofmotifs]]
		
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
		
def PDFstats(data,filename):
	
	data['rand'] = genRandMats(100)
	motifsRAND = findMotifs(data,'rand')
	
	filename = "result/" + filename + ".tex"
	
	with open(filename,'wb') as f:
		f.write(
		"""
		\\documentclass{article}
		\\usepackage{amsmath,fullpage,graphicx,fancyhdr,xcolor,colortbl}
		\\definecolor{yellow}{rgb}{1,1,0}
		\\definecolor{orange}{rgb}{1,0.647,0}
		\\definecolor{red}{rgb}{1,0,0}
		\\title{Motif Data}
		\\author{Graham Tremper}
		\\date{}
		\\fancyhead{}
		\\begin{document}
		""")
		for corr in ('corr','lcorr','lacorr'):
			motifsNL = findMotifs(data, ('NL',corr))
			motifsMCI = findMotifs(data, ('MCI',corr))
			motifsAD = findMotifs(data, ('AD',corr))
			
			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys()) )
							
			motifStats = []
			for key in allMotifs:
				c1 = stats.ttest_ind(motifsMCI[key], motifsNL[key])
				c2 = stats.ttest_ind(motifsAD[key], motifsNL[key])
				c3 = stats.ttest_ind(motifsMCI[key], motifsAD[key])
				c4 = stats.ttest_ind(motifsNL[key], motifsRAND[key])
				c5 = stats.ttest_ind(motifsMCI[key], motifsRAND[key])
				c6 = stats.ttest_ind(motifsAD[key], motifsRAND[key])
				motifStats.append((key,c1,c2,c3,c4,c5,c6))
			
			motifStats.sort(key=lambda x: motifsNL[x[0]].mean(),reverse=True)
						
			f.write(
			"""
			\\begin{table}[t]
			""" + "\\caption{Motif T-test results from "+corr+" data}" +
			"""
			\\centering
			\\begin{tabular}{|c|c|c|c|c|c|c|}
			\\hline
			\\rowcolor[gray]{0.85} 
			Key	& MCI to Norm & AD to Norm & MCI to AD & Norm to Rand & MCI to Rand & AD to Rand \\\\ \\hline
			""")
			for stat in motifStats:
				f.write( str(stat[0]) + " \\cellcolor[gray]{0.95}")
				for sign,col in stat[1:]:
					cell = " & {0:.3}".format(col)
					if sign > 0:
						cell += '(+)'
					else:
						cell += '(-)'
					
					if col <= 0.01:
						cell += " \\cellcolor{red} "
					elif col <= 0.05:
						cell += " \\cellcolor{orange}"
					elif col <= 0.1:
						cell += " \\cellcolor{yellow}"
					f.write(cell)
				f.write("\\\\ \\hline\n")
				
			f.write(
			"""
			\\end{tabular}
			\\end{table}
			""")
		
		f.write("\\end{document}\n")
	
	os.system("pdflatex -output-directory result " + filename)
	os.system("rm result/*.tex result/*.log result/*.aux")

if __name__ == '__main__':
	with open("aznorbert_corrsd.pkl","rb") as f:
		data = pickle.load(f)
	
	PDFstats(data,"Motif_Statistics")
	#motifStats(data)
	
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
	

	
	
	
