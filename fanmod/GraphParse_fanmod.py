import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import cPickle as pickle
import scipy.stats as stats
import re

USECACHE = False

#Also, I thought of a good test for your programs.
#Generate 100 random correlation matrices via
#x=np.random.rand(n,n)
#
#then create one group from
#x -= np.diag(np.diag(x))
#
#and the other from
#x = np.triu(x,1)
#x+= x.T
#
#and run all you analyses on these and compare.

def generateRandonGraphsStats():
	graphdict = {}
	
	for i in range(100):
		x = np.random.rand(100, 100)
		y = x - np.diag(np.diag(x))
		if i == 0:
			graphdict['diag'] = [y]
		else:
			graphdict['diag'].append(y)
		
		z = np.triu(x,1)
		z += z.T

		if i == 0:
			#print z
			graphdict['triu'] = [z]
		else:
			graphdict['triu'].append(z)
	plotMotifGraphs2(graphdict, 3, 10, 10)
	plotMotifGraphs2(graphdict, 3, 10, 10, useFanmod=False)
	

def plotMotifGraphs2(data,motifSize,degree,numofmotifs,usetotal=False, useFanmod=True):
	for corr in ('diag','triu'):
		nl=findMotifs(data, corr, motifSize, degree, usetotal, useFanmod)
		
		motifs = nl.items()
		motifs = sorted(motifs,key=lambda x:-x[1].mean())
		keys = [int(x[0]) for x in motifs[:numofmotifs]]
		
		meansNL = []
		stdNL = []
		
		while len(keys) < numofmotifs:
			keys.append(0)
		
		
		for key in keys:			
			if key in nl:
				meansNL.append(nl[key].mean())
				stdNL.append(nl[key].std())
			else:
				meansNL.append(0.0)
				stdNL.append(0.0)
		
		ind = np.arange(numofmotifs)
		width = 0.3 
		print meansNL
		NLplt = plt.bar(ind + 0.3, meansNL, width, color='b', yerr=stdNL, ecolor='y')

		plt.ylabel('Average number of motifs')
		plt.xlabel('Motif ID')
		graphname = 'Fanmod'
		if useFanmod == False:
			graphname = 'Kavosh'
		plt.title(graphname+': Motif size '+str(motifSize) +' distribution for '+corr+" with average degree "+str(degree))
		plt.xticks(ind+width+width/2., keys)
		plt.ylim(ymin=0.0)
#		plt.legend( (NLplt[0], MCIplt[0], ADplt[0]), ('NL', 'MCI', 'AD') )
		plt.grid(True)

		resultdir = 'result2'
		if useFanmod == False:
			resultdir = 'result'
		if usetotal:
			plt.savefig(resultdir + "/TotalMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		else:
			plt.savefig(resultdir + "/PercentMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		plt.clf()

def parseOutput(path):
	pattern = re.compile('%')
	returnlist = []
	f = open(path, 'r')
	total = None
	for line in f:
		myList = line.strip().split(",")
		if len(myList) == 1 and myList[0].find('original network') != -1:
			total = float(line.strip().split(" ")[0])
			#print total
		elif len(myList) == 7 and myList[0] != 'ID':
			estimate = float(re.sub(pattern, '', myList[2]))/100 * total
			#print estimate
			returnlist.append([int(myList[0]), estimate, float(re.sub(pattern, '', myList[2]))/100])
			
	return returnlist

#Out put a text file representation of graph
def outputGraph(graph,name="OUTPUT.txt"):
	G = nx.convert_node_labels_to_integers(graph,1)
	with open("result2/"+name,'wb') as f:
		nx.write_edgelist(G,f,data=False)

#Find motifs of motifSize in G
def findSingleMotif(G,motifSize):
	outputGraph(G)
	os.system("./fanmod_command_line_linux " +str(motifSize) + " 100000 1 result2/OUTPUT.txt 1 0 0 2 0 0 0 1 3 3 result2/MotifCount.txt 0 1")
	data = parseOutput("result2/MotifCount.txt")
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
	#plt.savefig("result/id-"+str(id)+"size-"+str(motifSize))
	plt.show()

# Main finding motifs routine	
def findMotifs(data,key,motifSize,degree,usetotal=False, useFanmod=True):
	#Check cache
	filename = str(key)+'s'+str(int(motifSize))+'d'+str(int(degree))+str(usetotal)+".pkl"
	cachedir = 'cache2/'
	if useFanmod == False:
		cachedir = 'cache/'
		
	if os.path.exists(cachedir+filename) and USECACHE:
		print "in cache"
		with open(cachedir+filename,"rb") as f:
			return pickle.load(f)
			
	graphs = data[key]
	motifs = {}
	numstring ="/"+str(len(graphs))
	counter = 1
	for G in graphs:
		#calculate threshold
		sortedWeights = np.sort(G,axis=None)
#		print G
#		print len(sortedWeights)
#		print str(-len(G)*degree-1)
		threshold = sortedWeights[-len(G)*degree]
		#Print progress
		sys.stdout.write("\r")
		sys.stdout.write("Motif Finding Progress: "+str(counter)+numstring)
		sys.stdout.write(" Threshold: "+str(threshold))
		sys.stdout.flush()
		counter+=1
		#Output graph to txt file
		graph = nx.DiGraph(G>threshold)
		graph = nx.convert_node_labels_to_integers(graph,1)
		
		resultdir = 'result2'
		if useFanmod == False:
			resultdir = 'result'
		with open(resultdir+'/OUTPUT.txt','wb') as f:
			nx.write_edgelist(graph,f,data=False)
		#Jenky way to use c++ motif finder in python
		if useFanmod:
			os.system("./fanmod_command_line_linux " +str(motifSize) + " 100000 1 result2/OUTPUT.txt 1 0 0 2 0 0 0 0 3 3 result2/MotifCount.txt 0 1")
			data = parseOutput("result2/MotifCount.txt")
		else:
			os.system("./Kavosh " + str(motifSize))
			data = np.loadtxt("result/MotifCount.txt",ndmin=2)
#			print data
			
		for iD,total,percent in data:
			if iD in motifs:
				motifs[int(iD)].append(total if usetotal else percent)
			else:
				motifs[int(iD)] = [total if usetotal else percent]
		
		
	print '\nMotifs Done!'
	
	#add zeros to graphs that didn't contain motifs
	for key,value in motifs.iteritems():
		numZero = len(graphs)-len(value)
		value.extend([0 for derp in xrange(numZero)])
		motifs[int(key)] = np.array(value)
	
	#add motifs to cache
	if USECACHE:
		with open(cachedir+filename,'wb') as f:
			pickle.dump(motifs,f)
			
	return motifs

def motifOrder(data,key,orderSize,motifSize,degree):			
	graphs = data[key]
	pattern = {}
	for G in graphs:
		#calculate threshold
		sortedWeights = np.sort(G,axis=None)
		threshold = sortedWeights[-len(G)*degree-1]
		#Output graph to txt file
		graph = nx.DiGraph(G>threshold)
		graph = nx.convert_node_labels_to_integers(graph,1)
		with open('result2/OUTPUT.txt','wb') as f:
			nx.write_edgelist(graph,f,data=False)
		#Jenky way to use c++ motif finder in python
		os.system("./fanmod_command_line_linux " +str(motifSize) + " 100000 1 result2/OUTPUT.txt 1 0 0 2 0 0 0 1 3 3 result2/MotifCount.txt 0 1")
		data = parseOutput("result2/MotifCount.txt")
		
		order = []
		for iD,total,percent in data:
			order.append((iD,total))
		keys = sorted(order,key=lambda x:-x[1])
		keys = [int(k[0]) for k in keys]
		pat = tuple(keys[:orderSize])
		pattern[pat] = pattern.setdefault(pat,0) + 1/float(len(graphs))
		
		
	total = sorted(pattern.items(), key = lambda x: -x[1])
	
	for key,value in total:
		print str(key)+": " + str(value)

		

def motifStats(data,motifSize,degree, usetotal=False):
	
	for corr in ('corr','lcorr','lacorr'):
		motifsNL = findMotifs(data,('NL',corr), motifSize, degree, usetotal)
		motifsMCI = findMotifs(data,('MCI',corr), motifSize, degree, usetotal)
		motifsAD = findMotifs(data,('AD',corr), motifSize, degree, usetotal)
		
		allMotifs = list(set(motifsNL.keys()) | set(motifsAD.keys()) | set(motifsMCI.keys()))
		
		datatype = "Total" if usetotal else "Percent"
		filename = "result2/{}_ks-stats_size-{}_deg-{}.txt".format(corr+datatype,motifSize,degree)
		with open(filename,'w') as f:
			f.write("{0:>10}{1:>15}{2:>15}{3:>15}{4:>15}{5:>15}\n".format('ID','MCI','AD','NORM NL','NORM MCI','NORM AD'))
			for key in allMotifs:
				NLdata = motifsNL.get(key,np.zeros(88))
				MCIdata = motifsMCI.get(key,np.zeros(88))
				ADdata = motifsAD.get(key,np.zeros(88))
				KSstatistic, MCIpvalue = stats.ks_2samp(MCIdata,NLdata)
				KSstatistic, ADpvalue = stats.ks_2samp(ADdata,NLdata)
				k2,NLnorm = stats.normaltest(NLdata)
				k2,MCInorm = stats.normaltest(MCIdata)
				k2,ADnorm = stats.normaltest(ADdata)
				if MCIpvalue<0.01 or ADpvalue<0.01:
					line = "*{0:>9}{1:15.3}{2:15.3}{3:15.3}{4:15.3}{5:15.3}\n"
				else:
					line = "{0:>10}{1:15.3}{2:15.3}{3:15.3}{4:15.3}{5:15.3}\n"
				f.write(line.format(str(int(key)),MCIpvalue,ADpvalue,NLnorm,MCInorm,ADnorm))

					
def plotMotifGraphs(data,motifSize,degree,numofmotifs,usetotal=False):
	for corr in ('corr','lcorr','lacorr'):
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
			if key in nl:
				meansNL.append(nl[key].mean())
				stdNL.append(nl[key].std())
			else:
				meansNL.append(0.0)
				stdNL.append(0.0)
			if key in mci:
				meansMCI.append(mci[key].mean())
				stdMCI.append(mci[key].std())
			else:
				meansMCI.append(0.0)
				stdMCI.append(0.0)
			if key in ad:
				meansAD.append(ad[key].mean())
				stdAD.append(ad[key].std())
			else:
				meansAD.append(0.0)
				stdAD.append(0.0)
		
		ind = np.arange(numofmotifs)
		width = 0.2 

		NLplt = plt.bar(ind, meansNL, width, color='b', yerr=stdNL, ecolor='y')
		MCIplt = plt.bar(ind+width, meansMCI, width, color='y', yerr=stdMCI, ecolor='b')
		ADplt = plt.bar(ind+width+width, meansAD, width, color='g', yerr=stdAD, ecolor='r')

		plt.ylabel('Average number of motifs')
		plt.xlabel('Motif ID')
		plt.title('Fanmod: Motif size '+str(motifSize) +' distribution for '+corr+" with average degree "+str(degree))
		plt.xticks(ind+width+width/2., keys)
		plt.ylim(ymin=0.0)
		plt.legend( (NLplt[0], MCIplt[0], ADplt[0]), ('NL', 'MCI', 'AD') )
		plt.grid(True)
		if usetotal:
			plt.savefig("result2/TotalMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		else:
			plt.savefig("result2/PercentMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		plt.clf()

def addRandom(data):
	for corr in ['corr', 'lcorr', 'lacorr']:
		for i in range(100):
			x = np.random.rand(88, 88)
			y = x - np.diag(np.diag(x))
			if i == 0:
				data[('RAN', corr)] = [y]
			else:
				data[('RAN', corr)].append(y)
	return data
	
def plotMotifGraphsWithRandom(data,motifSize,degree,numofmotifs,usetotal=False):
	data = addRandom(data);
	for corr in ('corr','lcorr','lacorr'):
		nl=findMotifs(data,('NL',corr), motifSize, degree,usetotal)
		mci=findMotifs(data,('MCI',corr), motifSize, degree,usetotal)
		ad=findMotifs(data,('AD',corr), motifSize, degree,usetotal)
		ran=findMotifs(data, ('RAN', corr), motifSize, degree, usetotal)
		
		motifs = nl.items()
		motifs = sorted(motifs,key=lambda x:-x[1].mean())
		keys = [int(x[0]) for x in motifs[:numofmotifs]]
		
		meansNL = []
		meansMCI = []
		meansAD = []
		meansRAN = []
		stdNL = []
		stdMCI = []
		stdAD = []
		stdRAN = []
		for key in keys:			
			if key in nl:
				meansNL.append(nl[key].mean())
				stdNL.append(nl[key].std())
			else:
				meansNL.append(0.0)
				stdNL.append(0.0)
			if key in mci:
				meansMCI.append(mci[key].mean())
				stdMCI.append(mci[key].std())
			else:
				meansMCI.append(0.0)
				stdMCI.append(0.0)
			if key in ad:
				meansAD.append(ad[key].mean())
				stdAD.append(ad[key].std())
			else:
				meansAD.append(0.0)
				stdAD.append(0.0)
			if key in ran:
				meansRAN.append(ran[key].mean())
				stdRAN.append(ran[key].std())
			else:
				meansRAN.append(0.0)
				stdRAN.append(0.0)
		
		ind = np.arange(numofmotifs)
		width = 0.15 

		NLplt = plt.bar(ind, meansNL, width, color='b', yerr=stdNL, ecolor='y')
		MCIplt = plt.bar(ind+width, meansMCI, width, color='y', yerr=stdMCI, ecolor='b')
		ADplt = plt.bar(ind+width+width, meansAD, width, color='g', yerr=stdAD, ecolor='r')
		RANplt = plt.bar(ind+width*3, meansRAN, width, color='m', yerr=stdRAN, ecolor='k')

		plt.ylabel('Average number of motifs')
		plt.xlabel('Motif ID')
		plt.title('Fanmod: Motif size '+str(motifSize) +' distribution for '+corr+" with average degree "+str(degree))
		plt.xticks(ind+width+width/2., keys)
		plt.ylim(ymin=0.0)
		plt.legend( (NLplt[0], MCIplt[0], ADplt[0], RANplt[0]), ('NL', 'MCI', 'AD', 'RANDOM') )
		plt.grid(True)
		if usetotal:
			plt.savefig("result2/TotalMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		else:
			plt.savefig("result2/PercentMotifDist-"+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		plt.clf()

if __name__ == '__main__':
	with open("aznorbert_corrsd.pkl","rb") as f:
		data = pickle.load(f)
		
	plotMotifGraphsWithRandom(data,3,10,10)
#	plotMotifGraphs(data,3,12,10)
#	plotMotifGraphs(data,3,15,10)
#	plotMotifGraphs(data,4,10,10)

#	print 'Normal'
#	motifOrder(data,('NL','corr'),4 , 3, 10)
#	print 'MCI'
#	motifOrder(data,('MCI','corr'),4, 3, 10)
#	print 'AD'
#	motifOrder(data,('AD','corr'),4, 3, 10)

#	keys = data.keys()
#	G = data[keys[0]][0]
#	degree = 10
#	sortedWeights = np.sort(G,axis=None)
#	threshold = sortedWeights[-len(G)*degree-1]
#
#	graph = nx.DiGraph(G>threshold)
#	findSingleMotif(graph, 3)

#	statTest(data,3,15)
#	statTest(data,4,10)
#	statTest(data,4,15)

	
	#generateRandonGraphsStats()
	
	
	
