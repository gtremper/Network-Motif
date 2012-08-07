#!/usr/bin/env python
# encoding: utf-8
"""
FinalMotif.py

Created by BNIAC on 2012-08-02.
"""

import sys
import os
import random
import math
import heapq
import json
import cPickle as pickle
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import graph_helper as gh
from scipy.sparse import dok_matrix
from itertools import izip


class MotifData:
	"Class containing motif data for a set of graphs"

	def __init__(self, data):
		self.subgraphs, self.data = data
		allkeys = set()
		for dic in self.data:
			allkeys.update(set(dic.keys()))
		self.keys = allkeys
			
	def __getitem__(self, motif):
		"Get array of number of motif for each patient"
		motif = unicode(motif)
		row = [d[motif] if motif in d else 0. for d in self.data]
		return np.array(row)
	
	def __contains__(self, item):
		return self.keys.__contains__(unicode(item))
	
	def __iter__(self):
		return iter(self.data)
	
	def __len__(self):
		return len(self.data)
	
	def iterSortedValues(self):
		"Iterate sorted motif values for each patient"
		return (sorted(g.values()) for g in iter(self.data))
	
	def topMotifs(self, num):
		"Return top num motifs"
		return heapq.nlargest(num, self.keys, key = lambda x: self[x].mean())
	
	def iterTotals(self):
		"Iterate through total number of motifs for each patient"
		for d,sub in izip(self.data,self.subgraphs):
			t = {}
			for motif,value in d.iteritems():
				t[int(i)] = int(value*sub+0.1)
			yield t

	def getPatient(self, pat):
		return self.data[pat]

	def getSubgraphs(self, pat):
		return self.subgraphs[pat]


def findMotifs(data,key,motifSize=3,degree=10,randGraphs=None, useCache=True):
	"""Main finding motifs routine"""

	#generate random matricies
	if key == "rand":
		graphs = genRandMats(100)
	else:
		graphs = data[key]

	#Check cache
	filename = "" if randGraphs is None else "RAND"
	filename += str(key)+'s'+str(int(motifSize))+'d'+str(int(degree))+".json"
	if os.path.exists('MotifCache/'+filename) and useCache:
		print "in cache"
		cachedata = json.load( open('MotifCache/'+filename,"rb"))
		return MotifData(cachedata)

	motifs = []
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
		if randGraphs is not None:
			graph = randGraphs[key][index]
		with open('result/OUTPUT.txt','wb') as f:
			f.write(str(len(graph)) + '\n')
			nx.write_edgelist(graph,f,data=False)

		#Jenky way to use c++ motif finder in python
		os.system("./Kavosh "+str(motifSize))
		with open("result/MotifCount.txt","rb") as f:
			subgraphs = float(f.next())
			data = np.loadtxt(f, ndmin=2)

		#Append data for this graph
		personMotifs = {}
		for iD,total in data:
			personMotifs[unicode(int(iD))] = total/subgraphs
		motifs.append((int(subgraphs),personMotifs))

	print '\nMotifs Done! Graphs Rejected: '+str(rejected)

	#add motifs to cache
	if useCache:
		if not os.path.isdir('MotifCache'):
			os.makedirs('MotifCache')
			
		json.dump(motifs, open('MotifCache/'+filename,'wb'), separators=(',',':'))

	return MotifData(motifs)


def convertIDToGraph(mid, motifSize, save=False):
	"""Draw graph with id and motifSize"""
	binary = bin(mid);
	adj = np.zeros(motifSize*motifSize)
	l = 0
	for x in xrange(1,motifSize*motifSize+1):
		if binary[-x+l] == 'b':
			break
		if (x-1) % (motifSize+1) == 0:
			l += 1
		else:
			adj[-x] = int(binary[-x+l])
	adj.shape = (motifSize,motifSize)
	graph = nx.to_networkx_graph(adj,create_using=nx.DiGraph())
	nx.draw_circular(graph)
	if save:
		plt.savefig("result/id-"+str(id)+"size-"+str(motifSize))
	else:
		plt.show()
	plt.clf()


def plotMotifGraphs(data,motifSize=3,degree=10,numofmotifs=10,usetotal=False):
	"""Draws graph compairing average motif count between samples in the data"""
	for corr in ('corr','lcorr','lacorr'):

		nl=findMotifs(data,('NL',corr), motifSize, degree)
		mci=findMotifs(data,('MCI',corr), motifSize, degree)
		ad=findMotifs(data,('AD',corr), motifSize, degree)
		convert = findMotifs(data,('CONVERT',corr), motifSize, degree)
		
		keys = nl.topMotifs(numofmotifs)

		meansNL = []
		meansMCI = []
		meansAD = []
		meansCONVERT = []
		stdNL = []
		stdMCI = []
		stdAD = []
		stdCONVERT = []
		for key in keys:
			meansNL.append(nl[key].mean() if key in nl else 0.)
			stdNL.append(nl[key].std() if key in nl else 0.)
			meansMCI.append(mci[key].mean() if key in mci else 0.)
			stdMCI.append(mci[key].std() if key in mci else 0.)
			meansAD.append(ad[key].mean() if key in mci else 0.)
			stdAD.append(ad[key].std() if key in mci else 0.)
			meansCONVERT.append(convert[key].mean() if key in convert else 0.)
			stdCONVERT.append(convert[key].std() if key in convert else 0.)

		ind = np.arange(numofmotifs)
		width = 0.2 

		NLplt = plt.bar(ind, meansNL, width, color='b', yerr=stdNL, ecolor='y')
		MCIplt = plt.bar(ind+width, meansMCI, width, color='y', yerr=stdMCI, ecolor='b')
		ADplt = plt.bar(ind+width+width, meansAD, width, color='g', yerr=stdAD, ecolor='r')
		CONVERTplt = plt.bar(ind+width+width+width, meansCONVERT, width, color='r', yerr=stdCONVERT, ecolor='y')

		plt.ylabel('Average number of motifs')
		plt.xlabel('Motif ID')
		plt.title('Motif size '+str(motifSize) +' distribution for '+corr+" with average degree "+str(degree))
		plt.xticks(ind+width+width/2., keys)
		plt.ylim(ymin=0.0)
		plt.legend( (NLplt[0], MCIplt[0], ADplt[0], CONVERTplt[0]), ('NL', 'MCI', 'AD', 'CONVERT') )
		plt.grid(True)
		header = 'result/MotifDistribution-'
		plt.savefig(header+corr+"_D-"+str(degree)+"_S-"+str(motifSize))
		plt.clf()


def PDFstats(data, filename, edgeSwap=False, motifSize=3, degree=10):
	"""Output a latex pdf of motif stats"""
	filename = "result/" + filename + ".tex"

	if not edgeSwap:
		motifsNLRAND = motifsMCIRAND = motifsADRAND = motifsCONVERTRAND = findMotifs(data,"rand",motifSize=motifSize,degree=degree)

	with open(filename,'wb') as f:
		f.write(
		"\\documentclass{article}\n"
		"\\usepackage{amsmath,fullpage,graphicx,fancyhdr,xcolor,colortbl,chngpage}\n"
		"\\usepackage[landscape]{geometry}"
		"\\definecolor{yellow}{RGB}{255,255,70}\n"
		"\\definecolor{orange}{RGB}{255,165,70}\n"
		"\\definecolor{red}{RGB}{255,70,70}\n"
		"\\title{Motif Data}\n"
		"\\author{Graham Tremper}\n"
		"\\date{}\n"
		"\\fancyhead{}\n"
		"\\begin{document}\n"
		)

		if edgeSwap:
			with open("SwapData"+str(degree)+".pkl","rb") as pic:
				randGraphs = pickle.load(pic)

		statistics = {}		
		for corr in ('corr','lcorr','lacorr'):
			print "Starting " + corr +"..."
			motifsNL = findMotifs(data, ('NL',corr), motifSize = motifSize, degree=degree)
			motifsMCI = findMotifs(data, ('MCI',corr), motifSize = motifSize, degree=degree)
			motifsAD = findMotifs(data, ('AD',corr), motifSize = motifSize, degree=degree)
			motifsCONVERT = findMotifs(data, ('CONVERT',corr), motifSize = motifSize, degree=degree)
			if edgeSwap:
				motifsNLRAND = findMotifs(data, ('NL',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				motifsMCIRAND = findMotifs(data, ('MCI',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				motifsADRAND = findMotifs(data, ('AD',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				motifsCONVERTRAND = findMotifs(data, ('CONVERT',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)

			allMotifs = list( motifsNL.keys
							& motifsAD.keys
							& motifsMCI.keys
							& motifsCONVERT.keys
							& motifsNLRAND.keys
							& motifsMCIRAND.keys
							& motifsADRAND.keys
							& motifsCONVERTRAND.keys )

			allMotifs = motifsNL.topMotifs(30)

			motifStats = []
			for key in allMotifs:
				norm = motifsNL[key]
				mci = motifsMCI[key]
				ad = motifsAD[key]
				conv = motifsCONVERT[key]
				c1 = stats.ttest_ind(norm, mci)
				c2 = stats.ttest_ind(norm, ad)
				c3 = stats.ttest_ind(norm, conv)
				c4 = stats.ttest_ind(mci, ad)
				c5 = stats.ttest_ind(mci, conv)
				c6 = stats.ttest_ind(ad, conv)
				c7 = stats.ttest_ind(norm, motifsNLRAND[key])
				c8 = stats.ttest_ind(mci, motifsMCIRAND[key])
				c9 = stats.ttest_ind(ad, motifsADRAND[key])
				c10 = stats.ttest_ind(conv, motifsCONVERTRAND[key])
				motifStats.append((key,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10))

			f.write(
			"\\begin{table}[t]\n"
			"\\begin{adjustwidth}{-1.5in}{-1.5in} "
			"\\caption{Motif T-test results from "+corr+" data with using edge swap}\n"
			"\\centering\n"
			"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\n"
			"\\hline\n"
			"\\rowcolor[gray]{0.85}\n"
			"Key & NL to MCI & NL to AD & NL to Conv & MCI to AD & MCI to Conv & AD to Conv & NL to Rand & MCI to Rand & AD to Rand & Conv to Rand \\\\ \\hline\n"
			)
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
			"\\end{tabular}\n"
			"\\end{adjustwidth}"
			"\\end{table}\n"
			)

		f.write("\\end{document}\n")

	os.system("pdflatex -output-directory result " + filename)
	os.system("rm result/*.log result/*.aux")


def createfakeGroups(data, motifSize):
	"Create new random groups by proportionally shuffling all the groups together"
	newdata ={}
	for corr in	 ('corr','lcorr','lacorr'):
		newdata[('NL', corr)] = []
		newdata[('AD', corr)] = []
		newdata[('MCI', corr)] = []
		newdata[('CONVERT', corr)] = []

		nlData = list(findMotifs(data,('NL',corr), motifSize=motifSize).data)
		adData = list(findMotifs(data,('AD',corr), motifSize=motifSize).data)
		mciData = list(findMotifs(data,('MCI',corr), motifSize=motifSize).data)
		convertData = list(findMotifs(data,('CONVERT',corr), motifSize=motifSize).data)

		adlen = len(adData) 
		nllen = len(nlData)
		mcilen = len(mciData)
		convertlen = len(convertData)
		total = adlen + nllen + mcilen + convertlen

		myList = [nllen, adlen, mcilen, convertlen]
		for i in xrange(4):
			n = int(float(myList[i])/float(total) * nllen)
			for j in xrange(n):
				if i == 0:
					g = nlData.pop(random.randrange(len(nlData)))
					newdata[('NL', corr)].append(g)
				if i == 1:
					g = adData.pop(random.randrange(len(adData)))
					newdata[('NL', corr)].append(g)
				if i == 2:
					g = mciData.pop(random.randrange(len(mciData)))
					newdata[('NL', corr)].append(g)
				if i == 3:
					g = convertData.pop(random.randrange(len(convertData)))
					newdata[('NL', corr)].append(g)

		for i in xrange(4):
			n = int(float(myList[i])/float(total) * adlen)
			for j in xrange(n):
				if i == 0:
					g = nlData.pop(random.randrange(len(nlData)))
					newdata[('AD', corr)].append(g)
				if i == 1:
					g = adData.pop(random.randrange(len(adData)))
					newdata[('AD', corr)].append(g)
				if i == 2:
					g = mciData.pop(random.randrange(len(mciData)))
					newdata[('AD', corr)].append(g)
				if i == 3:
					g = convertData.pop(random.randrange(len(convertData)))
					newdata[('AD', corr)].append(g)

		for i in xrange(4):
			n = int(float(myList[i])/float(total) * mcilen)
			for j in xrange(n):
				if i == 0:
					g = nlData.pop(random.randrange(len(nlData)))
					newdata[('MCI', corr)].append(g)
				if i == 1:
					g = adData.pop(random.randrange(len(adData)))
					newdata[('MCI', corr)].append(g)
				if i == 2:
					g = mciData.pop(random.randrange(len(mciData)))
					newdata[('MCI', corr)].append(g)
				if i == 3:
					g = convertData.pop(random.randrange(len(convertData)))
					newdata[('MCI', corr)].append(g)

		for i in xrange(4):
			n = int(float(myList[i])/float(total) * convertlen)
			for j in xrange(n):
				if i == 0:
					g = nlData.pop(random.randrange(len(nlData)))
					newdata[('CONVERT', corr)].append(g)
				if i == 1:
					g = adData.pop(random.randrange(len(adData)))
					newdata[('CONVERT', corr)].append(g)
				if i == 2:
					g = mciData.pop(random.randrange(len(mciData)))
					newdata[('CONVERT', corr)].append(g)
				if i == 3:
					g = convertData.pop(random.randrange(len(convertData)))
					newdata[('CONVERT', corr)].append(g)

		leftovers = nlData + adData + mciData + convertData
		while len(newdata['NL', corr]) < nllen:
			g = leftovers.pop(random.randrange(len(leftovers)))
			newdata[('NL', corr)].append(g)

		while len(newdata['AD', corr]) < adlen:
			g = leftovers.pop(random.randrange(len(leftovers)))
			newdata[('AD', corr)].append(g)

		while len(newdata['MCI', corr]) < mcilen:
			g = leftovers.pop(random.randrange(len(leftovers)))
			newdata[('MCI', corr)].append(g)

		while len(newdata['CONVERT', corr]) < convertlen:
			g = leftovers.pop(random.randrange(len(leftovers)))
			newdata[('CONVERT', corr)].append(g)

	return newdata

def PDFstatsShuf(data, filename, motifSize=3, degree=10):
	"""Output a latex pdf of motif stats with the patient groups shuffled"""
	filename = "result/" + filename + ".tex"
	shufData = createfakeGroups(data, motifSize=motifSize)

	with open(filename,'wb') as f:
		f.write(
		"\\documentclass{article}\n"
		"\\usepackage{amsmath,fullpage,graphicx,fancyhdr,xcolor,colortbl,chngpage}\n"
		"\\usepackage[landscape]{geometry}"
		"\\definecolor{yellow}{RGB}{255,255,70}\n"
		"\\definecolor{orange}{RGB}{255,165,70}\n"
		"\\definecolor{red}{RGB}{255,70,70}\n"
		"\\title{Motif Data}\n"
		"\\author{Graham Tremper}\n"
		"\\date{}\n"
		"\\fancyhead{}\n"
		"\\begin{document}\n"
		)

		statistics = {}		
		for corr in ('corr','lcorr','lacorr'):
			print "Starting " + corr +"..."
			motifsNL = findMotifs(data, ('NL',corr), motifSize = motifSize, degree=degree)
			motifsMCI = findMotifs(data, ('MCI',corr), motifSize = motifSize, degree=degree)
			motifsAD = findMotifs(data, ('AD',corr), motifSize = motifSize, degree=degree)
			motifsCONVERT = findMotifs(data, ('CONVERT',corr), motifSize = motifSize, degree=degree)

			motifsNLRAND = shufData[('NL',corr)]
			motifsMCIRAND = shufData[('MCI',corr)]
			motifsADRAND = shufData[('AD',corr)]
			motifsCONVERTRAND = shufData[('CONVERT',corr)]

			NLRANDkeys = set()
			for dic in motifsNLRAND:
				NLRANDkeys.update(set(dic.keys()))

			MCIRANDkeys = set()
			for dic in motifsMCIRAND:
				MCIRANDkeys.update(set(dic.keys()))

			ADRANDkeys = set()
			for dic in motifsADRAND:
				ADRANDkeys.update(set(dic.keys()))

			CONVERTRANDkeys = set()
			for dic in motifsCONVERTRAND:
				CONVERTRANDkeys.update(set(dic.keys()))

			allMotifs = list( motifsNL.keys
							& motifsAD.keys
							& motifsMCI.keys
							& motifsCONVERT.keys
							& NLRANDkeys
							& MCIRANDkeys
							& ADRANDkeys
							& CONVERTRANDkeys )

			allMotifs = motifsNL.topMotifs(30)

			motifStats = []
			for key in allMotifs:
				norm = motifsNL[key]
				mci = motifsMCI[key]
				ad = motifsAD[key]
				conv = motifsCONVERT[key]
				c1 = stats.ttest_ind(norm, mci)
				c2 = stats.ttest_ind(norm, ad)
				c3 = stats.ttest_ind(norm, conv)
				c4 = stats.ttest_ind(mci, ad)
				c5 = stats.ttest_ind(mci, conv)
				c6 = stats.ttest_ind(ad, conv)
				c7 = stats.ttest_ind(norm, [d[key] if key in d else 0. for d in motifsNLRAND])
				c8 = stats.ttest_ind(mci, [d[key] if key in d else 0. for d in motifsMCIRAND])
				c9 = stats.ttest_ind(ad, [d[key] if key in d else 0. for d in motifsADRAND])
				c10 = stats.ttest_ind(conv, [d[key] if key in d else 0. for d in motifsCONVERTRAND])
				motifStats.append((key,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10))

			motifStats.sort(key=lambda x: motifsNL[x[0]].mean(),reverse=True)

			f.write(
			"\\begin{table}[t]\n"
			"\\begin{adjustwidth}{-1.5in}{-1.5in} "
			"\\caption{Motif T-test results from "+corr+" data with using edge swap}\n"
			"\\centering\n"
			"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\n"
			"\\hline\n"
			"\\rowcolor[gray]{0.85}\n"
			"Key & NL to MCI & NL to AD & NL to Conv & MCI to AD & MCI to Conv & AD to Conv & NL to Rand & MCI to Rand & AD to Rand & Conv to Rand \\\\ \\hline\n"
			)
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
			"\\end{tabular}\n"
			"\\end{adjustwidth}"
			"\\end{table}\n"
			)

		f.write("\\end{document}\n")

	os.system("pdflatex -output-directory result " + filename)
	os.system("rm result/*.log result/*.aux")

def PDFdiststats(data, filename, edgeSwap=False, motifSize=3, degree=10):
	"""Output a latex pdf of motif distribution (entrophy, gini coeff, fatness) stats"""
	filename = "result/" + filename + ".tex"

	if not edgeSwap:
		motifsNLRAND = motifsMCIRAND = motifsADRAND = motifsCONVERTRAND = findMotifs(data,"rand",motifSize=motifSize,degree=degree)

	with open(filename,'wb') as f:
		f.write(
		"\\documentclass{article}\n"
		"\\usepackage{amsmath,fullpage,graphicx,fancyhdr,xcolor,colortbl,chngpage}\n"
		"\\usepackage[landscape]{geometry}"
		"\\definecolor{yellow}{RGB}{255,255,70}\n"
		"\\definecolor{orange}{RGB}{255,165,70}\n"
		"\\definecolor{red}{RGB}{255,70,70}\n"
		"\\title{Motif Data}\n"
		"\\author{Graham Tremper}\n"
		"\\date{}\n"
		"\\fancyhead{}\n"
		"\\begin{document}\n"
		)

		if edgeSwap:
			with open("SwapData"+str(degree)+".pkl","rb") as pic:
				randGraphs = pickle.load(pic)

		statistics = {}		
		for corr in ('corr','lcorr','lacorr'):
			print "Starting " + corr +"..."
			motifsNL = findMotifs(data, ('NL',corr), motifSize = motifSize, degree=degree)
			NLd = diststats(motifsNL)
			motifsMCI = findMotifs(data, ('MCI',corr), motifSize = motifSize, degree=degree)
			MCId = diststats(motifsMCI)
			motifsAD = findMotifs(data, ('AD',corr), motifSize = motifSize, degree=degree)
			ADd = diststats(motifsAD)
			motifsCONVERT = findMotifs(data, ('CONVERT',corr), motifSize = motifSize, degree=degree)
			CONVERTd = diststats(motifsCONVERT)
			if edgeSwap:
				motifsNLRAND = findMotifs(data, ('NL',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				NLRANDd = diststats(motifsNLRAND)
				motifsMCIRAND = findMotifs(data, ('MCI',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				MCIRANDd = diststats(motifsMCIRAND)
				motifsADRAND = findMotifs(data, ('AD',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				ADRANDd = diststats(motifsADRAND)
				motifsCONVERTRAND = findMotifs(data, ('CONVERT',corr), motifSize = motifSize, degree=degree, randGraphs=randGraphs)
				CONVERTRANDd = diststats(motifsCONVERTRAND)

			motifStats = []
			for pos,key in enumerate(('Entrophy', 'Gini Coeff', 'Fatness')):
				c1 = stats.ttest_ind(NLd[pos], MCId[pos])
				c2 = stats.ttest_ind(NLd[pos], ADd[pos])
				c3 = stats.ttest_ind(NLd[pos], CONVERTd[pos])
				c4 = stats.ttest_ind(MCId[pos], ADd[pos])
				c5 = stats.ttest_ind(MCId[pos], CONVERTd[pos])
				c6 = stats.ttest_ind(ADd[pos], CONVERTd[pos])
				c7 = stats.ttest_ind(NLd[pos], NLRANDd[pos])
				c8 = stats.ttest_ind(MCId[pos], MCIRANDd[pos])
				c9 = stats.ttest_ind(ADd[pos], ADRANDd[pos])
				c10 = stats.ttest_ind(CONVERTd[pos], CONVERTRANDd[pos])
				motifStats.append((key,c1,c2,c3,c4,c5,c6))#,c7,c8,c9,c10))

			f.write(
			"\\begin{table}[t]\n"
			"\\begin{adjustwidth}{-2in}{-2in} "
			"\\caption{Motif Distribution T-test results from "+corr+" data with using edge swap}\n"
			"\\centering\n"
			"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\n"
			"\\hline\n"
			"\\rowcolor[gray]{0.85}\n"
			"Measure & NL to MCI & NL to AD & NL to Conv & MCI to AD & MCI to Conv & AD to Conv & NL to Rand & MCI to Rand & AD to Rand & Conv to Rand \\\\ \\hline\n"
			)
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
			"\\end{tabular}\n"
			"\\end{adjustwidth}"
			"\\end{table}\n"
			)

		f.write("\\end{document}\n")

	os.system("pdflatex -output-directory result " + filename)
	os.system("rm result/*.log result/*.aux")
	
"""""""""""""""""HELPER FUNCTIONS"""""""""""""""""""""
	

def diststats(graphdict):
	""" Helper Function to calculate entrophy, gini coeff, fatness"""
	listofentrophy = []
	listofgini = []
	listoffatness = []
	for graph in graphdict:
		listofentrophy.append(findentrophy(graph))
		listofgini.append(findgini(graph))
		listoffatness.append(findfatness(graph))
	listofentrophy = np.array(listofentrophy)
	listofgini = np.array(listofgini)
	listoffatness = np.array(listoffatness) 
	return	(listofentrophy, listofgini, listoffatness)

def findentrophy(x):
	sum = 0
	for value in x:
		sum += math.log(value) * value
	return -sum

def findgini(x):
	N = len(x)
	B = sum( xi * (N-i) for i,xi in enumerate(x) ) / (N*sum(x))
	return 1 + (1./N) - 2*B

def findfatness(x):
	x.sort(reverse=True)
	N = min(int(len(x)/5), 1)
	return sum(x[:N])/sum(x[N:])


def main():
	pass


if __name__ == '__main__':
	main()

