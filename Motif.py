#!/usr/bin/env python
# encoding: utf-8
"""
Motif.py

Created by Graham Tremper on 2012-07-10.
"""

import math
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess
import cPickle as pickle
import scipy.stats as stats
import graph_helper as gh
import random
from collections import defaultdict
from copy import deepcopy

USECACHE = False

def convertIDToGraph(id, motifSize, save=False):
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

def genRandMats(num):
	"""Generate random adjacency matricies"""
	mats = []
	for i in xrange(num):
		x = np.random.rand(88,88)
		x -= np.diag(np.diag(x))
		mats.append(x)
	return mats

def findMotifs(data,key,motifSize=3,degree=10,swap=False):
	"""Main finding motifs routine"""

	#Check cache
	filename = str(key)+'s'+str(int(motifSize))+'d'+str(int(degree))+".pkl"
	if os.path.exists('cache/'+filename) and USECACHE:
		print "in cache"
		with open('cache/'+filename,"rb") as f:
			return pickle.load(f)
	
	#generate random matricies
	if key == "rand":
		graphs = genRandMats(100)
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
		if swap:
			graph = gh.randomize_graph(graph, 1000)
		with open('result/OUTPUT.txt','wb') as f:
			f.write(str(len(graph)) + '\n')
			nx.write_edgelist(graph,f,data=False)
		
		#Jenky way to use c++ motif finder in python
		subprocess.call(['./Kavosh', str(motifSize)])
		data = np.loadtxt("result/MotifCount.txt",ndmin=2)
		for iD,total,percent in data:
			motifs[int(iD)].append(percent)
			
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

def plotMotifGraphs(data,motifSize=3,degree=10,numofmotifs=10,usetotal=False):
	"""Draws graph compairing average motif count between samples in the data"""
	for corr in ('corr','lcorr','lacorr'):

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

def PDFstats(data, filename, edgeSwap=False):
	"""Output a latex pdf of motif stats"""
	filename = "result/" + filename + ".tex"

	if not edgeSwap:
		motifsNLRAND = motifsMCIRAND = motifsADRAND = findMotifs(data,"rand")

	with open(filename,'wb') as f:
		f.write(
		"\\documentclass{article}\n"
		"\\usepackage{amsmath,fullpage,graphicx,fancyhdr,xcolor,colortbl}\n"
		"\\definecolor{yellow}{RGB}{255,255,70}\n"
		"\\definecolor{orange}{RGB}{255,165,70}\n"
		"\\definecolor{red}{RGB}{255,70,70}\n"
		"\\title{Motif Data}\n"
		"\\author{Graham Tremper}\n"
		"\\date{}\n"
		"\\fancyhead{}\n"
		"\\begin{document}\n"
		)
		for corr in ('corr','lcorr','lacorr'):
			motifsNL = findMotifs(data, ('NL',corr))
			motifsMCI = findMotifs(data, ('MCI',corr))
			motifsAD = findMotifs(data, ('AD',corr))
			if edgeSwap:
				motifsNLRAND = findMotifs(data, ('NL',corr), randomize=True)
				motifsMCIRAND = findMotifs(data, ('MCI',corr), randomize=True)
				motifsADRAND = findMotifs(data, ('AD',corr), randomize=True)

			allMotifs = list( set(motifsNL.keys())
							& set(motifsAD.keys())
							& set(motifsMCI.keys()) )

			motifStats = []
			for key in allMotifs:
				c1 = stats.ttest_ind(motifsMCI[key], motifsNL[key])
				c2 = stats.ttest_ind(motifsAD[key], motifsNL[key])
				c3 = stats.ttest_ind(motifsMCI[key], motifsAD[key])
				c4 = stats.ttest_ind(motifsNL[key], motifsNLRAND[key])
				c5 = stats.ttest_ind(motifsMCI[key], motifsMCIRAND[key])
				c6 = stats.ttest_ind(motifsAD[key], motifsADRAND[key])
				motifStats.append((key,c1,c2,c3,c4,c5,c6))

			motifStats.sort(key=lambda x: motifsNL[x[0]].mean(),reverse=True)

			f.write(
			"\\begin{table}[t]\n"
			"\\caption{Motif T-test results from "+corr+" data with using random matrices}\n"
			"\\centering\n"
			"\\vspace{2pt}\n"
			"\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
			"\\hline\n"
			"\\rowcolor[gray]{0.85}\n"
			"Key & MCI to Norm & AD to Norm & MCI to AD & Norm to Rand & MCI to Rand & AD to Rand \\\\ \\hline\n"
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
			"\\end{table}\n"
			)

		f.write("\\end{document}\n")

	subprocess.call(['pdflatex','-output-directory result',filename])
	subprocess.call(['rm','result/*.log','result/*.aux"'])
	#os.system("pdflatex -output-directory result " + filename)
	#os.system("rm result/*.log result/*.aux")


def main():
	with open("aznorbert_corrsd.pkl","rb") as f:
		data = pickle.load(f)	
	"""
	for i,G in enumerate(data[('NL','corr')]):
		print i
		sortedWeights = np.sort(G,axis=None)
		threshold = sortedWeights[-len(G)*10-1]

		graph = nx.DiGraph(G>threshold)
		graph = gh.randomize_graph(graph, 1000)
	"""
	
	PDFstats(data,"Motif_Statistics_Mats",False)


if __name__ == '__main__':
	main()
