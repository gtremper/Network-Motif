#mynetalgs.py
"""File for measures of graphs"""

import os
import sys
import time
import numpy as np
import scipy
import networkx as nx
import matplotlib.pyplot as plt
import cPickle 
import random
import pickle
from scipy import stats
from greedy import *


# import community - do this later
#from mpl_toolkits.mplot3d import Axes3D
#import mayavi.mlab as mmlab

'''  my wrappers to get single valued stats from clean graphs 

*** NOTE: some  routines expect a connected graph so clean (giant component)
first 

'''

#computes all measures    
def myallmeasures(G,directed=True ,test = 0):
    snames = [] #list of measures
    svals = []   #list of statistics values
    if directed:
        comps = nx.strongly_connected_components(G)
    else:
	    comps = nx.connected_component_subgraphs(G)
        
    gcomp = comps[0]  #largest components
    
    snames += ['StronglyConnected']
    svals += [ 1.0*len(gcomp)/G.number_of_nodes()]
    
    Ggiant = nx.subgraph(G,gcomp)
    
    svals += mydegrees(G)
    snames += ['avgoutdeg','stdoutdeg','fatoutdeg','avgindeg','stdindeg','fatindeg']
    
    if test > 0: #only do simplest measure
        return (snames,svals)
    #the rest of the measures
    
    svals += mycentrality(G)
    snames += ['avgcent','stdcent','fatcent']
    
    svals += mydiameter(Ggiant)
    snames += ['diameter']
    
    svals += myavgpathlength(Ggiant)
    snames += ['avgpathlength']
    
    svals += mypagerank(G)
    snames += ['stdpagerank','fatpagerank']
    
    svals += mylocalclustering(G)
    snames += ['avglocalclust','stdlocalclust','fatlocalclust']
    
    svals += myglobalclust(G)
    snames += ['globalclustering']

    svals += linkrank(G)
    snames += ['linkrank','LR groups']

    svals += modularity(G)
    snames += ['modularity','M groups']
    
    return (snames,svals)             
        

#given a dictionary compute the fatness (80-20 rule)
def fatness(d):
    d.sort()
    i80 = int(0.8*len(d))  #top 80% in size
    #fatness, if uniform then 1, if 80-20 rule then 3.2
    if np.sum(d) == 0:
        return None #if all values sum up to 0, returns None.
    else:
        fat = (1.0/0.2)*np.sum(d[i80:])/np.sum(d)
        return fat
    
#find the average, std and fattness of degree distribution
def mydegrees(G):
    if(nx.is_directed(G)):
        dd = G.out_degree()
        di = G.in_degree()
    else:
        dd = G.degree()
        di = dd;
    d = dd.values()
    i = di.values()
    avgdeg = np.average(d)
    stddeg = np.std(d)
    fatdeg = fatness(d)
    avgdegi = np.average(i)
    stddegi = np.std(i)
    fatdegi = fatness(i)
    return [avgdeg,stddeg,fatdeg,avgdegi,stddegi,fatdegi]

#find the avg, std, and fattness of the in_degree distrubition
def myindegrees(G):
    dd = G.in_degree()
    d = dd.values()
    avgdeg = np.average(d)
    stddeg = np.std(d)
    fatdeg = fatness(d)
    return [avgdeg,stddeg,fatdeg]

#find the avg, std, and fattness of the out_degree distrubition
def myoutdegrees(G):
    dd = G.out_degree()
    d = dd.values()
    avgdeg = np.average(d)
    stddeg = np.std(d)
    fatdeg = fatness(d)
    return [avgdeg,stddeg,fatdeg]

#local clustering (NOT directed triangles)
#using definition in : http://en.wikipedia.org/wiki/Clustering_coefficient
def mylocalclustering(G):
    if not(isinstance(G,nx.DiGraph)):
        dval = nx.clustering(G).values()
        return [np.average(dval),np.std(dval),fatness(dval)]
    uG = G.to_undirected()
    d = dict()
    for node in G.nodes():
        k = nx.degree(uG)[node] #doesn't count double edges twice!
        prod = 1.0*k*(k-1)
        l = []
        for e in G.edges():
			#if an edge's two nodes are neighbors of this node, then add to l.
            if e[0] in uG.neighbors(node) and e[1] in uG.neighbors(node):
                l+=[e]
        if prod != 0:
            d[node] = 1.0*len(l)/(prod) #count each edge individually.
        else:
            d[node] = 0.0 #no clustering of this node only has one or zero neighbors
    dval = d.values()
    avgdeg = np.average(dval)
    stddeg = np.std(dval)
    fatdeg = fatness(dval)
    return [avgdeg,stddeg,fatdeg]



#computes the global clustering coefficient for a directed graph
#@argument corr is the DiGraph
def myglobalclust(corr):
    if not isinstance(corr,nx.DiGraph):
        return [nx.transitivity(corr)]  #3*triangles/triads
    corr = np.array(nx.to_numpy_matrix(corr))
    mat = np.dot(corr,corr)
    paths = np.sum(mat) - np.trace(mat)
    mat = np.dot(mat,corr)
    loops = np.trace(mat)
    if paths == 0:
        return [0]
    else:
        return [float(loops)/paths]
    
        
#needs to fix
#find the average, std and fattness of centrality distribution
def mycentrality(G):
    dd = nx.betweenness_centrality(G)
    d = dd.values()
    avgcent = np.average(d)
    stdcent = np.std(d)
    fatcent = fatness(d)
    return [avgcent,stdcent,fatcent]

#find the average, std and fattness of rich club coefficients
def myrichclubcoefficient(G):
    dd = nx.rich_club_coefficient(G,normalized=False)
    d=dd.values()
    avgrich = np.average(d)
    stdrich = np.std(d)
    fatrich = fatness(d)
    return [avgrich,stdrich,fatrich]
        
    
#diameter of largest component - modified file in nx.distance_measures
def mydiameter(G):
    diameter = nx.diameter(G)
    return [diameter]
       
#average pathlength
#This function is usually called on Giant Component of graphs. If this is true, then except clause is never reached.
def myavgpathlength(G):
    try:
        apl =  nx.average_shortest_path_length(G)
        return [apl]
    except nx.NetworkXError as e:
        #this means graph is not connected
        if isinstance(G,nx.DiGraph):
		    return [nx.average_shortest_path_length(nx.strongly_connected_component_subgraphs(G)[0])]
        else:
            return [nx.average_shortest_path_length(nx.connected_component_subgraphs(G)[0])]
    except ZeroDivisionError as e:
        return [1]     

#average clustering (triangles) - what do to for Directed?
def myclustering(G):
    clust = nx.average_clustering(G)
    return [clust]
    
#smallworld - what do to for Directed?
def mysmallworld(G):
    smallworld = myclustering(G)[0]/myavgpathlength(G)[0]
    return [smallworld]
        
#pagerank, omitted the average page rank (wasn't useful.)
def mypagerank(G):
    dd=nx.pagerank_numpy(G)
    d = []
    for nd in G.nodes():
        d += [dd[nd]]
    avgpr = np.average(d)
    stdpr = np.std(d)
    fatpr = fatness(d)
    return [stdpr,fatpr]




