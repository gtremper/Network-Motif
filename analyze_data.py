'''
Created on Jun 7, 2012

@author: jason
'''

import networkx as nx
import sys
import matplotlib.pyplot as plt
import cPickle
import copy
import alzTest

def minMaxValues(subdata):
  min = 1
  max = 0
  for j in range(len(subdata)):
    for k in range(len(subdata[j])):
      for l in range(len(subdata[j][k])):
        if subdata[j][k][l] < min and subdata[j][k][l] != 0:
          min = subdata[j][k][l]
        if subdata[j][k][l] > max and subdata[j][k][l] != 1:
          max = subdata[j][k][l]
  return (min, max)

def createGraphs(subdata, threshold):
  graphlist = []
  for j in range(len(subdata)):
    for k in range(len(subdata[j])):
      for l in range(len(subdata[j][k])):
        if subdata[j][k][l] > threshold:
          subdata[j][k][l] = 1
        else:
          subdata[j][k][l] = 0
    G = nx.DiGraph(subdata[j])
    graphlist.append(G)
  return graphlist

def createGraphs2(subdata, deg):
  avg = 0
  for j in range(len(subdata)):
    matrix = subdata[j]
    threshold = alzTest.findThresh(matrix, deg)
    for k in range(len(subdata[j])):
      for l in range(len(subdata[j][k])):
        if subdata[j][k][l] > threshold:
          subdata[j][k][l] = 1
        else:
          subdata[j][k][l] = 0
    G = nx.DiGraph(subdata[j])
    listoflists = nx.kosaraju_strongly_connected_components(G)
    avg += len(listoflists[0])
  return float(avg/len(subdata))

def determineThreshold(subdata, deg):
  sum = 0
  for matrix in subdata:
    sum += alzTest.findThresh(matrix, deg)
  return float(sum/len(subdata))
  
def degvssizeAvg(data):
  degreedict = {}
  for deg in range(1, 100, 5):
    subdata = copy.deepcopy(data[('NL', 'corr')])
    threshold = determineThreshold(subdata, deg)
    print "degree: " + str(deg) + " threshold: " + str(threshold)     
    myList = createGraphs(subdata, threshold)
    avg = 0
    for graph in myList:
      listoflists = nx.kosaraju_strongly_connected_components(graph)
      avg += len(listoflists[0])
    avg = float(avg/len(myList))
    degreedict[deg] = avg
    
  plt.plot(degreedict.keys(), degreedict.values(), 'bo')
  plt.show()

def degvssizeInd(data):
  degreedict = {}
  for deg in range(1, 100, 5):
    subdata = copy.deepcopy(data[('NL', 'corr')])
    print "degree: " + str(deg)
    avgSize = createGraphs2(subdata, deg)
    degreedict[deg] = avgSize
  
  plt.plot(degreedict.keys(), degreedict.values(), 'bo')
  plt.show()

if __name__ == '__main__':
    f = file("aznorbert_corrsd.pkl")
    myPickle = cPickle.Unpickler(f) 
    data = myPickle.load()
    degvssizeInd(data);
    
    
    
    
