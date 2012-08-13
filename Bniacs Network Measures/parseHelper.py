#parseHelper.py
"""Helper file for parseNewData.py"""
from corrData import *

#input: a graph
#output: computes the all measures of graph and zips the measure names and the measure values together.    
def zipMeasures(graph):
    stats = myallmeasures(graph)
    return zip(stats[0],stats[1])
    
    
#input: a list of stats
#output: ((measure strings),([stats for m1],[stats for m2]....[stats for mn])
def avgMeasuretoList(stats):
    temp = []
    MEASURES = []
    for i in range(len(stats[0])):
        MEASURES += [stats[0][i][0]]
        temp += [[s[i][1] for s in stats]]
    convertedStats = [tuple(MEASURES),tuple(temp)]
    return convertedStats
    
#normalStats computes the original stats of a list of corrs.
#input: a list of correlation matrices, t: threshhold (average degree of each graph), directed = True if graphs are directed, false otherwise.
#output: a list zipped measurements of the original thresholded matrices (without swaps)
def normalStats(corrs,t = 10.0,directed = True):
    n = 0
    s = []
    for c in corrs:
        if directed:
            print(n)
            n += 1
            g = nx.DiGraph(c > findThresh(c,t))
            s += [zipMeasures(g)]
        else:
            print(n)
            n += 1
            g = nx.Graph(c)
            ary = np.array(nx.to_numpy_matrix(g))
            newcorr = ary > findThresh(ary, 10.0)
            graph = nx.Graph(newcorr)
            s += [zipMeasures(graph)]
    return s
    
#mapNormalStats maps the normalStats function to a list of list of corrs.
def mapNormalStats(lst,t=10.0,directed=True):
    return [normalStats(e,t,directed) for e in lst]
    
    
    
    