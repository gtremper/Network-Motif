#corrData.py

#we will be given a dictionary with ID and other stuff. 
#corrData adds the measures to the dictionary for each patient.


#This file parses that dictionary and computes the measures for each graph.
#from corr import *
from cPickle import *
from mynetalgs import *

#findthresh
#input: a correlation matrix and desired average 
#output: a threshold value
def findThresh(corr, deg=10.0):
    scorr = np.sort(np.ravel(corr))
    scorr = scorr[::-1]
    thresh = int(deg*len(corr)); #this gets avg out degree deg.
    return scorr[thresh]

#input: dictionary: key = patient ID value = (dict for measures,dict for matrices,"patient type")
#default: the matrices are for directed graphs, set second arg to False if you want _undirected_ data.
#updates the dictionary so that it includes measures
def update_dict(d,directed = True):
    n = 1
    length = len(d)
    for key, value in d.iteritems():
        print("corrData.update_dict progress: {0} of {1}".format(n,d))
        if directed ==True:
            gr = nx.DiGraph
        else:
            gr = nx.Graph
        n += 1
        matrices = value[1]
        data_dict = value[0]
        measures, corrdata = myallmeasures(gr(matrices["tcorr"]))
        measures, lcorrdata = myallmeasures(gr(matrices["tlcorr"]))
        measures, lacorrdata = myallmeasures(gr(matrices["tlacorr"]))
        #measures, pcorrdata = myallmeasures(gr(matrices["tpcorr"]))
        #measures, pcorrdata = myallmeasures(gr(matrices["tpacorr"]))
        measures = zip(measures, corrdata, lcorrdata, lacorrdata) #,pcorrdata,pacorr)
        for m, c, lc , lac in measures:
            data_dict[(m,"corr")] = c
            data_dict[(m,"lcorr")] = lc
            data_dict[(m,"lacorr")] = lac
            #data_dict[(m,"pcorr")] = pc
            #data_dict[(m,"pacorr")] = pac
    print("update_dict function call finished.")

#main opens a file that contains a measure-less dictionary
#input: original file name(string), new filename (string)
#computes the measures and dumps a new dictionary onto a file w/ a different name.
def main(orig_fname = "aznorbert_corrsd_undirected.pkl", new_fname = "aznorbert_corrsd_new_measures_undirected.pkl"):
    with open(orig_fname,"rb") as f:#opens the file
        data = load(f)
    update_dict(data)
    with open(new_fname,"wb") as f: 
        dump(data,f)                #dumps the new updated dictionary into new_fname.
    print("corrData.main finished.")

#randomUndirected creates n randomly generated undirected graphs, computes measures, and returns the stats in a list.
#return value is a list of ([measure1,measure2,...,measure n],[stat1,stat2,...stat n])
def randomUndirected(n=100): 
    stats = []
    for i in range(n):
        x = np.random.rand(88,88)
        x = np.triu(x,1)
        x += x.T
        g = nx.Graph(x > findThresh(x))
        stats.append(myallmeasures(g))
        stat_dict = dict()
        for n in range(len(stats[0][0])):
            l = []
            for s in stats:
                l += [s[1][n]]
            stat_dict[stats[0][0][n]] = l
    return stats
   
#imports the .pkl file for dictionary ID.  
def importMS(f_name = ".pkl"):
    print("corrData: finished importing")
    d = load(open(f_name,"rb"))
    return d



   
    

