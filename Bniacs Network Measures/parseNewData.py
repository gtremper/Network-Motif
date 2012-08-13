#parseNewData.py

#parses the new Data for MS patients.

from corrData import *
import networkx as nx
from mynetalgs import myallmeasures
from parseHelper import normalStats,mapNormalStats,avgMeasuretoList
from corrData import *
from nalz_test import *

TYPES = ['corr','lcorr','lacorr','pcorr','pacorr']
GROUPTYPES = ["PP","CIS","RR","SP"]
FNAME = "ms.pkl"

#input: dictionary
#output: dictionary with each string in TYPES as key and a 3 list- list of corrs (for each type)
def filterDict(dic,groups = GROUPTYPES): #change groups to the three types of groups
    new_dic = dict()
    for t in TYPES:
        print(t)
        new_dic[t] = [dic[(g,t)] for g in groups]
    return new_dic
    
    
#input a filtered dictionary (the output of filterDict), it writes files into format that
#can be passed into nalz_test.


#NOTE: The filter_dic has key: 'corr' 'lcorr', etc
#                     and value: [corr1,corr2,corr3...corrN] (or [lcorr1,lcorr2,...])
# 10,50,236,30
def ConvertFormat(filter_dic):
    print("start computing stats for original (Directed) ")
    for t in TYPES:
        print("starting mapNormalStats on {}".format(t))
        f6 = open("ms_stats_{0}_directed.pkl".format(t),"wb")
        dump(mapNormalStats(filter_dic[t]),f6)
        f6.close()
        
        f7 = open("ms_stats_{0}_directed.pkl".format(t),"rb")
        normal_d = load(f7)
    
        print("starting avgMeasuretoList")
        f8 = open("ms_convertedFormat_{0}_directed.pkl".format(t),"wb")
        dump([avgMeasuretoList(d)for d in normal_d],f8)
        f8.close()
        
        print("finished {} (Directed)".format(t))
        
    print("start computing stats for original (unDirected) ")
    for t in TYPES:
        print("starting mapNormalStats on {}".format(t))
        f6 = open("ms_stats_{0}_undirected.pkl".format(t),"wb")
        dump(mapNormalStats(filter_dic[t],10.0,False),f6)
        f6.close()
        
        f7 = open("ms_stats_{0}_undirected.pkl".format(t),"rb")
        normal_d = load(f7)
    
        print("starting avgMeasuretoList")
        f8 = open("ms_convertedFormat_{0}_undirected.pkl".format(t),"wb")
        dump([avgMeasuretoList(d)for d in normal_d],f8)
        f8.close()
        
        print("finished {}".format(t))
        
def getFilteredDict(d):
    return filterDict(d)

#converts original dictionary into files for t-testing chart.
def main1():
    newD = getFilteredDict(load(open(FNAME,"rb")))
    print("file opening")
    file = open("convertedDict.pkl","wb")
    dump(newD,file)
    print("dumped contents in file and closed.")
    file.close()
    
    print("starting ConvertFormat process")
    file = open("convertedDict.pkl","rb")
    new_dict = load(file)
    ConvertFormat(new_dict)
    file.close()
    print("done with ConvertFormat Process")
    
    
#creates and saves t-testing chart in same directory.    
def main2():
    
    
    print("creating t-test charts")
    
    for t in TYPES:
        r = load(open("prand_D10.pkl","rb"))
        f = open("ms_convertedFormat_{0}_directed.pkl".format(t),"rb") #directed
        data = load(f)
        f.close()
        mycompare(data,GROUPTYPES,r,"t-test_{0}_directed".format(t))
        
        r = load(open("uprand_D10.pkl","rb"))
        f = open("ms_convertedFormat_{0}_undirected.pkl".format(t),"rb") #undirected
        data = load(f)
        f.close()
        mycompare(data,GROUPTYPES,r,"t-test_{0}_undirected".format(t))
        
    print("done with creating t-test charts.")
    
if __name__ != "__main__": 
    print("starting main1...") # when we finish main1(), you can comment it out
    #main1()                    # and proceed with main2 (if main1 files already created)
    print("main1 finished.")
    print("starting main2...")
    #main2()
    print("main2 finished.")
    
        
        
        

        
        
    
    
    
    
    
                 
            
    