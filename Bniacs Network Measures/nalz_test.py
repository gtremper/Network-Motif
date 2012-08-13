from numpy import *
import networkx as nx
import matplotlib as mp
from cPickle import *
from mynetalgs import mydegrees, myallmeasures, myglobalclust, mypagerank
import scipy.stats as scistat
import scipy.special as sci
from random import sample, randint
import itertools as it

def GetData(filename = "aznorbert_corrsd_new.pkl"):
    global data #****
    global keys #****
    FILE = file(filename)
    data = Unpickler(FILE).load()
    keys = data.keys()
    #print keys

""" helper function for compare that returns a string
indicating whether data is less than or greater than """
def plusminus(input):
    if input < 0:
        return '(-)'
    return '(+)'

"""helper function for compare that returns the significance
of the p-value as a negative number"""
def importance(input):
    if input < .01:
        return -3   #invalid t-test values yield BLACK.
    elif input < .05:
        return -2
    elif input < .1:
        return -1
    return 0

"""helper function for compare that returns the appropriate color
for the appropriate significant level"""
def numColor(input):
    if input < 0.0:
        return (1,1,1)
    if input < .01:
        return (1,0,0)
    elif input < .05:
        return (1,.65,0)
    elif input < .1:
        return (1,1,0)
    return (1,1,1)

#formats an inputed number appropriate for a table in function compare
def numformat(num):
    if num < 0:
        return str("NAN")
    if num < 1.0/1000:
        return str(("%.02e"%num))
    return str(("%.04f"%num))
    
def RunNetAlgs():
    for key, value in data.iteritems():
        print key
        mats = value[1]
        corr = myallmeasures(nx.DiGraph(mats['tcorr']))
        lcorr = myallmeasures(nx.DiGraph(mats['tlcorr']))
        lacorr = myallmeasures(nx.DiGraph(mats['tlacorr']))
        for i in range(len(corr[0])):
            measure = corr[0][i]
            value[0][(measure,'corr')] = corr[1][i]
            value[0][(measure,'lcorr')] = lcorr[1][i]
            value[0][(measure,'lacorr')] = lacorr[1][i]

#picks out just the measure value dictionary in the data
def extractData():
    alldata = []
    for key, value in data.iteritems():
        alldata.append(value[0])
    dump(alldata, open("computedData_aznorbert.pkl","wb"))

#Takes the computed measured data in the entire data structure and
#converts it to the old format of the measure data. Includes all
#graph types for a patient type 'name': [corr,lcorr,lacorr]
#Also sorts by measure name, thus if all measures exist, the measures
#should be in the same order.
def convert(group,name):
    group = group.values()
    group = [i[0] for i in group if i[2] == name]
    corr = [[],[]]
    lcorr = [[],[]]
    lacorr = [[],[]]
    for elem in group:
        for key, val in elem.iteritems():
            if val == None:
                continue
            if key[1] == 'corr':
                if key[0] in corr[0]:
                    corr[1][corr[0].index(key[0])].append(val)
                else:
                    corr[0].append(key[0])
                    corr[1].append([val])

            if key[1] == 'lcorr':
                if key[0] in lcorr[0]:
                    lcorr[1][lcorr[0].index(key[0])].append(val)
                else:
                    lcorr[0].append(key[0])
                    lcorr[1].append([val])

            if key[1] == 'lacorr':
                if key[0] in lacorr[0]:
                    lacorr[1][lacorr[0].index(key[0])].append(val)
                else:
                    lacorr[0].append(key[0])
                    lacorr[1].append([val])

    corr = zip(corr[0],corr[1])
    lcorr = zip(lcorr[0],lcorr[1])
    lacorr = zip(lacorr[0],lacorr[1])
    g = lambda n: n[0]
    corr.sort(key = g)
    lcorr.sort(key = g)
    lacorr.sort(key = g)
    return [zip(*corr),zip(*lcorr),zip(*lacorr)]
            
""" -Prints the graph for a set of patients in a graph type
    -Data is input in format: [[measure list],[list of list of data]]
    must be in the implied order
    -Note, this uses a slightly modified version of table.py,
    resulting graph may not look pretty in some cases """    
#more generalized version of compare. Can input a list of groups and ran.        
def mycompare(lst = [],labels = ["A","B","C"],ran = '',title=''):
    subsets = list(it.combinations(list(range(len(lst))),2)) #subsets of 2
    WIDE = len(subsets) + len(labels) + 1 #the number of cells in a row
    
    ans1 = lst[0]
    mp.pyplot.figure()  
    rans1 = ran
    
    print ' '
    g = lambda x: x[-1]
    fin = []
    col = ['']
    col += ['{0} vs {1}'.format(labels[g1],labels[g2]) for g1,g2 in subsets]
    col += ['{} vs rand'.format(l) for l in labels]
    
    color = [[(1,1,1) for c in range(WIDE)] for r in range(len(ans1[0]))]
    for i in range(len(ans1[1])):
        print lst[0][0][i],lst[1][0][i],lst[2][0][i],lst[3][0][i], rans1[0][i]
        ans = []
        for g1,g2 in subsets:
            try:
                ans += [scistat.ttest_ind(lst[g1][1][i],lst[g2][1][i])]
            except:
                print("Encountered NAN, set to NEGATIVE")
                ans += [(0,-1)]
        for l in lst:
            try:
                ans += [scistat.ttest_ind(l[1][i],rans1[1][i])]
            except:
                print("Encountered NAN, set to NEGATIVE")
                ans += [(0,-1)]
            
        color[i] = [(.8,.8,.8)] #for measure axis
        [color[i].append(numColor(a[1])) for a in ans]
        imp = [importance(a[1]) for a in ans]
        imp = sum(imp)

        done = [ans1[0][i]]
        [done.append(numformat(a[1]) + plusminus(a[0])) for a in ans]
        done = done
        done.append(imp) #for sorting purposes
        color[i].append(imp) #for sorting purposes
        fin.append(done)

    fin.sort(key = g)
    color.sort(key = g)
    color = [i[0:WIDE] for i in color] #remove importance value
    fin = [i[0:WIDE] for i in fin] #remove importance value
    mp.pyplot.title(title)
    mp.pyplot.xticks([])
    mp.pyplot.yticks([])
    for a in fin:
        print a
    table = mp.pyplot.table(cellText = fin,loc = 'center',cellColours = color,
                            colLabels = col)
    mp.pyplot.subplots_adjust(bottom = .01, left = .01, right = .99, top = .9)
    k = table.properties()
    k = k.values()
    a = k[-1]
    b = a.keys()
    for n in b:
        if n[1] != 0:
            a[n].set_width(1.0/WIDE) #changed so cells will fit exactly on screen
            if n[0] == 0:
                a[n].set_facecolor((.8,.8,.8))
    fig = mp.pyplot.gcf()
    fig.set_size_inches(18.5,10.5)
    mp.pyplot.savefig(title) #saves the figure.

#if __name__ != '__main__':

    #d = load(open("ms_stats\ms_convertedFormat_corr_directed.pkl","rb")) #look in ms_stats for more files.
    #r = load(open("stats_randDirected.pkl","rb")) #for undireced, use "stats_randUnDirected"
    
    #mycompare(d,['NL','MCI','AD'],r,'T-TEST lacorr directed')

