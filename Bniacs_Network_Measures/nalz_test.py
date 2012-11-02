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
for the appropriate significance level"""
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

#Runs the measures for mynetalgs on the data
#The data should be in the dictionary format
#described in the README
def RunNetAlgs(dumpData = True):
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
    if dumpData:
        dump(data,open('computed_data','wb'))
        dump(convert('AD'),open('final_AD.pkl','wb'))
        dump(convert('NL'),open('final_NL.pkl','wb'))
        dump(convert('MCI'),open('final_MCI.pkl','wb'))
        dump(convert('CONVERT'),open('final_CONVERT.pkl','wb'))

#Does the same as RunNetAlgs but can work on a range of degree
#thresholding for the graphs. Also will work do undirected
#Also dumps all the data but only in format specified for mycompare
def DegreeNetAlgs(undirected = False,startdeg = 10,enddeg = 10):
    for i in range(startdeg,enddeg+1):
        print i
        count = 0
        for key, value in data.iteritems():
            count += 1
            print count, i
            mats = value[1]

            corr = mats['corr']
            lcorr = mats['lcorr']
            lacorr = mats['lacorr']

            if undirected:
                corr += transpose(corr)
                lcorr += transpose(lcorr)
                lacorr += transpose(lacorr)
            
            corr = corr > findThresh(corr,i)
            lcorr = lcorr > findThresh(lcorr,i)
            lacorr = lacorr > findThresh(lacorr,i)
            
            corr = myallmeasures(nx.DiGraph(corr))
            lcorr = myallmeasures(nx.DiGraph(lcorr))
            lacorr = myallmeasures(nx.DiGraph(lacorr))
            for k in range(len(corr[0])):
                measure = corr[0][k]
                value[0][(measure,'corr')] = corr[1][k]
                value[0][(measure,'lcorr')] = lcorr[1][k]
                value[0][(measure,'lacorr')] = lacorr[1][k]
        dump(convert('AD',i),open('totalD_AD_D'+str(i)+'.pkl','wb'))
        dump(convert('NL',i),open('totalD_NL_D'+str(i)+'.pkl','wb'))
        dump(convert('MCI',i),open('totalD_MCI_D'+str(i)+'.pkl','wb'))
        dump(convert('CONVERT',i),open('totalD_CONVERT_D'+str(i)+'.pkl','wb'))

#Runs the measures on 100 random networks (specified direction by 'undirected')
#And dumps the data in the format needed for mycompare
def RandomDegrees(undirected = False,startdeg = 10, enddeg = 10):
    for des in range(startdeg,enddeg + 1):
        print des
        count = 0
        anslist = []
        orderedans = []
        for nn in range(100):
            x = random.rand(88,88)
            x -= diag(diag(x))
            if undirected:
                x = triu(x,1)
                x += x.T
            x = x > findThresh(x, des)
            patientGraph = nx.DiGraph(x)
            m = myallmeasures(patientGraph)
            anslist.append(m[1])
            if nn % 10 == 0:
                print nn, des
        mnames = m[0]; 
        for j in range(len(mnames)):
            measure = [i[j] for i in anslist]
            measure = [i for i in measure if i != None]
            orderedans.append(measure)    
        corr = zip(mnames,orderedans)
        g = lambda n: n[0]
        corr.sort(key = g)
        dump([mnames,orderedans], open("ufinal_rand_D"+str(des)+".pkl","wb"))
        
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
#Note: must have run 'GetData' at some point before running this
def convert(name, deg = 10):
    group = data.values()
    group = [i[0] for i in group if i[2] == name]
    corr = [[],[]]
    lcorr = [[],[]]
    lacorr = [[],[]]
    for elem in group:
        passcorr = 1
        passlcorr = 1
        passlacorr = 1
        #if one wants to exclude all graphs that couldn't
        #attain the avgdeg, "deg", uncomment the below code
##        if elem[('avgindeg','corr')] < deg:
##            passcorr = 0
##            print elem[('avgindeg','corr')]
##        if elem[('avgindeg','lcorr')] < deg:
##            passlcorr = 0
##            print elem[('avgindeg','lcorr')]
##        if elem[('avgindeg','lacorr')] < deg:
##            passlacorr = 0
##            print elem[('avgindeg','lacorr')]
        for key, val in elem.iteritems():
            if val == None:
                continue
            if passcorr and key[1] == 'corr':
                if key[0] in corr[0]:
                    corr[1][corr[0].index(key[0])].append(val)
                else:
                    corr[0].append(key[0])
                    corr[1].append([val])

            if passlcorr and key[1] == 'lcorr':
                if key[0] in lcorr[0]:
                    lcorr[1][lcorr[0].index(key[0])].append(val)
                else:
                    lcorr[0].append(key[0])
                    lcorr[1].append([val])

            if passlacorr and key[1] == 'lacorr':
                if key[0] in lacorr[0]:
                    lacorr[1][lacorr[0].index(key[0])].append(val)
                else:
                    lacorr[0].append(key[0])
                    lacorr[1].append([val])
    print corr[0][0]
    print len(corr[1][0])
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

