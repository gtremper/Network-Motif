import numpy as np
import csv
import math
import cPickle as pickle
from collections import defaultdict
import string

graphsize = 88

def parse(name):
    myReader = csv.reader(open(name, 'rb'), delimiter=',')
    indexdict = {}
    pdd = {}
    firstRow = True
    for row in myReader:
        if firstRow == True:
            for i in range(len(row)):
                indexdict[row[i]] = i
            firstRow = False
            continue
        
        name = row[indexdict['SubjectCode']]
        time = row[indexdict['Timepoint']]
        type = row[indexdict['DXCURREN']]
            
        convert = row[indexdict['DXCONTYP']]
        if convert != '':
            type = convert
        if 'Normal Control' in type:
            type = string.replace(type, 'Normal Control', 'NL')
        
        if name not in pdd.keys():
            pdd[name] = [type, {}]
        pdd[name][1][time] = row[indexdict['RightHippocampus']:indexdict['RightInsula']+1] 
    return pdd


def createCorrMatrixes(measuredict):
    if set(['1','3','5']).issubset(measuredict.keys()):
        r_e = []
        r_f = []
        for i in range(graphsize):
            if measuredict['1'][i] == 'NA' or measuredict['3'][i] == 'NA':
                r_e.append(-1.0)
            else:
                r_e.append(float(measuredict['1'][i])/float(measuredict['3'][i]))
            if measuredict['3'][i] == 'NA' or measuredict['5'][i] == 'NA':
                r_f.append(-1.0)
            else:
                r_f.append(float(measuredict['3'][i])/float(measuredict['5'][i]))
                
        corrmat = np.zeros(shape=(graphsize, graphsize))
        for i in range(graphsize):
            for j in range(graphsize):
                if r_e[i] == -1.0 or r_f[j] == -1.0:
                    corrmat[i][j] = 0
                else:
                    corrmat[i][j] = 1 - math.pow((r_e[i] - r_f[j]), 2)/100
                    if corrmat[i][j] < 0:
                        print r_e[i], r_f[j]
                        
        lcorrboolean = np.zeros(shape=(graphsize, graphsize))
        
        for i in range(graphsize):
            for j in range(graphsize):
                if r_e[i] > 1.0:
                    lcorrboolean[i][j] = 1.0
        lcorrmat = np.multiply(corrmat, lcorrboolean)
        
        lacorrboolean = np.zeros(shape=(graphsize, graphsize))
        
        for i in range(graphsize):
            for j in range(graphsize):
                if r_f[j] > r_e[j]:
                    lacorrboolean[i][j] = 1.0
        lacorrmat = np.multiply(lcorrmat, lacorrboolean)
        return (corrmat, lcorrmat, lacorrmat) 
    else:
        return None

        
        
def createData():
    data = defaultdict(list)
    pdd = parse('outVolumes2.csv')
    for element in pdd.values():
        corrmatrixes = createCorrMatrixes(element[1])
        if corrmatrixes != None:
            c, l, la = corrmatrixes
            data[(element[0], 'corr')].append(c)
            data[(element[0], 'lcorr')].append(l)
            data[(element[0], 'lacorr')].append(la)
    
    sum = 0
    for key in data.keys():
        sum += len(data[key])
        print str(key) + " " + str(len(data[key]))
    print 'Total Patients: ' + str(sum/3)
    
    with open('aznorbert_corrsd_new.pkl','wb') as f:
        pickle.dump(data,f)

if __name__ == '__main__':
    createData()