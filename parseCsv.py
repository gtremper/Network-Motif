import numpy as np
import csv
import math
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import string
import copy

GRAPHSIZE = 88

def parse(name):
	myReader = csv.reader(open(name, 'rb'), delimiter=',')
	indexdict = {}
	for index, name in enumerate(next(myReader)):
		indexdict[name] = index
	pdd = {}
	for row in myReader:
		name = row[indexdict['SubjectCode']]
		time = row[indexdict['Timepoint']]
		patienttype = row[indexdict['DXCURREN']]	
		converted = row[indexdict['DXCONTYP']]
		
		if converted:
			patienttype = converted
		if 'Normal Control' in patienttype:
			patienttype = string.replace(patienttype, 'Normal Control', 'NL')
		
		if name not in pdd.keys():
			pdd[name] = (patienttype, {})
		pdd[name][1][time] = row[indexdict['RightHippocampus']:indexdict['RightInsula']+1] 
	return pdd


def createCorrMatrixes(measuredict):
	"Create graph adjacency matricies for measuredict"
	if not set(['1','3','5']).issubset(measuredict.keys()):
		return False
	
	r_e = []
	r_f = []
	for i in xrange(GRAPHSIZE):
		if measuredict['1'][i] == 'NA' or measuredict['3'][i] == 'NA':
			return False
            #r_e.append(-1.0)
		else:
			rate = float(measuredict['1'][i])/float(measuredict['3'][i])
			rate = max(min(rate,2),0.5)
			r_e.append(rate)
		if measuredict['3'][i] == 'NA' or measuredict['5'][i] == 'NA':
			return False
            #r_f.append(-1.0)
		else:
			rate = float(measuredict['3'][i])/float(measuredict['5'][i])
			rate = max(min(rate,2),0.5)
			r_f.append(rate)
			
	corrmat = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	for i, early in enumerate(r_e):
		for j, final in enumerate(r_f):
			if i == j or early == -1.0 or final == -1.0:
				corrmat[i][j] = 0
			else:
				corrmat[i][j] = 1 - abs(early - final)/float(early + final)
				if corrmat[i][j] < 0:
					print "r_e[i], r_f[j]"
					
	boolean = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	for i, early in enumerate(r_e):
		for j in xrange(GRAPHSIZE):
			if early > 1.0:
				boolean[i][j] = 1.0
	lcorrmat = np.multiply(corrmat, boolean)
	
	boolean = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	for i in xrange(GRAPHSIZE):
		for j, rates in enumerate(izip(r_e, r_f)):
			early, final = rates
			if final > early:
				boolean[i][j] = 1.0
	lacorrmat = np.multiply(lcorrmat, boolean)
	
	return (corrmat, lcorrmat, lacorrmat)

def tresh(G):
    sorted_weights = np.sort(G, axis=None)
    threshold = sorted_weights[-GRAPHSIZE*10-1]
    return np.matrix(G>threshold)
    
		
def createFinalData():
    data = {}
    pdd = parse('outVolumes2.csv')
    for patientname, patientdata in pdd.iteritems():
        patienttype = patientdata[0]
        if patienttype in ['NL', 'AD', 'MCI', 'MCI to AD']:
            if patienttype == 'MCI to AD':
                patienttype = 'CONVERT'
        else:
            continue
        corrmatrixes = createCorrMatrixes(patientdata[1])
        if corrmatrixes == False:
            continue
        tcorrmatrixes = [tresh(x) for x in corrmatrixes]
#            for i in xrange(GRAPHSIZE):
#                for j in xrange(GRAPHSIZE):
#                    if matrix[i][j] > threshold:
#                        matrix[i][j] = 1
#                    else:
#                        matrix[i][j] = 0
        matrixdict = {}
        corr, lcorr, lacorr = corrmatrixes
        tcorr, tlcorr, tlacorr = tcorrmatrixes
        matrixdict['corr'] = corr
        matrixdict['lcorr'] = lcorr
        matrixdict['lacorr'] = lacorr
        matrixdict['tcorr'] = tcorr
        matrixdict['tlcorr'] = tlcorr
        matrixdict['tlacorr'] = tlacorr
        data[patientname] = ({}, matrixdict, patienttype)
        
    with open('aznorbert_corrsd_new.pkl','wb') as f:
        pickle.dump(data,f)
    
    return data

                               
def createData():
	data = defaultdict(list)
	pdd = parse('outVolumes2.csv')
	for patienttype, patientdata in pdd.values():
		corrmatrixes = createCorrMatrixes(patientdata)
		if corrmatrixes:
			c, l, la = corrmatrixes
			data[(patienttype, 'corr')].append(c)
			data[(patienttype, 'lcorr')].append(l)
			data[(patienttype, 'lacorr')].append(la)
	
	totalsum = 0
	for key in data.keys():
		totalsum += len(data[key])
		print str(key) + " " + str(len(data[key]))
	print 'Total Patients: ' + str(totalsum/3)
	
	with open('aznorbert_corrsd_new.pkl','wb') as f:
		pickle.dump(data,f)

if __name__ == '__main__':
	createData()