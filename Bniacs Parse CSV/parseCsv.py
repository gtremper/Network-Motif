import numpy as np
import csv
import math
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import string
import random
import copy
from scipy.stats import pearsonr

def parse_new(name, pdd={}):
	global GRAPHSIZE
	myReader = csv.reader(open(name, 'rb'), delimiter=',')
	indexdict = {}
	for index, name in enumerate(next(myReader)):
		if index >= 5:
			name = name.split('_')
			indexdict['_'.join(name[1:])] = index
		else:
			indexdict[name] = index
	for i, row in enumerate(myReader):
		if i == 0:
			GRAPHSIZE = 2*len(row[indexdict['superiorparietal_thickness']:indexdict['precuneus_thickness']+1])
			print "GRAPHSIZE: ", GRAPHSIZE
		raw = row[indexdict['freesurfer-id']].split('_')

		name = raw[0]
		if name in pdd.keys() and pdd[name] is None:
			continue
		time = int(raw[2])
		patienttype = row[indexdict['disease-subtype']]
		if len(patienttype) == 0:
			pdd[name] = None
			continue
#		converted = row[indexdict['DXCONTYP']]
#		
#		if converted:
#			patienttype = converted
#		if 'Normal Control' in patienttype:
#			patienttype = string.replace(patienttype, 'Normal Control', 'NL')
		
		if name not in pdd.keys():
			pdd[name] = (patienttype, defaultdict(list))
		pdd[name][1][time].extend(row[indexdict['superiorparietal_thickness']:indexdict['precuneus_thickness']+1])
	return pdd


def createCorrMatrixes_new(measuredict):
	"Create graph adjacency matricies for measuredict"
	if not set([1,2,3,4,5]).issubset(measuredict.keys()):
		return False
	
	r_e = []
	r_f = []
	for i in xrange(GRAPHSIZE):
		if 0 in [float(measuredict[k+1][i]) for k in xrange(5)]:
			return False
			
		rate = float(measuredict[1][i])/float(measuredict[3][i])
		rate = max(min(rate,2),0.5)
		r_e.append(rate)
		
		rate = float(measuredict[3][i])/float(measuredict[5][i])
		rate = max(min(rate,2),0.5)
		r_f.append(rate)
			
	corrmat = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	for i, early in enumerate(r_e):
		for j, final in enumerate(r_f):
			if i == j:
				continue
			else:
				corrmat[i][j] = 1 - abs(early - final)/float(early + final)
				if corrmat[i][j] < 0:
					print "This shouldn't be possible"
					
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
	
	pcorr = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	pacorr = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	
	for i in xrange(GRAPHSIZE):
		for j in xrange(GRAPHSIZE):
			if i == j:
				continue
			x = [float(measuredict[k+1][i]) for k in xrange(5)]
			y = [float(measuredict[k+1][j]) for k in xrange(5)]
			p,p2 = pearsonr(x,y)
			pcorr[i][j] = max(p, 0)
			pacorr[i][j] = abs(p)
	
	return (corrmat, lcorrmat, lacorrmat, pcorr, pacorr)

def createFinalDataNew(filename):
	data = defaultdict(list)
	pdd = parse_new('ms-lh-thk-5-tps-clinical.csv')
	pdd = parse_new('ms-rh-thk-5-tps-clinical.csv', pdd)
	
	
	for i, patient in enumerate(pdd.keys()):
		if pdd[patient] is None:
			print 'Patient Parse Fail:', i+1
			continue
		type = pdd[patient][0]
		corrmatrices = createCorrMatrixes_new(pdd[patient][1])
		if corrmatrices == False:
			print 'Patient Parse Fail:', i+1
			continue
		else:
			print 'Patient Parse Success:', i+1
			c, l, la, p, pa = corrmatrices
			tc, tl, tla, tp, tpa = [treshold(x) for x in corrmatrices]
		data[(type, 'corr')].append(c)
		data[(type, 'lcorr')].append(l)
		data[(type, 'lacorr')].append(la)
		data[(type, 'pcorr')].append(p)
		data[(type, 'pacorr')].append(pa)
		data[(type, 'tcorr')].append(tc)
		data[(type, 'tlcorr')].append(tl)
		data[(type, 'tlacorr')].append(tla)
		data[(type, 'tpcorr')].append(tp)
		data[(type, 'tpacorr')].append(tpa)
	
	with open(filename,'wb') as f:
		pickle.dump(data,f)
	
	for key in data.keys():
		print key, len(data[key])
		
	return data

def parse(name):
	global GRAPHSIZE
	myReader = csv.reader(open(name, 'rb'), delimiter=',')
	indexdict = {}
	for index, name in enumerate(next(myReader)):
		indexdict[name] = index
	pdd = {}
	
	for i, row in enumerate(myReader):
		if i == 0:
			GRAPHSIZE = len(row[indexdict['RightHippocampus']:indexdict['RightInsula']+1])
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
		if measuredict['1'][i] == 'NA' or measuredict['3'][i] == 'NA' or measuredict['5'][i]=='NA':
			return False
			
		rate = float(measuredict['1'][i])/float(measuredict['3'][i])
		rate = max(min(rate,2),0.5)
		r_e.append(rate)
		
		rate = float(measuredict['3'][i])/float(measuredict['5'][i])
		rate = max(min(rate,2),0.5)
		r_f.append(rate)
			
	corrmat = np.zeros(shape=(GRAPHSIZE, GRAPHSIZE))
	for i, early in enumerate(r_e):
		for j, final in enumerate(r_f):
			if i == j:
				continue
			else:
				corrmat[i][j] = 1 - abs(early - final)/float(early + final)
				if corrmat[i][j] < 0:
					print "This shouldn't be possible"
					
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


def treshold(G):
	sorted_weights = np.sort(G, axis=None)
	thresh = sorted_weights[-GRAPHSIZE*10-1]
	return np.matrix(G>thresh)
		
def createFinalData(undirected = False, filename="aznorbert_corrsd_new.pkl"):
	
	data = {}
	pdd = parse('outVolumes2.csv')
	for patientname, patientdata in pdd.iteritems():
		patienttype = patientdata[0]
		if patienttype in ('NL', 'AD', 'MCI', 'MCI to AD'):
			if patienttype == 'MCI to AD':
				patienttype = 'CONVERT'
		else:
			continue
		corrmatrixes = createCorrMatrixes(patientdata[1])
		if not corrmatrixes:
			continue
		if undirected:
			corrmatrixes = [x + x.T for x in corrmatrixes]
		tcorrmatrixes = [treshold(x) for x in corrmatrixes]
		corr, lcorr, lacorr = corrmatrixes
		tcorr, tlcorr, tlacorr = tcorrmatrixes
		matrixdict = {}
		matrixdict['corr'] = corr
		matrixdict['lcorr'] = lcorr
		matrixdict['lacorr'] = lacorr
		matrixdict['tcorr'] = tcorr
		matrixdict['tlcorr'] = tlcorr
		matrixdict['tlacorr'] = tlacorr
		data[patientname] = ({}, matrixdict, patienttype)
		
	with open(filename,'wb') as f:
		pickle.dump(data,f)
	
	return data

							   
def createData():
	data = defaultdict(list)
	pdd = parse('outVolumes2.csv')
	for patienttype, patientdata in pdd.values():
		if patienttype == "NL to AD":
			continue
		if patienttype == "NL to MCI":
			continue
		if patienttype == "MCI to AD":
			patienttype = "CONVERT"
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
	#createFinalData(False,"aznorbert_corrsd.pkl")
	createFinalDataNew('ms.pkl')





