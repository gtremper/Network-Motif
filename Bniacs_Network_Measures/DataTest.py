"""Test Pickle files (patient ID dictionary) to make sure they are the same"""

import numpy as np
import cPickle as pickle


#tests two Patient ID dictionaries and see if they match.
def testData():
	print "Opening first pickle file"
	with open("aznorbert_corrsd_new_measures.pkl","rb") as f: #first file name
		data1 = pickle.load(f)
	
	print "Opening second pickle file"                      
	with open("aznorbert_corrsd_new_measures.pkl","rb") as f: #second file name
		data2 = pickle.load(f)
		
	globalmax = 0
	
	print "Comparing..."
	for patient, tup1 in data1.iteritems():
		tup2 = data2[patient]
		measure1, mats1, type1 = tup1
		measure2, mats2, type2 = tup2
		if type1 != type2:
			print "Types wrong for "+str(patient)
			return False
		
		for key, value1 in mats1.iteritems():
			value2 = mats2[key]
			if not np.array_equal(value1, value2):
				print "Mats different for "+str(patient)
				return False
		
		for key, value1 in measure1.iteritems():
			value2 = measure2[key]
			
			if value1 is None:
				if value2 is None:
					continue
				else :
					print "One value is None"
					return False
					
			if value1 == 0.0:
				if value2 == 0.0:
					continue
				else :
					print "One value is zero"
					return False
					
			maxv = max(value1,value2)
			minv = min(value1,value2)
			globalmax = max(globalmax,maxv/minv)
	print "Globalmax: "+str(globalmax)
	return True

if testData():
	print "SUCCESS!"

			
		