#doubleCheck.py

#checks for NAN/None values in the Patient ID file with stats
from cPickle import *
import numpy as np


def main(dict_fname = "aznorbert_corrsd_new_measures_undirected.pkl"):
    file = open(dict_fname, "rb")
    d = load(file)
    file.close()
    nones = []
    for k,v in d.iteritems():
        for k1,v1 in v[0].iteritems():
            if v1 == None:
                print(k,k1,v1)
                nones += [(k,k1,v1)]

#result yields empty lists if there are no none-values.
main()                
        
        
