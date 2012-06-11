from numpy import *
import networkx as nx


#returns threshold such that corr > threshold yeilds average deg of 'deg'
def findThresh(corr, deg):
    scorr = sort(ravel(corr))
    scorr = scorr[::-1]
    thresh = int(deg*88/2);
    return scorr[thresh]