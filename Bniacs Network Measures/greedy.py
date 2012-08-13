#greedy.py
from numpy import *
import networkx as nx
import random

def greedyMax(J,c,Q):
    B = copy(J)

    oldQ = -1
    groups = copy(c)
    old = zeros(len(c))
    while (False in (groups == old)) and oldQ != Q:
        oldQ = Q
        old = copy(groups)
        dc = zeros(len(c))
        while (False in (dc == c)):
            dc = copy(c)
            for i in range(len(c)):
                g = unique(c)
                dQ = 0
                switch = c[i]
                others = [j for j in g if j != c[i]]
                for j in others:
                    gain = c == j
                    loss = c == c[i]
                    loss[i] = False
                    dQt = sum(B[gain,i]) + sum(B[i,gain])
                    dQt = dQt - sum(B[loss,i]) - sum(B[i,loss])
                    if dQt > dQ + .000001:
                        dQ = dQt
                        switch = j
                Q += dQ
                c[i] = switch
        #print c,Q  
        clean = unique(c)
        for m in range(len(clean)):
            c[c == clean[m]] = m
        temp = copy(groups)
        for i in range(len(c)):
            groups[temp == i] = c[i]
        B2 = zeros((len(clean),len(clean)))
        
        for i in range(len(clean)):
            for j in range(len(clean)):
                temp = J[groups == i,:]
                temp = temp[:,groups == j]
                B2[i][j] = sum(temp)
        c = clean
        B = B2
    return [Q]

def linkrank(G):
    c = arange(len(G.nodes()))
    goo = nx.google_matrix(G)
    goo = array(goo)
    m = nx.pagerank_numpy(G)
    m = m.items()
    m = [i[1] for i in m]
    m = array([m])
    m = m.T
    L = tile(m,[1,len(goo)])*goo
    Q = 0
    mm = tile(m,[1,len(goo)])*tile(m.T,[len(goo),1])
    Qlr = L - mm
    return greedyMax(Qlr,c,0)

def modularity(G):
    c = arange(len(G.nodes()))
    A = array(nx.to_numpy_matrix(G))
    M = sum(A)
    outv = array([sum(A,axis=1)])
    outv = tile(outv.T,[1,len(A)])
    inv = array([sum(A,axis=0)])
    inv = tile(inv,[len(A),1])
    B = A - outv*inv/M
    J = B/M
    Q = sum(trace(J))
    return greedyMax(J,c,Q)