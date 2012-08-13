#sepGroups.py


from cPickle import *


#parses the patient ID dictionary, and seperates into 4 groups. (CONV,NL,AD,MCI)
def pType(d):
    conv,ad,mci,nl = [],[],[],[]
    for k,v in d.iteritems():
        if v[2] == "AD":
            ad.append(v[1])
        elif v[2] == "MCI":
            mci.append(v[1])
        elif v[2] == "NL":
            nl.append(v[1])
        else: #they are converts
            conv.append(v[1])
    return conv,ad,mci,nl

    
#parses each type of patient, seperates into three groups (corr,lcorr,lacorr)
def cType(pType):
    corr,lcorr,lacorr = [],[],[]
    for d in pType:
        corr.append(d["corr"])
        lcorr.append(d["lcorr"])
        lacorr.append(d["lacorr"])
    return [corr,lcorr,lacorr]
   

#saveFile saves data onto a pickle file named {label}.pkl.   
def saveFile(data,label):
    file = open("{}.pkl".format(label),"wb")
    dump(data,file)

    
#opens the dictionary file
#parses the dictionary and saves each group(ad,conv,mci,nl) into four files.
#for each particular group, there is a three element list.
#adCorrs[0] --> corr
#adCorrs[1] --> lcorr
#adCorrs[2] --> lacorr
#same goes for adCorrs
def sepGroups(dict_fname = "aznorbert_corrsd_new.pkl"):
    file = open(dict_fname,"rb")
    d = load(file)
    conv,ad,mci,nl = pType(d)
    convCorrs = cType(conv)
    adCorrs = cType(ad)
    mciCorrs = cType(mci)
    nlCorrs = cType(nl)
    saveFile(adCorrs,"adCorrs")
    saveFile(convCorrs,"convCorrs")
    saveFile(mciCorrs,"mciCorrs")
    saveFile(nlCorrs,"nlCorrs")
    print("done")
    
    
if __name__ != "__main__":
    # sepGroups("aznorbert_corrsd_new.pkl") -not in directory. 
   
    
    
    

