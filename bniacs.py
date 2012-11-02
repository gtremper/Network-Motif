from Bniacs_Motifs.FinalMotif import *
from Bniacs_Parse_CSV.parseCsv import *
from Bniacs_Network_Measures.nalz_test import *
from Bniacs_Network_Measures.parseNewData import *

def main():
  with open("aznorbert_corrsd_new.pkl","rb") as f:
    data = pickle.load(f)
#  print data.keys()
  for key in data.keys():
    for size in (3,4,5):
      findMotifs(data,key,motifSize=size,degree=10,randGraphs=None, useCache=True, printMotifs=True)
  
if __name__ == '__main__':
  main()