from Bniacs_Motifs.FinalMotif import *
from Bniacs_Parse_CSV.parseCsv import *
from Bniacs_Network_Measures.nalz_test import *
from Bniacs_Network_Measures.parseNewData import *

def main():
  with open("aznorbert_corrsd_new.pkl","rb") as f:
    data = pickle.load(f)
  findMotifs(data,('AD', 'corr'),motifSize=3,degree=10,randGraphs=None, useCache=True, printMotifs=True)
  
if __name__ == '__main__':
  main()