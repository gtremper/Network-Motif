FANMOD command-line version 

input parameters (numbered in order):
algorithm options
 1- subgraph (motif) size [default(as appears in original FANMOD GUI): 3]
 2- # of samples used to determine approx. # of subgraphs [100000]
 3- full enumeration? 1(yes)/0(no) [1]
 sampling probabilities will be input at end of command line, as their number may vary!

input file
 4- infile name
 5- directed? 1(yes)/0(no) [1]
 6- colored vertices? 1(yes)/0(no) [0]
 7- colored edges? 1(yes)/0(no) [0]

random networks
 8- random type: 0(no regard)/1(global const)/2(local const) [2]
 9- regard vertex colors? 1(yes)/0(no) [0]
 10- regard edge colors? 1(yes)/0(no) [0]
 11- reestimate subgraph number? 1(yes)/0(no) [0]
 12- # of random networks [1000]
 13- # of exchanges per edge [3]
 14- # of exchange attempts per edge [3]

output file
 15- outfile name
 16- outfile format 1(ASCII - human readable)/0(CSV - for easy parsing) [1]
 17- create dumpfile? 1(yes)/0(no) [0]

sampling probabilities (if NOT full enumeration)
 18->25 (depending on subgraph size) [.5]

Each run also creates a log file in addition to the regular OUT file created by FANMOD.
This log file (<outputfilename>.log) includes run times for each analyzed network (original and randomized networks) or run-preventing errors.
A 'dump' file (<outputfilename>.dump) that contains a list of all subnetworks found in the input network as a list (including 
their adjacency matrix and the participating vertices) can also be created.	



October 2006
Tomer Benyamini
Yoav Teboulle
Tel-Aviv University