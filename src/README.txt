Python wrapper for Kavosh motif finding algorithm

There's very little documentation for this so the following is just what I've discovered from the code and the example files. I've modified the code so that the input is much simpler. To run, simply run ./Kavosh sizeOfMotif. It will use the edgelist titled result/OUTPUT.txt. Our python code creates this edgelist in the findMotif() function

------ Compile -------
just run "make" in this directory. An executable called "Kavosh" is created


------- Kavosh Format ---------
These are commented out in the c++ code because they're unnecessary for our purposes

There are a few flags to put after "./Kavosh" in order to run it correctly

-h Display help
-i Input filename
-o Output path
-r Number of random graphs
-s Motif size


------- File Format -------
first line: number of nodes
remaining lines: List of edges. Source first, then space, then target

Example: A 5 node cycle graph is represented by the following
5
1 2
2 3
3 4
4 5
5 1


------ Running with a networkx graph ------
Since this is just kind of hacked together, you'll have to put Kavosh executable with all of your python files.

findMotifs in GraphParse.py takes a networkx graph G and a motifSize and parses and runs that graph on the C++ code



-------- OUTPUT --------

1 filee are outputted

MotifCount.txt
	The first line is the total number of subgraphs explored. What follows is a list of motif IDs and total number of that motif. The findMotif function reads this in to python so you never need to see it. Just make sure you have a "result" directory becasue that is the temparary file directory for findMotif().
	





