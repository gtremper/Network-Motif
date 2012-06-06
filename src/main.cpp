#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include <graph.h>
#include <getopt.h>
#include <randomGenerator.h>
#include <iostream>

#define EPOC 10
//#define Debug

using namespace std;

void GEN ( int n, int k, int root, int level, int reminder, int m);
void NEG ( int n, int k, int root, int level, int reminder, int m);

//1st dimension is the layer, while the 2nd is the list of selected nodes in each layer. subgraph[k][0] is the number of nodes in layer k.
//vertex subgraph[subgraphSize][subgraphSize+1];
vertex **subgraph;
int subgraphSize = -1, num_random_graphs = 0;

bool *Visited;

// This array collects the valid children at each layer.
unsigned int **childSet;

// This array maps the childSet to the subgraph[level]
unsigned int **Index;  

//total number of enumerated subgraphs.
long long subgraphCounter; 

//g stores the input graph
Graph *g;

long num;

bool isRand;


/****************************************************************
****************************************************************/

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: Kavosh options[inputfile...] \n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input filename\tInput filename.\n"
		 "\t-o	--output path\t\tOutput directory.\n"
		 "\t-r 	--random number \tNumber of random graphs (default = 0).\n"
		 "\t-s 	--size motifsize \tMotif size.\n" );
	     
    exit (exit_code);
}

/****************************************************************
****************************************************************/

bool ReadData(const char *path) {
	register int i, j;
	int graphSize;
	FILE * inFile = fopen(path, "r");
	
	if (inFile == NULL) {
		printf("Error opening %s file.\n", path);
		return false;
	}
	
	if (!feof(inFile))
		fscanf(inFile, "%d\n", &graphSize);	
	
	Visited = new bool[graphSize+1];
	for(i = 1; i <= graphSize; i++)
		Visited[i] = false;
	
	
	g = new Graph(graphSize, subgraphSize);
	while (!feof(inFile)) {
		fscanf(inFile, "%d %d\n", &i, &j);
		if(i == j) continue;
		g->addEdge(i, j);
	}
	
	g->Finalize();
	
	//g->Print();
	
	fclose(inFile);

	childSet = new unsigned int*[subgraphSize];
	Index    = new unsigned int*[subgraphSize];
	
	for(i = 0; i < subgraphSize; i++) {
		childSet[i] = new unsigned int[g->MaxDegree() * subgraphSize + 1];
		Index[i] = new unsigned int[g->MaxDegree() * subgraphSize + 1];
	}

	return true;
}

/***************************************************************
 * This function finds the valid children in each given level. *
***************************************************************/

void initChildSet(int root, int level) {
	register int *N;
	register int i, j;
	const int *parents = subgraph[level-1];
		
	childSet[level][0] = 0;
	for(i = 1; i <= parents[0]; i++) {
		N = g->getNeighbours(parents[i]);
		for(j = 1; j <= N[0] && root <= N[j]; j++) {
			if(!Visited[N[j]]) {
				Visited[N[j]] = true;
				childSet[level][0]++;
				childSet[level][childSet[level][0]] = N[j];				
			}
		}		
	}
}

/****************************************************************************************
 * This function Explores the root vertex and generates subgraphs containing that node. *
****************************************************************************************/

void Explore(vertex root, int level, int reminder) {
	register int i, j, k; // k is the number of selected nodes in the current level.

#ifdef Debug
	printf("************************************\n");
	printf("*****   Exploring level %3d  *******\n", level);
	printf("************************************\n");
#endif

	if (reminder == 0) { //reminder == 0 assures level <= subgraphSize
		subgraphCounter++;
		
#ifdef Debug
		printf("--> Subgraph Number %d: \n", subgraphCounter);
		for(i = 0; i < level; i++) {
			printf("Level %d: ", i);
			for(k = 1; k <= subgraph[i][0]; k++) {
				printf("%d ", subgraph[i][k]);
			}
			printf("\n");
		}
		printf("\n");
		printf("------------------------------------\n");
#endif

		g->Classify(subgraph, level);
		
		return;
	}
	
	
	initChildSet(root, level); 
	
#ifdef Debug
	printf("Valid Children in level %d:\nN = { ", level);
	for(k = 1; k <= childSet[level][0]; k++) {
		printf("%d ", childSet[level][k]);
	}		
	printf("}\n");
#endif

	for(k = 1; k <= reminder; k++) {
		if( childSet[level][0] < k ) { //There is not enough child to choose m from.
			for(i = 1; i <= childSet[level][0]; i++) {
				Visited[childSet[level][i]] = false;
			}
			return;
		}

#ifdef Debug
		printf("Selecting %d node(s) from level %d\n", k, level);
		printf("Initial Selection = { ");
		for(i = 1; i <= k; i++) {
			printf("%d ", childSet[level][i]);
		}
		printf("}\n");
#endif
		subgraph[level][0] = k;
		for(i = 1; i <= k; i++) {
			subgraph[level][i] = childSet[level][i];
			Index[level][i] = i;
		}	
		
		Explore(root, level + 1, reminder - k);
	    GEN( childSet[level][0], k, root, level, reminder, k);

#ifdef Debug
		printf("************************************\n");
		printf("*****    Back to level %3d   *******\n", level);
		printf("************************************\n");
#endif
	}
	
	for(i = 1; i <= childSet[level][0]; i++) {
		Visited[childSet[level][i]] = false;
	}
	subgraph[level][0] = 0;
	return;
}

/***************************************************************************************************
 * The following three functions generate all C(n, k) in Gray Code order, Adopted from Rusky code. *
***************************************************************************************************/

void swap( int i, int j, int root, int level, int reminder, int m) {
#ifdef Debug
	printf("Switch %d with %d in level %d\n", childSet[level][j], childSet[level][i], level);
#endif

	Index[level][i] = Index[level][j];
	subgraph[level][Index[level][i]] = childSet[level][i];	
	Explore(root, level + 1, reminder - m);
}

/****************************************************************
****************************************************************/

void GEN( int n, int k, int root, int level, int reminder, int m) {	
	if (k > 0 && k < n) {
    	GEN( n-1, k, root, level, reminder, m);
		if (k == 1) 
			swap( n, n-1, root, level, reminder, m);  
		else 
			swap( n, k-1, root, level, reminder, m);
    
		NEG( n-1, k-1, root, level, reminder, m);
    }
}

/****************************************************************
****************************************************************/

void NEG( int n, int k, int root, int level, int reminder, int m) {	
	if (k > 0 && k < n) {
    	GEN( n-1, k-1, root, level, reminder, m);
    	
		if (k == 1) 
			swap( n-1, n, root, level, reminder, m);  
		else 
			swap( k-1, n, root, level, reminder, m);
    
		NEG( n-1, k, root, level, reminder, m);
	}
	
}

/***********************************************************************************
 * This function enumerates the subgraphs related to each vertex of inpur network. *
***********************************************************************************/
void Enumerate() {
	register int v;

	for (v = 1; v <= g->Size(); v++)
	{
#ifdef Debug
		printf("+ Exploring Node %d ...\n", v);
#endif Debug
		subgraph[0][0] = 1;
		subgraph[0][1] = v;
		
		Visited[v] = true;
		Explore(v, 1, subgraphSize - 1);
		Visited[v] = false;
	}	
}

/****************************************************************
****************************************************************/

int main(int argc, char *argv[]) {
	double total_random_time = 0;
	double random_time , main_time;
	clock_t start_random_time, end_random_time;

	register int i, j;
	//long long subgraphCounterMain;
	generator gen;
	int next_option;
	const char *const short_options = "hi:o:r:s:";
 	const struct option long_options[] = {
		{"help",   0, NULL, 'h'},
		{"input",  1, NULL, 'i'},
		{"output", 1, NULL, 'o'},
		{"random", 1, NULL, 'r'},
		{"size",   1, NULL, 's'},		
		{NULL,     0, NULL,  0 }		
	};
	
	char *program_name;
    char input_filename[256], output_directory[256];

    int verbose = 0;
    strcpy(output_directory, "result");

    program_name = argv[0];
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_filename, optarg);
	    		break;

			case 'o':
				strcpy(output_directory, optarg);
	    		break;
			
			case 'r':
				 num_random_graphs = atoi(optarg);
	    		break;

			case 's':
				subgraphSize = atoi(optarg);
	    		break;
			
			case '?':
	    		print_usage (stdout, 1);
				
			case -1:		/* Done with options. */
			    break;
			
			default:		/* Something else: unexpected. */
                print_usage (stderr, -1);
		}
    } while (next_option != -1);


	if (input_filename == NULL) {
		fprintf(stderr, "Input Argument Error: Please specify an input filename using \"-i <FILE NAME>\".\n");
        print_usage (stderr, -1);
	}
	
	if (subgraphSize == -1) {
		fprintf(stderr, "Input Argument Error: Please specify a motif size using \"-s <MOTIF SIZE>\".\n");
        print_usage (stderr, -1);
	}

	subgraph = new int*[subgraphSize];
	for (int i = 0; i < subgraphSize; i++)
		subgraph[i] = new int[subgraphSize+1];
	
	num = 0;	
	printf("Motif Size: %d\n", subgraphSize);
	printf("Input Graph: %s\n", input_filename);
		
	if (!ReadData(input_filename))
		return 0;
	
    //g->setPath(output_directory);

	clock_t startTime = clock();
	//for main graph
	//isRand = false;
	subgraphCounter = 0;
	Enumerate();
	g->AllocateCounter();
	printf("Total Number of Subgraphs: %d\n", subgraphCounter);	
	g->Extract();
	
	clock_t end_main_time = clock();
	main_time = difftime(end_main_time, startTime)/(double)CLOCKS_PER_SEC;
	printf("Time Used for main graph: %f\n", main_time);
	
	//This function was added because we aren't using random graphs
	g->outputMotifResults(output_directory);
	
	/*
	subgraphCounterMain = subgraphCounter;
	//for random graphs
    srand(time(NULL));
	isRand = true;
    long long boz = 0;
	for (i = 1; i <= num_random_graphs; i++) {
		clog <<"\r"<<"Random Graphs Completed: "<< i << "/"<<num_random_graphs <<flush;
		//start_random_time = clock();
		gen.genRandGraph_Edge(g);
		subgraphCounter = 0;
		Enumerate();
		g->Extract();
		//printf("Total Number of Subgraphs in Random graph %d: %d\n", i, subgraphCounter);	
		//boz  += subgraphCounter;
		//end_random_time = clock();
		//random_time = difftime(end_random_time, start_random_time)/(double)CLOCKS_PER_SEC;
		//printf("Time Used for random %d : %f\n", i, random_time);
		//total_random_time += random_time;
	}
	printf("\n");
	
	//printf("Avg. RAndom = %lld\n", (long long)((float)boz/num_random_graphs));
	if (0 < num_random_graphs)
		g->calculateZSCORE(num_random_graphs, subgraphCounterMain, output_directory);
	*/
		
	for(i = 0; i < subgraphSize; i++) {
		delete [] Index[i];
		delete [] childSet[i];
	}
	
	for (int i = 0; i < subgraphSize; i++)
		delete [] subgraph[i];

    delete [] subgraph;
	delete [] Index;
	delete [] childSet;
    delete [] Visited;
	delete g;
	
	clock_t endTime = clock();
	double total_time = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;
	printf("Time Used: %f\n", total_time); 
	printf("=============================\n");	
	
	//printf("Time Used for random : %f\n", total_random_time); 		
	
	return 0;
}
