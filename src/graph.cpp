#include "graph.h"
#include <iostream>

graph canon[MAXN * MAXM];
graph * nauty_g;
int* lab;
int* ptn;
int* orbits;


static DEFAULTOPTIONS(options);
statsblk(stats);
setword workspace[160*MAXM];

double * C_main;
double * C_rand;
double * mean;
double * var;
double* Score;
long * ID;
int idxID;
int head;
double enumerated_class;

FILE * o; 

extern bool isRand;
/****************************************************************
****************************************************************/

bool sortCmp(int i, int j) {
	return (j < i);
}

/****************************************************************
****************************************************************/

Graph::Graph(const int n, const int k) {
	register int i, j;
	subgraphSize = k;
	T = new tree(k);
	M = ((subgraphSize + WORDSIZE - 1) / WORDSIZE);
	nauty_g = new graph[subgraphSize * MAXM];
	lab = new int[subgraphSize];
	ptn = new int[subgraphSize];
	orbits = new int[subgraphSize];

	head = 1;
    nV = n;
	nE = 0;
	nEd = 0;
	E_temp.resize(nV+1);
	h = sizeof(Entry) << 3;

	options.writeautoms = FALSE;
	options.getcanon = TRUE;
	options.defaultptn = TRUE;
	options.digraph = TRUE;

	nauty_check(WORDSIZE, M, subgraphSize, NAUTYVERSIONID);
	
	rowSize = (int)ceil((float)nV/8);
	adjMat = new char[rowSize*(nV+1)+1];
	
	for(i = 1; i <= nV; i++) {
		for(j = 0; j <= (nV >> 3); j++) {
			adjMat[i*rowSize + j] = 0;
		}
	}	
}

/****************************************************************
****************************************************************/

Graph::~Graph() {
    //printf("Destroying graph ...\n");
	delete [] adjMat;
    for(int i = 1; i <= nV; i++) {
        delete [] E[i];    
    }
    delete [] E;
    delete [] lab;
    delete [] ptn	;
    delete [] orbits;
	delete [] nauty_g;
    delete [] ID;
    delete [] mean;
    delete [] var;
    delete [] Score;
    delete [] C_rand;
    delete [] C_main;
    delete [] degree;
 //   delete T;
    fclose(am);
}


/****************************************************************
****************************************************************/

void Graph::setPath(char *path) { 
    char file[256];
    sprintf(file, "%s/adjMatrix.txt", path);
    printf("Opening %s ...\n", file); 
    fflush(stdout);
	am = fopen(file, "w+");
    if(!am) {
        printf("Can't open %s\n", path);
        sprintf(file, "result/adjMatrix.txt", path);
        printf("Writing to %s instead ...\n", file); 
    	am = fopen(file, "w+");
        if (!am) {
            printf("Error again ... Sorry!\n");
            exit(-1);
        }
    }
}


/****************************************************************
****************************************************************/

void Graph::Print() {
	register int i, j;
	
	fflush(stdout);
	for(i = 1; i <= nV; i++) {
		printf("Node %d: (%d) ", i, E[i][0]);
		for(j = 1; j <= E[i][0]; j++) {
			printf("%d ", E[i][j]);
		}
		printf("\n");
	}
	printf("---------------------------------------------------------------\n");
	
}

/****************************************************************
****************************************************************/

void Graph::addEdgeAdjMat(vertex u, vertex v) {
	adjMat[(rowSize*u) + (v>>3)] |= (1 << (v % 8));
}

/****************************************************************
****************************************************************/

void Graph::deleteEdgeAdjMat(vertex u, vertex v) {
	adjMat[(rowSize*u) + (v>>3)] &= (~(1 << (v % 8)));
}

/****************************************************************
****************************************************************/

void Graph::swapEdge(vertex v, int ind, vertex u) {	
	vertex c = E[v][ind];
	
	if(u < E[v][ind]) {
		ind++;
		while(ind <= E[v][0] && u < E[v][ind]) {
			E[v][ind-1] = E[v][ind];
			ind++;
		}
		if(ind <= E[v][0] && u == E[v][ind]) {
			while(ind <= E[v][0]) {
				E[v][ind-1] = E[v][ind];
				ind++;
			}
			E[v][0]--;
		}
		else
			E[v][ind-1] = u;
	}
	else {
		ind--;
		while(ind > 0 && u > E[v][ind]) {
			E[v][ind+1] = E[v][ind];
			ind--;
		}
		if(ind > 0 && u == E[v][ind]) {
			ind++;
			ind++;
			while(ind <= E[v][0]) {
				E[v][ind-1] = E[v][ind];
				ind++;
			}
			E[v][0]--;
		}
		else
			E[v][ind+1] = u;
	}	
}

/****************************************************************
****************************************************************/

void Graph::addEdge(vertex u, vertex v) {
	nE++;
	E_temp[u].push_back(v);
	E_temp[v].push_back(u);
	adjMat[(rowSize*u) + (v>>3)] |= (1 << (v % 8));
}

/****************************************************************
****************************************************************/

bool Graph::isConnected(vertex u, vertex v) {
	return adjMat[(rowSize*u) + (v>>3)] & (1 << (v % 8));
}

/****************************************************************
****************************************************************/

int* Graph::getNeighbours(vertex v) {
	return E[v];
	
}

/****************************************************************
****************************************************************/

int Graph::get_vertex() {
	int ind = rand() % nEd;
	//printf("nEd = %d, ind = %d, degree = %d\n", nEd, ind, degree[ind]);
	return degree[ind];
	
}

/****************************************************************
****************************************************************/


void Graph::Finalize() {
	register int i, j, max = 0;
	vector<int>::iterator it;
	vector<int> degs;
	
	E = new int*[nV+1];
	degree = new int[nE];
	int degInd = 0;
	
	for(i = 1; i <= nV; i++) {
		sort(E_temp[i].begin(), E_temp[i].end(), sortCmp);
		it = unique(E_temp[i].begin(), E_temp[i].end());		
		E_temp[i].resize(it - E_temp[i].begin());
		
		degs.push_back(E_temp[i].size());

		if(max < E_temp[i].size())
			max = E_temp[i].size();
		
		E[i] = new int[E_temp[i].size() + 1];
		E[i][0] = E_temp[i].size();
		for(j = 0; j < E_temp[i].size(); j++) {
			E[i][j+1] = E_temp[i][j];
			if(isConnected(i, E_temp[i][j])) {
				degree[degInd] = i;
				degInd++;
				nEd++;
			}
		}
		E_temp[i].resize(0);
		E_temp[i].clear();		
	}
	maxDegree = max;
	/*
	for(int jj = 0; jj <degInd; jj++) {
		printf("Deg[%d] = %d\n", jj, degree[jj]);
	}
	
	int *temp = new int[max+1];
	for(int jj = 0; jj < max; jj++) {
		temp[jj] = 0;	
	}
	
	for(int jj =0; jj < degs.size(); jj++) {
		temp[degs[jj]]++;
	}

	for(int jj = 1; jj <= max; jj++) {
		printf("Degree = %d, Count = %d\n", jj, temp[jj]);
	}

    delete []temp;

	printf("Number of Nodes: %d\n", nV);
	printf("Number of Edges: %d\n", nE);
	printf("Maximum Degree: %d\n", maxDegree);
	*/
}

/****************************************************************
****************************************************************/

void Graph::Classify(vertex **subgraph, int level) {
	register int i = 0, j, l, k;
    set *gv;
	int tempSubgraph[subgraphSize];
	
	for (l = 0; l < level; l++) {
		for(k = 1; k <= subgraph[l][0]; k++) {
			tempSubgraph[i++] = subgraph[l][k];
		}
	}
	
	for (i = 0; i < subgraphSize; i++) {
		gv = GRAPHROW(nauty_g, i, M);
		EMPTYSET(gv, M);		
		for(j = 0; j < subgraphSize; j++) {
			if(i == j)
				continue;
			if (isConnected(tempSubgraph[i], tempSubgraph[j])) 
			{
				ADDELEMENT(gv, j);
			}
		}
	}
		
	nauty(nauty_g, lab, ptn, NULL, orbits, &options, &stats, 
		  workspace, 160*MAXM, M, subgraphSize, canon);
	
	T->init_cur_node();
	
	if (!isRand) {
		for (i = 0; i < subgraphSize-1; i++) {	
			for(j = 0; j < subgraphSize; j++) {
				if(i == j)
					continue;
				if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
					T->insert_one_main();
				else
					T->insert_zero_main();
			}
		}
		
		for(j = 0; j < subgraphSize-2; j++) {
			if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
				T->insert_one_main();
			else
				T->insert_zero_main();				
		}
		
		if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
			T->update_one_main(1);
		else 
			T->update_zero_main(1);
	}
	else {
		for (i = 0; i < subgraphSize-1; i++) {	
			for(j = 0; j < subgraphSize; j++) {
				if(i == j)
					continue;
				if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
					T->insert_one_rand();
				else
					T->insert_zero_rand();
			}
		}
		
		for(j = 0; j < subgraphSize-2; j++) {
			if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
				T->insert_one_rand();
			else
				T->insert_zero_rand();				
		}
		
		if(isConnected(tempSubgraph[lab[i]], tempSubgraph[lab[j]])) 
			T->update_one_rand(1);
		else 
			T->update_zero_rand(1);
	}
}

/****************************************************************
****************************************************************/

void Graph::AllocateCounter() {
	int class_num = T->get_leafnum();
	C_main = new double[class_num + 1];
	C_rand = new double[class_num];
	C_main[0] = 0;
	
	mean = new double[class_num];
	var = new double[class_num];
	Score = new double[class_num];
	
	register int i;
	for(i = 0; i < class_num; i++) {
		mean[i] = 0;
		var[i] = 0;
	}
	
	idxID = 0;
	ID = new long[class_num];
}

/****************************************************************
****************************************************************/

void Graph::DFS(Node * cur) {
	if(!cur->left && !cur->right) {
		Leaf * leaf = (Leaf *)cur;
		if(leaf->count > 0) {
			C_rand[head] = leaf->count;
			head++;
			leaf->count = 0;
		}
		return;
	}
	if(cur->left)
		DFS(cur->left);
	if(cur->right)
		DFS(cur->right);
} 

/****************************************************************
****************************************************************/

void Graph::print_adjMatrix(char * str) {
	register int i, j;
	int l = 0;
	int index = 0;
	int maxpow = subgraphSize * subgraphSize -1;

	for(i = 0; i < subgraphSize; i++) {
		for(j = 0; j < subgraphSize; j++) {
			if(i == j) {
				l++;
				//fprintf(am,"0");
			}
			else
			{
				index = i*(subgraphSize)+(j-l);
				//fprintf(am, "%d", str[index]);
				if (str[index] == 1)
					ID[idxID] += (long) (pow(2, (maxpow - (i*subgraphSize+j))));
			}
		}
		//fprintf(am,"\n");
	}
	//fprintf(am,"ID: %d", ID[idxID]);
	//fprintf(am,"\n\n");
	idxID++;
}

/****************************************************************
****************************************************************/

void Graph::DFSmain(Node * cur, char * str, int lev) {
	if(!cur->left && !cur->right) {
		print_adjMatrix(str);
		Leaf * leaf = (Leaf *)cur;
		head++;
		C_main[head] = leaf->count;
		leaf->count = 0;
		C_main[0]++;
		return;
	}
	if(cur->left) {
		str[lev] = 0;
		DFSmain(cur->left, str, lev+1);
	}
	if(cur->right) {
		str[lev] = 1;
		DFSmain(cur->right, str, lev+1);
	}
}

/****************************************************************
****************************************************************/

void Graph::Extract() {
	
	register int i, j;
	int class_num = T->get_leafnum();
	char * adj_str = new char[subgraphSize*(subgraphSize-1)];
	Node * current = T->return_root();
	
	head = 0;
	if(isRand) 
		DFS(current);
	
	else {
		DFSmain(current, adj_str, 0);
		enumerated_class = C_main[0];	
		//printf("Number of Non-isomorphic Classes: %.0f\n", enumerated_class);
	}
	
	if(isRand) {
		j = 0;
		for(i = 0; i < class_num; i++) {
			mean[j] += C_rand[i];
			var[j] += (C_rand[i]*C_rand[i]);
			j++;
		}
	}
}

/****************************************************************
****************************************************************/
void Graph::outputMotifResults(int subgraphCounter, char* path){
	FILE * cm;
	char file[256];
	sprintf(file, "%s/MotifCount.txt", path);
	cm = fopen(file, "w+");
	for(int i=0; i<T->get_leafnum(); i++){
		fprintf(cm,"%d %.0f %f\n",ID[i],C_main[i+1],C_main[i+1]/subgraphCounter);
	}
	fclose(cm);
}


void Graph::calculateZSCORE(int RAND, int subgraphCounter, char *path) {
	//cm2 is the same info as cm but in an easier to parse format
	FILE * cm;
	FILE * cm2;
	int i , j;
	for (i = 0; i < T->get_leafnum(); i++) {
		mean[i] = mean[i]/RAND;
		var[i] = sqrt((var[i]-(RAND*(mean[i]*mean[i])))/RAND);
		
		if(var[i] != 0)
			Score[i] = (C_main[i+1] - mean[i])/var[i];
		else
			Score[i] = -1;
	}
    
    char file[256];
	char file2[256];
    sprintf(file, "%s/ZScore.txt", path);
	sprintf(file2, "%s/ZS.txt", path);
    printf("Writing ZScores to %s ...\n", file); 
	cm = fopen(file, "w+");
	cm2 = fopen(file2,"w+");
    if(!cm) {
        printf("Can't open %s\n", path);
        sprintf(file, "result/ZScore.txt", path);
        printf("Writing to %s instead ...\n", file); 
    	cm = fopen(file, "w+");
        if (!cm) {
            printf("Error again ... Sorry!\n");
            exit(-1);
        }
    }

	fprintf(cm, "TOTAL NUMBER OF CLASSES:: %.0f\n\n", enumerated_class);
	fprintf(cm, "ID\t\tNUM IN REAL \t\t MEAN IN RANDOM \t VAR IN RANDOM \t\t ZSCORE\n");
	for (i = 0; i < T->get_leafnum(); i++) {
		if (var[i] != 0) 
			fprintf(cm, "%d\t\t %f%% \t\t %f%% \t\t %f \t\t %f\n", ID[i], C_main[i+1]/subgraphCounter*100, mean[i]/subgraphCounter*100, var[i], Score[i]);
		if (var[i] == 0) 
			fprintf(cm, "**%d\t\t %f%% \t\t %f%% \t\t %f \t\t %f \n", ID[i], C_main[i+1]/subgraphCounter*100, mean[i]/subgraphCounter*100, Score[i], double((C_main[i+1] - mean[i])));
	}
	fclose (cm);
	
	for(i=0; i<T->get_leafnum(); i++){
		if (var[i] != 0){
			fprintf(cm2,"%d %.3f\n",ID[i],Score[i]);
		}
	}
	fclose (cm2);
}
