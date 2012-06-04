#ifndef GRAPH_H
#define GRAPH_H
#define MAXN 100

#include <vector>
#include <algorithm>
#include <math.h>
#include <string>
#include <nauty.h>
#include "ZeroOneTree.h"


using namespace std;

typedef int vertex;
typedef long long Entry;

class Graph {
public:
    Graph(const int n, const int k);
    ~Graph();
    void  setPath(char *path);
	void addEdgeAdjMat(vertex u, vertex v);
	void deleteEdgeAdjMat(vertex u, vertex v);
	void swapEdge(vertex v, int ind, vertex u);
	void  addEdge(vertex u, vertex v);
	int* getNeighbours(vertex v);
	bool isConnected(vertex u, vertex v);
	void Print();
	void Finalize();
	int Size() { return nV; }
	int MaxDegree() { return maxDegree; }
	void Classify(vertex ** subgraph, int level);
	void AllocateCounter();
	void DFS(Node * cur);
	void DFSmain(Node * cur, char * str, int lev);
	void Extract();
	void calculateZSCORE(int RAND, int subgraphCounter, char *path);
	void print_adjMatrix(char * str);
	int get_vertex();
	
private:
    FILE *am;
	vector< vector<vertex> > E_temp;
	int subgraphSize;
	int M;
	int *degree;
	int **E;
	int h;
	char *adjMat;
	int rowSize;
    int nV;
	int nE;
	int nEd;
	int maxDegree;
	tree *T;
};
#endif //GRAPH_H
