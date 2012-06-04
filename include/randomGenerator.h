#ifndef	RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#define numOFexchange 3

#include "graph.h"
#include <string>

using namespace std;

class generator {
	protected:
	
	public:
		generator();
		int binarySearch(int *E, int v, int w, int l, int h);
		void genRandGraph_Edge(Graph *g);
};
	
#endif //RANDOMGENERATOR_H
