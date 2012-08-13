#include "randomGenerator.h"

generator::generator() {}

int generator::binarySearch(int *E, int v, int w, int l, int h) {
	int m;
	while(h >= l) {
		m = (l+h)/2;
		if(E[m] == w)
			return m;
		else {
			if(E[m] < w)
				h = m - 1;
			else
				l = m + 1;
		}
	}
	return m;
}

void generator::genRandGraph_Edge(Graph * g) {
	register int i, j;
	int len = g->Size();
	int p, q, a, b, c, d, c_ind, d_ind;
	
	int *Na, *Nb;
	int success = 0;

	//printf("Len = %d\n", len);
	for(j = 0; j < numOFexchange; j++) {
		for (a = 1; a <= len; a++) {
			Na = g->getNeighbours(a);
			do {
				b = g->get_vertex();	
			} while (a == b);
			
			Nb = g->getNeighbours(b);
			//printf("++++ in a ... b = %d\n", b);
			//fflush(stdout);
			//printf("++++ in a ... Nb[0] = %d\n", Nb[0]);
			
			p = q = -1;
			
			for (i = 0; i < numOFexchange; i++) {
				c_ind = (rand() % Na[0]) + 1;
				if(g->isConnected(a, Na[c_ind]) && !g->isConnected(b, Na[c_ind]) && Na[c_ind] != b && !g->isConnected(Na[c_ind], a)) {
					c = Na[c_ind];
					p = c_ind;
					break;
				}
			}
			
			if (p == -1)
				continue;
			
			for (i = 0; i < numOFexchange; i++) {
				//printf("---> Size = %d\n", Nb[0]); 
				//fflush(stdout);
				d_ind = (rand() % Nb[0]) + 1;
				if(g->isConnected(b, Nb[d_ind])	&& !g->isConnected(a, Nb[d_ind]) && Nb[d_ind] != a && !g->isConnected(Nb[d_ind], b)) {
					d = Nb[d_ind];
					q = d_ind;
					break;
				}
			}
			
			if (q == -1)
				continue;
			
			success++;
			
			g->deleteEdgeAdjMat(a,c);
			g->deleteEdgeAdjMat(b,d);
			g->addEdgeAdjMat(a,d);
			g->addEdgeAdjMat(b,c);
			
			g->swapEdge(a, p, d);
			g->swapEdge(b, q, c);
		
			p = binarySearch(g->getNeighbours(d),d,b,1,g->getNeighbours(d)[0]);
			g->swapEdge(d, p, a);
			p = binarySearch(g->getNeighbours(c),c,a,1,g->getNeighbours(c)[0]);
			g->swapEdge(c, p, b);			
			
		}
	}
	
}
