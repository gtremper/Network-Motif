#ifndef GRAPH64_HPP
#define GRAPH64_HPP

#include <iostream>
using std::cout; using std::endl;
#include <algorithm>
using std::sort;
#include "nautyinit.h"

typedef unsigned long uint32;
typedef unsigned long long uint64;

typedef struct
{
    short matrix[64];
	short size;
	short g_N;
	short g_M;
	bool has_vertex_colors;
	bool has_edge_colors;
	unsigned short num_vertex_colors;
	unsigned short num_edge_colors;
	unsigned short num_vertex_bits;
	unsigned short num_edge_bits;
	bool directed;
	uint64 codestamp;
	graph nauty_canon[MAXN * MAXM];
	graph nauty_g[MAXN * MAXM];
	int lab[MAXN], ptn[MAXN], orbits[MAXN];
} graph64;

typedef uint64 graphcode64;

void init_graph(graph64 &g, short size, unsigned short num_vcolors, 
                unsigned short num_ecolors, bool directed);

graphcode64 toHashCode(graph64 &g); 

void readHashCode(graph64 &g, graphcode64 gc);

graphcode64 getGraphID(graph64 &g, graphcode64 gc); 

inline void delete_edge(graph64 &g, short source, short target) {
	g.matrix[(source<<3)|target] = 0; //g.matrix[u*8+v] = 0;
	DELELEMENT( ( GRAPHROW(g.nauty_g, source, g.g_M) ) , target);
}

inline void delete_vertex(graph64 &g, short idx) {
	g.matrix[(idx<<3)|idx] = 0;
}

inline void delete_element(graph64 &g, short row, short col) {
	g.matrix[(row<<3)|col] = 0;
	DELELEMENT( ( GRAPHROW(g.nauty_g, row, g.g_M) ) , col);
}

inline void set_edge(graph64 &g, short source, short target) {
	g.matrix[(source<<3)|target] = 1;
	ADDELEMENT( ( GRAPHROW(g.nauty_g, source, g.g_M) ) , target);
}

inline void set_edge(graph64 &g, short source, short target, short color) {
	g.matrix[(source<<3)|target] = color;
	ADDELEMENT( ( GRAPHROW(g.nauty_g, source, g.g_M) ) , target);
}

inline void set_element(graph64 &g, short source, short target, short color) {
	g.matrix[(source<<3)|target] = color;
	ADDELEMENT( ( GRAPHROW(g.nauty_g, source, g.g_M) ) , target);
}

inline void color_vertex(graph64 &g, short idx, short color) {
	g.matrix[(idx<<3)|idx] = color;
}

inline short get_element(const graph64 &g, const short row, const short col) {
	return g.matrix[(row<<3)|col];
}

#endif
