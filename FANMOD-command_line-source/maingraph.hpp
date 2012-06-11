#ifndef MAINGRAPH_HPP
#define MAINGRAPH_HPP


#include <iostream>
using std::endl;
#include <stdexcept>
using std::domain_error;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
//using std::ofstream;
#include <iostream>
#include <sstream>
using std::stringstream;

extern "C" {
#include <ctime>
}

#include "random.hpp"
#include "hashmap.h"
#include "nautyinit.h"
#include "graph64.hpp"

typedef unsigned long long edge;
typedef unsigned long vertex;
typedef short edgetype;

const edgetype NOEDGE_UV = 0;
const edgetype DIR_U_T_V = 1;
const edgetype DIR_V_T_U = 2;
const edgetype UNDIR_U_V = 3;
const short INDEG = 0;
const short OUDEG = 1;
const vertex NILLVERTEX = 0xFFFFFFFFUL;
const edge NILLEDGE = 0xFFFFFFFFFFFFFFFFULL;
const short G_N_MAX = 8;

inline edge new_edge(vertex u, vertex v) { return uint64(u) << 32 | uint64(v); }

inline edge edge_code(vertex u, vertex v) 
{ return (u < v) ? (uint64(u)<<32 | uint64(v)) : (uint64(v)<<32 | uint64(u)); }

inline vertex edge_get_u(edge e) { return vertex(e >> 32); }

inline vertex edge_get_v(edge e) { return vertex(e & 0xFFFFFFFFULL); }

typedef struct
{
    hash_map < edge, edgetype > edges;
    unsigned long n, m; // Number vertices, Number edges
    unsigned long num_nodes;
    unsigned long num_lonely_nodes;
    unsigned long num_dir_edges;
    unsigned long num_undir_edges;
    bool directed;
    unsigned short num_vertex_colors;
    unsigned short num_edge_colors;
    
    short *vertex_colors;
    unsigned long *num_neighbours;
    vertex **neighbours;
//	unsigned long *num_larger_neighbours;
//	unsigned long maxnumlargerneighbors;
	unsigned long maxnumneighbours;	
    unsigned long *v_util;
    
}  maingraph;

typedef struct {
    graphcode64 gc;
    unsigned int vertices[8];
} subgraph;

inline subgraph get_subgraph(vertex* v_subgraph, short k, graphcode64 gc) {
    subgraph s;
    s.gc = gc;
    s.vertices[0] = 0;
    s.vertices[1] = 0;
    for (short i = 0; i!= k; ++i)
        s.vertices[i] = v_subgraph[i+1] ;
    return s;        
}

inline vertex getvertex(subgraph s, short i) {
    //vertex v = (s.vertices[i] >> ((i%4)*12)) & 0x000000000000FFFFULL;
    return s.vertices[i];
}

// These are constants for the random types
const short NO_REGARD = 0;    // No regard to bidirectional edges
const short GLOBAL_CONST = 1; // Keeping bidir edges globally const
const short LOCAL_CONST = 2;  // Keeping bidir edges const for each vertex

inline edgetype rm_colors(const edgetype input){
    return input & 3; // 3 = bit-mask 11
}

inline edgetype get_uv(const edgetype input){
    return input&0x1D; // 0x1D = 00011101
}

inline edgetype get_vu(const edgetype input){
    return input&0xE2;  // 0xE2 = 11100010
}

inline edgetype reverse(const edgetype et)
{
    return ((et & 2) >> 1) | ((et & 1) << 1) | ((et & 0x1C) << 3) | ((et & 0xE0) >> 3);
}

inline short get_color_u_v(const edgetype et)
{
    return ((et&0x1C) >> 2); // 0x1C = 00011100
}

inline short get_color_v_u(const edgetype et)
{
    return ((et&0xE0) >> 5); // 0xE0 = 11100000
}

inline void set_color_u_v(edgetype &et, const unsigned short color)
{
    et = (et&(0xE3))|(color << 2);       // 0xE3 = 11100011
}

inline void set_color_v_u(edgetype &et, const unsigned short color)
{
    et = (et&(0x1F))|(color << 5);       // 0x1F = 00011111 
}

inline edgetype colored_u_v(const edgetype et, const unsigned short color)
{
    return (et&(0xE3))|(color << 2);       // 0xE3 = 11100011
}

inline edgetype colored_v_u(const edgetype et, const unsigned short color)
{
    return (et&(0x1F))|(color << 5);       // 0x1F = 00011111 
}

inline void set_undir(edgetype &et, const bool use_u_v)
{
    et = UNDIR_U_V | (use_u_v?(et&0x1C << 3 | et&0x1C):(et&0xE0 >> 3 | et&0xE0));
}

inline edgetype colored_edge(const vertex u, const vertex v, const short color_uv, 
                         const short color_vu) {
    short ex_uv = color_uv > 0 ? 1 : 0;
    short ex_vu = color_vu > 0 ? 1 : 0;    
    edgetype ret = (u<v) ? (color_vu << 5) | (color_uv << 2) | (ex_vu << 1) | (ex_uv) 
                         : (color_uv << 5) | (color_vu << 2) | (ex_uv << 1) | (ex_vu);
    return ret;                      
}

inline void replace_neighbour(maingraph & maing, const vertex source, 
                              const vertex old, const vertex neu) {
    for (unsigned long i = 0; i != maing.num_neighbours[source]; ++i) {
        if (maing.neighbours[source][i] == old) {
           maing.neighbours[source][i] = neu;
           break;   
        }
    }
}

inline void add_neighbour(maingraph & maing, const vertex source, const vertex neu) {
    maing.neighbours[source][maing.num_neighbours[source]] = neu;
    maing.num_neighbours[source]++; 
}

inline void del_neighbour(maingraph & maing, const vertex source, const vertex old) {
    for (unsigned long i = 0; i != maing.num_neighbours[source]; ++i) {
        if (maing.neighbours[source][i] == old) {
           maing.neighbours[source][i] = maing.neighbours[source][maing.num_neighbours[source]-1];
           break;   
        }
    }
    maing.num_neighbours[source]--; 
}

#include "EdgeContainer.hpp"


// Maingraph functions:
bool read_graph(maingraph & maing, string filename, bool directed, 
                bool has_vertex_colors, bool has_edge_colors, short G_N,
                bool gen_dumpfile, string outputfile);
                 
void build_graph(maingraph & maing);

uint64 est_tree_size(const maingraph & g, long* v_extension, 
                     uint64 TREESMPLS, short G_N, randlib::rand &rand);
                     
double sampling(const maingraph & maing, long* v_extension, short G_N, 
                bool fullenumeration, const double* prob,
                const uint64 equiv100p, const int & perc_number,  
                hash_map < graphcode64, uint64 > & result_graphs, 
                uint64 & count_subgr, 
                randlib::rand &rand, bool gen_dumpfile, vector<subgraph>& subgraphdump);
                
double randomize_graph(maingraph & maing, short random_type,
                     int num_exchanges, int num_tries, 
                     const bool vertex_color_matters, const bool edge_color_matters,
                     EdgeContainer & EC, long & total_tries, long & total_success,
                     randlib::rand &rand);

// For the percentage event:
#define ID_PERCENT_REACHED 103

#endif
