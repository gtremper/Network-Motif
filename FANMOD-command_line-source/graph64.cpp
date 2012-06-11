#include "graph64.hpp"

static DEFAULTOPTIONS(options);
statsblk(stats);
setword nauty_workspace[160*MAXM];
set *nauty_gv;
static const uint64 ZRO_V_BIT = 0x00UL;
static const uint64 TWO_V_BIT = 0x01UL;
static const uint64 THR_V_BIT = 0x02UL;
static const uint64 FOU_V_BIT = 0x03UL;
static const uint64 V_BIT_MSK = 0x03UL;
static const uint64 ONE_E_BIT = 0x04UL;
static const uint64 TWO_E_BIT = 0x08UL;
static const uint64 THR_E_BIT = 0x0CUL;
static const uint64 E_BIT_MSK = 0x0CUL;

void init_graph(graph64 &g, short size, unsigned short num_vcolors, 
                unsigned short num_ecolors, bool directed) {
	for (int i = 0; i!=64 ; ++i) {
		g.matrix[i] = 0;
	}
	g.size = size;
	g.has_vertex_colors = (num_vcolors > 1);
	g.has_edge_colors = (num_ecolors > 1);
	g.num_vertex_colors = num_vcolors;
	g.num_edge_colors = num_ecolors;
	g.directed = directed;
	
	g.codestamp = ZRO_V_BIT;
	
	g.num_vertex_bits = 0;
	if (num_vcolors > 1)   {g.num_vertex_bits = 2; g.codestamp = TWO_V_BIT; }
	if (num_vcolors > 4)   {g.num_vertex_bits = 3; g.codestamp = THR_V_BIT; }	
	if (num_vcolors > 8)   {g.num_vertex_bits = 4; g.codestamp = FOU_V_BIT; }	
	   
	if (num_ecolors ==1)   {g.num_edge_bits = 1; g.codestamp |= ONE_E_BIT; }    
	if (num_ecolors > 1)   {g.num_edge_bits = 2; g.codestamp |= TWO_E_BIT; }
	if (num_ecolors > 3)   {g.num_edge_bits = 3; g.codestamp |= THR_E_BIT; }	
   
	
	g.g_N = (g.has_edge_colors) ? size*size : size;  // Graph is uncolored at edges
	g.g_M = (g.g_N + WORDSIZE - 1) / WORDSIZE;
	for (int i = 0; i != g.g_N; ++i) 
		EMPTYSET( ( GRAPHROW(g.nauty_g, i, g.g_M) ) , g.g_M);
	options.writeautoms = FALSE;
	options.getcanon = TRUE;

	options.defaultptn = (g.has_edge_colors || g.has_vertex_colors) ? FALSE : TRUE;
	
	if (directed) {
		options.digraph = TRUE;
		options.invarproc = adjacencies;
		options.mininvarlevel = 1;
		options.maxinvarlevel = 10;
	}

	nauty_check(WORDSIZE, g.g_M, g.g_N, NAUTYVERSIONID);

}

graphcode64 toHashCode(graph64 &g) {

	graphcode64 ret = 0ULL;
	
	if ((!g.has_vertex_colors) && (!g.has_edge_colors)) { //graph is not colored
	
		nauty(g.nauty_g, g.lab, g.ptn, NILSET,g. orbits, &options, &stats, 
			  nauty_workspace, 160*MAXM, g.g_M, g.g_N, g.nauty_canon);
			  
		for (int a = 0; a != g.size; ++a) {
			for (int b = 0; b != g.size; ++b) {
                if (a!=b) {
				   ret <<= 1;
   				   ret |= get_element(g,g.lab[a],g.lab[b]);
                }
			}
		}

		ret <<= 4;
        ret |= g.codestamp;
		
	} else { // graph is colored

		// convert for nauty
		int index = 0;
		uint32 sortarray[MAXN];

		// Build vertex partition
		if (g.has_vertex_colors) {
			while (index != g.size) {
				sortarray[index] = (0UL | index) | (g.matrix[index*8 +index] << 16);
				++index;
			}
		} else {
			while (index != g.size) {
				sortarray[index] = index;
				++index;
			}
		}
		register unsigned int edgecolor_ij;
		//Build edge partition
		if (g.has_edge_colors) {
			for (int i = 0; i != g.size; ++i) {
				for (int j = 0; j != g.size; ++j) {
					if (i != j)
					{
						edgecolor_ij = get_element(g,i,j);
						if (edgecolor_ij > 1) {
							//add vertex to naugraph
							 ADDELEMENT( ( GRAPHROW(g.nauty_g, i, g.g_M) ) , index);
							 ADDELEMENT( ( GRAPHROW(g.nauty_g, index, g.g_M) ) , j);
							sortarray[index] = (0UL | index) | (j<<8) | (i<<12) | (edgecolor_ij << 20);
							++index;
						}
					}
				}
			}
		}
		//Build fitting nauty partition
		sort(sortarray,sortarray+index);
		for (int i = 0; i != index; ++i)
		{
			g.lab[i] = sortarray[i] & 0x000000FF;
			if (i != index-1) {
				if ((sortarray[i] & 0x00FF0000) != (sortarray[i+1] & 0x00FF0000))
					g.ptn[i] = 0;
				else
					g.ptn[i] = 1;
			}
		}
		g.ptn[index-1] = 0; //ok since index always nonzero

		//perform nauty
		nauty(g.nauty_g, g.lab, g.ptn, NILSET, g.orbits, &options, &stats, 
			  nauty_workspace, 160*MAXM, g.g_M, index, g.nauty_canon);

		//1. g.lab contains permutation
		
		//2. Use this to get vertices/edges from g.matrix
		for (int a = 0; a != g.size; ++a) {
			for (int b = 0; b != g.size; ++b) {
				ret <<= (a==b) ? g.num_vertex_bits : g.num_edge_bits;
				ret |= get_element(g,g.lab[a],g.lab[b]);
			}
		}
        ret <<= 4;
        
        ret |= g.codestamp;
          		
		//3. Delete added vertices from the naugraph
		unsigned short source, target, inter;
		uint32 edge;
		for (int i = g.size; i != index; ++i)
		{
			edge = sortarray[i];
			source = (edge>>12) & 0x0F;
			target = (edge>> 8) & 0x0F;
			inter  = edge & 0xFF;
			DELELEMENT( ( GRAPHROW(g.nauty_g, source, g.g_M) ) , inter);
			DELELEMENT( ( GRAPHROW(g.nauty_g, inter, g.g_M) ) , target);
		}



	}
	return ret;
}


//
// TODO
//
void readHashCode(graph64 &g, graphcode64 gc) {
     unsigned short vbits = 0, ebits = 1 , vmsk = 0, emsk = 0;
     switch (gc & V_BIT_MSK) {
            case ZRO_V_BIT : vbits = 0; vmsk =  0; break;
            case TWO_V_BIT : vbits = 2; vmsk =  3; break;
            case THR_V_BIT : vbits = 3; vmsk =  7; break;                        
            case FOU_V_BIT : vbits = 4; vmsk = 15; break;            
     }
     switch (gc & E_BIT_MSK) {
            case ONE_E_BIT : ebits = 1; emsk = 1; break;
            case TWO_E_BIT : ebits = 2; emsk = 3; break;
            case THR_E_BIT : ebits = 3; emsk = 7; break;                        
     }     
     
     gc >>= 4; // remove codestamp
     
     for (int a = g.size-1; a >= 0; --a) {
         for (int b = g.size-1; b >= 0; --b) {
             if (a==b) {
                set_element(g, a, b, gc & vmsk); 
                gc >>= vbits;         
             } else {
                set_element(g, a, b, gc & emsk);
                gc >>= ebits;                    
             }
         }
     }     
}


graphcode64 getGraphID(graph64 &g, graphcode64 gc) {
     unsigned short vbits = 0, ebits = 1 , vmsk = 0, emsk = 0;
     switch (gc & V_BIT_MSK) {
            case ZRO_V_BIT : vbits = 0; vmsk =  0; break;
            case TWO_V_BIT : vbits = 2; vmsk =  3; break;
            case THR_V_BIT : vbits = 3; vmsk =  7; break;                        
            case FOU_V_BIT : vbits = 4; vmsk = 15; break;            
     }
     switch (gc & E_BIT_MSK) {
            case ONE_E_BIT : ebits = 1; emsk = 1; break;
            case TWO_E_BIT : ebits = 2; emsk = 3; break;
            case THR_E_BIT : ebits = 3; emsk = 7; break;                        
     }     
     
     gc >>= 4; // remove codestamp
     
    graph nau_c[MAXN * MAXM];
	graph nau_g[MAXN * MAXM];
	short gn = g.size;
	short gm = (gn + WORDSIZE - 1) / WORDSIZE;
	for (int i = 0; i != g.g_N; ++i) 
		EMPTYSET( ( GRAPHROW(nau_g, i, gm) ) , gm);	

     for (int a = g.size-1; a >= 0; --a) {
         for (int b = g.size-1; b >= 0; --b) {
             if (a==b) {
                gc >>= vbits;         
             } else {
				if ((gc & emsk) > 0)
					ADDELEMENT( ( GRAPHROW(nau_g, a, gm) ) , b);
                gc >>= ebits;
             }
         }
     }  
     
     options.defaultptn = TRUE;
     
 	 nauty(nau_g, g.lab, g.ptn, NILSET, g.orbits, &options, &stats, 
			  nauty_workspace, 160*MAXM, gm, gn, nau_c);
     
     options.defaultptn = (g.has_edge_colors || g.has_vertex_colors) ? FALSE : TRUE;
     
     graphcode64 ret = 0;     
     
	for (int a = 0; a != g.size; ++a) {
		for (int b = 0; b != g.size; ++b) {
		   ret <<= 1;
		   if ( ISELEMENT( ( GRAPHROW(nau_c, a, gm) ) , b) ) 
		       ret |= 1;
		}
	}     
     
     
	 return ret;
}
