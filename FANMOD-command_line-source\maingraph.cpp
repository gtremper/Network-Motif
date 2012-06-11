#include "maingraph.hpp"
using namespace std;

/*********************
 *  UTILITIES        *
 *********************/

// Tries to find the specified value in the array and replaces it.
// Important: assumes that the value exists!
void find_and_replace(edge* array, edge find, edge replace){
    long i = 0;
    while (array[i]!=find)
		++i;
    array[i]= replace;
}


/*********************
 *  MAIN-Definitions *
 *********************/

// read a graph from "filename" and return it.
bool read_graph(maingraph & maing, string filename, bool directed, 
                bool has_vertex_colors, bool has_edge_colors, short G_N,
                bool gen_dumpfile, string outputfile)
{
	vector <long> num_neighbours; // stores the number of neighbours of a vertex
    vector <long> num_undir; // stores the number of undirected edges from a vertex
    vector <short> vertex_colors;
    ifstream in(filename.c_str());
    const unsigned short NO_COLOR = 255;
    std::fstream logfile((outputfile+".log").c_str(),ios::app | ios::out);
    logfile.close();
	if (in.good()) {
	// Initialize variables
	maing.n = 0;
	maing.num_nodes = 0;
	maing.num_lonely_nodes = 0;
    maing.num_vertex_colors = 0;
    maing.num_edge_colors = 1;
    maing.directed = directed;
    unsigned short color1 = NO_COLOR, color2 = NO_COLOR, color3 = NO_COLOR, 
                   edge_color = NO_COLOR, color_u = NO_COLOR, color_v = NO_COLOR;
    unsigned short edge_color_limit = 100, vertex_color_limit = 100;
    unsigned short curr_colors = 0; // number of colors in the current line
    if (has_vertex_colors){
        switch (G_N){
            case 3: 
            case 4: vertex_color_limit = 15; break;
            case 5: vertex_color_limit = 15; break;
            case 6: vertex_color_limit = 15; break;
            case 7: vertex_color_limit = 3; break;
            case 8: vertex_color_limit = 0; break;
        }
    }
    if (has_edge_colors){
        switch (G_N){
            case 3: 
            case 4: edge_color_limit = 7; vertex_color_limit = 15; break;
            case 5: edge_color_limit = 3; vertex_color_limit = 7; break;
            case 6: edge_color_limit = 3; vertex_color_limit = 0; break;
            case 7: edge_color_limit = 1; vertex_color_limit = 0; break;
            case 8: edge_color_limit = 1; vertex_color_limit = 0; break;
        }
    }
    //wxLogMessage(wxString()<<edge_color_limit << " "<<vertex_color_limit);
    vertex u;
	vertex v;
	unsigned long linenumber = 0;
	edge e;
	// Parse input File
	string line;
	stringstream ss;
	while( std::getline( in, line ) ) {
	    line += " ";
     	ss.str(line);
	    ss >> u >> v;
	    if (!ss.fail()) {
			curr_colors = 0;
            if (ss >> color1) ++curr_colors;
            if (ss >> color2) ++curr_colors;
            if (ss >> color3) ++curr_colors;
			ss.clear(); //Reset Stream

			vertex maxuv = (u>v) ? u : v;
            if (maxuv > maing.n) {
				maing.n = maxuv;
                if (gen_dumpfile && (maxuv > 65535)){
                    // print to log file
                    logfile.open((outputfile+".log").c_str(),ios::app);
					logfile<<"ERROR: In Line "<< linenumber << ", vertex number " 
                        << maxuv << " was specified.\nLargest vertex number allowed with " 
                        << " an activated subgraph dump is 65535. Exceeded vertex limit"<<endl;
					logfile.close();
                    return false;
                }
				if (num_neighbours.size() <= maing.n) {
					num_neighbours.resize(maing.n+1);
                    num_undir.resize(maing.n+1);
                    if (has_vertex_colors) {
                        unsigned long oldsize = vertex_colors.size();
                        vertex_colors.resize(maing.n+1);
                        for (unsigned long i = oldsize; i != vertex_colors.size(); ++i)
                            vertex_colors[i] = NO_COLOR;
                    }
				}
			}

			if (u != v) {
				// which color belongs to which graph-part
                switch (curr_colors){
                    case 0: edge_color = NO_COLOR; color_u = NO_COLOR; 
                            color_v = NO_COLOR; break;
                    case 1: edge_color = color1; color_u = NO_COLOR; 
                            color_v = NO_COLOR; break;
                    case 2: edge_color = NO_COLOR; color_u = color1; 
                            color_v = color2; break;
                    case 3: edge_color = color3; color_u = color1; 
                            color_v = color2; break;
                }
                // set number of colors
                if (has_edge_colors && edge_color != NO_COLOR && edge_color > maing.num_edge_colors)
                    maing.num_edge_colors = edge_color;
                if (has_vertex_colors && color_u != NO_COLOR && color_u > maing.num_vertex_colors)
                    maing.num_vertex_colors = color_u;
                if (has_vertex_colors && color_v != NO_COLOR && color_v > maing.num_vertex_colors)
                    maing.num_vertex_colors = color_v;
                // check color limits
                if (has_edge_colors && edge_color == 0){
					// print to log file
                    logfile.open((outputfile+".log").c_str(),ios::app);
                    logfile<<"ERROR: In Line "<< linenumber << ", edge-color 0 was given.\n"
                           <<"This color is not allowed - edges must always have a color >=1." 
                           <<"Exceeded color limits"<<endl;
                    logfile.close();
					return false;
                }                
                if (has_edge_colors && maing.num_edge_colors > edge_color_limit){
                    // print to log file
					logfile.open((outputfile+".log").c_str(),ios::app);
                    logfile<<"ERROR: In Line "<< linenumber <<", an edge-color >" 
                        << edge_color_limit << " was given.\nLargest edge-color allowed with " 
                        << "your settings is " << edge_color_limit << "." 
                        << "Exceeded color limits" << endl;
                    logfile.close();
					return false;
                }
                if (has_vertex_colors && maing.num_vertex_colors > vertex_color_limit){
                    // print to log file
                    logfile.open((outputfile+".log").c_str(),ios::app);
                    logfile<<"ERROR: In Line "<< linenumber <<", a vertex-color >" 
                        << vertex_color_limit << " was given.\nLargest vertex-color allowed with " 
                        << "your settings is " << edge_color_limit << ". Exceeded color limits" << endl;
                    logfile.close();
					return false;
                }
                // set vertex colors
                if (!has_edge_colors)
                   edge_color = 1;
                if (has_vertex_colors && color_u != NO_COLOR) vertex_colors[u]= color_u;
                if (has_vertex_colors && color_v != NO_COLOR) vertex_colors[v]= color_v;
                if (u < v) {
                    e = new_edge(u, v);
                    if (maing.edges.find(e) == maing.edges.end()) {
						num_neighbours[u]++;
						num_neighbours[v]++;
						maing.edges[e] = DIR_U_T_V;
                        if(edge_color == NO_COLOR)
                            set_color_u_v(maing.edges[e],1); 
                        else 
                            set_color_u_v(maing.edges[e],edge_color);
					} else {
						maing.edges[e] |= DIR_U_T_V;
                        if(edge_color == NO_COLOR) {
                            if (get_color_u_v == 0) 
                               set_color_u_v(maing.edges[e],1); 
                        } else { 
                            set_color_u_v(maing.edges[e],edge_color);
                        }
					}
					if (!directed) {
						set_undir(maing.edges[e], (edge_color!=NO_COLOR));
						// use the color of u_v (pass true to set_undir), if edge_color is not NO_COLOR
					}
				} else { // v < u
					e = new_edge(v, u);
					if (maing.edges.find(e) == maing.edges.end()) {
						num_neighbours[u]++;
						num_neighbours[v]++;
						maing.edges[e] = DIR_V_T_U;
                        if(edge_color == NO_COLOR)
                            set_color_v_u(maing.edges[e],1); 
                        else 
                            set_color_v_u(maing.edges[e],edge_color);

					} else {
						maing.edges[e] |= DIR_V_T_U;
                        if(edge_color == NO_COLOR) {
                            if (get_color_v_u == 0) 
                               set_color_v_u(maing.edges[e],1); 
                        } else {
                            set_color_v_u(maing.edges[e],edge_color);
                        }
					}
					if (!directed) {
						set_undir(maing.edges[e], (edge_color==NO_COLOR));
                        // use the color of v_u (pass false to set_undir) if edge_color is not NO_COLOR
					}
				}
                if (rm_colors(maing.edges[e]) == UNDIR_U_V) {
                    ++num_undir[u];
                    ++num_undir[v];
                }
			}
	    } else { // if ss.fail
			// print to log file
			logfile.open((outputfile+".log").c_str(),ios::app);
            logfile<<"ERROR: Wrong format for input file \'"<< filename.c_str() << "\' in line " << linenumber << "\nExpected line format is <int> <int> [<int>] [<int>] [<int>]. Error reading file"<<endl;
            logfile.close();
			return false;
		}
		++linenumber;
	}
    } else { // if !in.good
		// print to log file
		logfile.open((outputfile+".log").c_str(),ios::app);
		logfile<<"ERROR: Unable to open file \'"<< filename.c_str()<< "\' for input. Error reading file"<<endl;
		logfile.close();
        return false;
    }
    ++maing.n;
	in.close();
    // Allocates the arrays of the maingraph    
    maing.num_neighbours = new unsigned long[maing.n+1];
    maing.neighbours = new vertex *[maing.n+1];
    if (has_vertex_colors)
       maing.vertex_colors = new short[maing.n+1];	
    //maing.num_larger_neighbours = new unsigned long[maing.n+1];
    maing.v_util = new unsigned long[maing.n+1];
    for (vertex i = 0; i != maing.n; ++i) {
      if (has_vertex_colors) maing.vertex_colors[i] = (vertex_colors[i] != NO_COLOR) ? vertex_colors[i] : 0; // copy vertex colors
      if (maing.directed)
          maing.neighbours[i] = new vertex[num_neighbours[i]+num_undir[i]];
      else
          maing.neighbours[i] = new vertex[num_neighbours[i]];
    }
    ++maing.num_vertex_colors;
    //++maing.num_edge_colors;    <<- STARTS AT 1 AND NOT 0                     
    // the vectors are no longer needed
    vertex_colors.clear();
    num_neighbours.clear();
    num_undir.clear();
    return true;
}

void build_graph(maingraph & maing)
{
    for (vertex i = 0; i != maing.n; ++i) {
      maing.v_util[i] = NILLVERTEX;          // reset util fields;
	  //maing.num_larger_neighbours[i] = 0;
      maing.num_neighbours[i] = 0;           // num_neighbours might have changed during randomize
    }
    maing.m = maing.edges.size();       // Number of edges might have changed during "no regard"-randomization
	//maing.maxnumlargerneighbors = 0;
	//maing.numlargerneighbors = 0;	
	maing.maxnumneighbours = 0;
    maing.num_undir_edges = 0;
    maing.num_dir_edges = 0;

    for (hash_map < edge, edgetype >::const_iterator iter = maing.edges.begin();
	 iter != maing.edges.end(); ++iter) {
		vertex u = edge_get_u(iter->first);
		vertex v = edge_get_v(iter->first);
		edgetype etype = rm_colors(iter->second);
        if (etype == UNDIR_U_V) {
			++maing.num_undir_edges;
        } else {
			++maing.num_dir_edges;
		}
		(maing.neighbours[u])[maing.num_neighbours[u]] = v;
		++maing.num_neighbours[u];
		if (maing.num_neighbours[u] > maing.maxnumneighbours)
		   maing.maxnumneighbours = maing.num_neighbours[u];
		(maing.neighbours[v])[maing.num_neighbours[v]] = u;
		++maing.num_neighbours[v];
		if (maing.num_neighbours[v] > maing.maxnumneighbours)
		   maing.maxnumneighbours = maing.num_neighbours[v];		
		//if (u>v) {
		//	++maing.num_larger_neighbours[v];
		//	if (maing.num_larger_neighbours[v]> maing.maxnumlargerneighbors)
		//	    ++maing.maxnumlargerneighbors;
		//} else { // v>u
		//	++maing.num_larger_neighbours[u];
		//	if (maing.num_larger_neighbours[u]> maing.maxnumlargerneighbors)
		//	    ++maing.maxnumlargerneighbors;
		//}
    }

	//Sort the neighbour array for every vertex
	for (unsigned int i = 0; i != maing.n; ++i)
		std::sort(maing.neighbours[i], maing.neighbours[i] + maing.num_neighbours[i]);

	return;
}

// estimates the size of the subgraph tree using TREESMPLS samples.
uint64 est_tree_size(const maingraph & g, long* v_extension, uint64 TREESMPLS, 
                     short G_N, randlib::rand &rand)
{
	register vertex* v_subgraph = new vertex[G_N+1];
	register long counter;
	register vertex SOURCE_VERTEX; 
	register vertex FILL_VERTEX;
	register long* min_scope = new long[G_N+1];
	register long* max_scope = new long[G_N+1];
	register long* scope_place = new long[G_N+1];
	register short depth;
    uint64 APPROX_TREE_SIZE = 0;
	for (uint64 th = 0; th!=TREESMPLS; ++th)  {
		vertex v = ++rand % g.n;
		uint64 approx = 1;
		bool flag = false;
		if (g.num_neighbours[v] > 0) {
			depth = 1;
			min_scope[depth] = 0;
			max_scope[depth] = 1;
			scope_place[depth] = 0;
			v_extension[0] = v;
			v_subgraph[depth-1] = NILLVERTEX;
			g.v_util[v] = v;
			while (depth != 0)  {
				if (flag) { //scope_place[depth] == max_scope[depth] || flag)	{
					if (min_scope[depth] != max_scope[depth]) {
						counter = max_scope[depth] - 1;
						while (counter >= 0 && g.v_util[v_extension[counter]] == v_subgraph[depth-1]) {				
							g.v_util[v_extension[counter]] = NILLVERTEX;
							--counter;
						}
					}
					--depth;
				} else {
					if (depth == G_N) { 
						flag = true; 
						while (scope_place[depth] != max_scope[depth]) {
							v_subgraph[depth] = v_extension[scope_place[depth]];
							++scope_place[depth];
						}
					} else {  
						if (scope_place[depth] == max_scope[depth]) {
							approx = 0;
							flag = true;
						} else {
							scope_place[depth] += (++rand) % (max_scope[depth]-scope_place[depth]); 
							SOURCE_VERTEX = v_extension[scope_place[depth]];
							++scope_place[depth];
							v_subgraph[depth] = SOURCE_VERTEX;
							min_scope[depth+1] = max_scope[depth];
							scope_place[depth+1] = scope_place[depth];
							max_scope[depth+1] = min_scope[depth+1];
							++depth;
							counter = g.num_neighbours[SOURCE_VERTEX] - 1;
							if (counter > 0)   {
								FILL_VERTEX = g.neighbours[SOURCE_VERTEX][counter];
							}
							while (counter > -1 && g.neighbours[SOURCE_VERTEX][counter] > v)	{
								FILL_VERTEX = g.neighbours[SOURCE_VERTEX][counter];
								if (g.v_util[FILL_VERTEX] == NILLVERTEX)	{
									v_extension[max_scope[depth]] = FILL_VERTEX;
									g.v_util[FILL_VERTEX] = SOURCE_VERTEX;
									++max_scope[depth];
								}
								--counter;
							}
							approx *= (max_scope[depth] - scope_place[depth]);
						}
					}
				}
			}
			g.v_util[v] = NILLVERTEX;
			APPROX_TREE_SIZE += approx;
		}

	}
    delete[] v_subgraph;
    delete[] min_scope;
    delete[] max_scope;
    delete[] scope_place;

    return APPROX_TREE_SIZE;
}

// Samples / Enumerates the graph returns results in result_graphs, and the time directly.
double sampling(const maingraph & maing, long* v_extension, short G_N, 
                bool fullenumeration, const double* prob, const uint64 equiv100p,
                const int & perc_number, hash_map < graphcode64, uint64 > & result_graphs,
                uint64 & count_subgr, 
                randlib::rand &rand, bool gen_dumpfile, vector<subgraph>& subgraphdump)
{

	// Init workgraph
	graph64 g;
    init_graph(g,G_N,maing.num_vertex_colors,maing.num_edge_colors,maing.directed); 

    // Init for Statusbar
    register int perc_index = 0;
	
    // Init for main loop
	count_subgr = 0;
	register vertex v_subgraph[G_N_MAX+1];
	register long counter;
	register vertex SOURCE_VERTEX; 
	register vertex FILL_VERTEX;
	register long min_scope[G_N_MAX+1];
	register long max_scope[G_N_MAX+1];
	register long scope_place[G_N_MAX+1];
	register short depth;
    const bool vcolors = (maing.num_vertex_colors > 1);
    const bool ecolors = (maing.num_edge_colors > 1);
    
  	register vertex* start_vertices = new vertex[maing.n+2];//##%!
	register short scope_place_point[G_N_MAX+1];//##%!
	register unsigned long* scope_place_loc[G_N_MAX+1];

    //for dumpfile
    //vector<subgraph> subgraphdump; //stores all found subgraphs in three
                                 //consecutive uint64

	// Beginning actual sampling
    for (int i = 0; i != G_N+1; ++i) {
		scope_place_loc[i] = new unsigned long[maing.maxnumneighbours*G_N];
	}

	//scope_place_loc[8][maxnumlargerneighbors*G_N];//##%!
	//gen_selection(0,tst,0.174,arr,rand);
	if (fullenumeration) {
		for (vertex v = 0; v != maing.n; ++v)
			start_vertices[v] = v;
		start_vertices[maing.n] = maing.n;
	} else {
		gen_selection(0,maing.n,prob[0],start_vertices,rand);
	}
	
	clock_t start_time(clock());
	for (int th = 0; th!=1; ++th) // this line for time measurements only 
	{

	int idx = 0;
	while (start_vertices[idx]!=maing.n)
	{	
		vertex v = start_vertices[idx];
		++idx;
		//if (fullenumeration || ((++rand)%MAXPROB) <= prob[0])///////////////////////////////////////////////
		{
		
		depth = 1;
		min_scope[depth] = 0;
		max_scope[depth] = 1;
		scope_place[depth] = 0;
		v_extension[0] = v;
		v_subgraph[depth-1] = NILLVERTEX;
		maing.v_util[v] = v;
		scope_place_point[depth] = 0;
		scope_place_loc[1][0] = 0;
		scope_place_loc[1][1] = 1;

		while (depth != 0)  {
			if (scope_place[depth] == max_scope[depth])	{ // go lower
				if (min_scope[depth] != max_scope[depth]) {
					counter = max_scope[depth] - 1;
					while (counter >= 0 && maing.v_util[v_extension[counter]] == v_subgraph[depth-1]) {
						maing.v_util[v_extension[counter]] = NILLVERTEX;
						--counter;
					}
				}
				--depth;
			} else {
				if (depth == G_N) { 
					while (scope_place[depth] != max_scope[depth])
					{
						{
						++count_subgr;
						v_subgraph[depth] = v_extension[scope_place[depth]];
						if (vcolors)
						   color_vertex(g,depth-1,maing.vertex_colors[v_subgraph[depth]]);

						// edit nauty-graph 
						for (int i = 0; i != depth-1; ++i) {
							vertex uc = v_subgraph[i+1];
							vertex vc = v_subgraph[depth];
							edge e_check = edge_code(uc, vc);
							delete_edge(g,depth-1,i);
							delete_edge(g,i,depth-1);
							if (maing.edges.find(e_check) != maing.edges.end()) {
								edgetype ec = maing.edges.find(e_check)->second;
								if (ecolors)
								{ // colored egdes
									if (uc < vc) {
										switch (rm_colors(ec)) {
											case DIR_U_T_V : set_edge(g,i,depth-1,get_color_u_v(ec));	break;
											case DIR_V_T_U : set_edge(g,depth-1,i,get_color_v_u(ec));	break;
											case UNDIR_U_V : set_edge(g,i,depth-1,get_color_u_v(ec));
															 set_edge(g,depth-1,i,get_color_v_u(ec));	break;
										}
									} else {
										switch (rm_colors(ec)) {
											case DIR_U_T_V : set_edge(g,depth-1,i,get_color_u_v(ec));	break;
											case DIR_V_T_U : set_edge(g,i,depth-1,get_color_v_u(ec));	break;
											case UNDIR_U_V : set_edge(g,depth-1,i,get_color_u_v(ec));
															 set_edge(g,i,depth-1,get_color_v_u(ec));	break;
										}
									}
                                } else { //no colored edges
									if (uc < vc) {
										switch (rm_colors(ec)) {
											case DIR_U_T_V : set_edge(g,i,depth-1);	break;
											case DIR_V_T_U : set_edge(g,depth-1,i);	break;
											case UNDIR_U_V : set_edge(g,i,depth-1);
															 set_edge(g,depth-1,i);	break;
										}
									} else {
										switch (rm_colors(ec)) {
											case DIR_U_T_V : set_edge(g,depth-1,i);	break;
											case DIR_V_T_U : set_edge(g,i,depth-1);	break;
											case UNDIR_U_V : set_edge(g,i,depth-1);
															 set_edge(g,depth-1,i);	break;
										}
									}                                
                                }
								
							}
							
						}
						// end edit nauty-graph

						// Canonical labelling with nauty
						graphcode64 hashg = toHashCode(g);
						result_graphs[ hashg ]++;
						//Dump subgraph if requested
                        if (gen_dumpfile)
						   subgraphdump.push_back(get_subgraph(v_subgraph, G_N, hashg));
						// end nauty

						if (fullenumeration)
						{
							++scope_place[depth];
						} else {
							++scope_place_point[depth];
							scope_place[depth] = scope_place_loc[depth][scope_place_point[depth]];						
						}
						}

					}


//----------		// Progress indicator event
                    if (perc_number > 0) // send events only when necessary
                    while (count_subgr > (equiv100p * perc_index * perc_number / 100) && perc_index < 100 / perc_number) {
						++perc_index;
					}

				} else {  //go deeper
					SOURCE_VERTEX = v_extension[scope_place[depth]];
					if (fullenumeration) {
						++scope_place[depth];
					} else {
						++scope_place_point[depth];
						scope_place[depth] = scope_place_loc[depth][scope_place_point[depth]];
					}

					{//!!

					v_subgraph[depth] = SOURCE_VERTEX;
					if (vcolors)
					   color_vertex(g,depth-1,maing.vertex_colors[SOURCE_VERTEX]);
					
					// edit nauty-graph 
					for (int i = 0; i != depth-1; ++i) {
						vertex uc = v_subgraph[i+1];
						vertex vc = SOURCE_VERTEX;
						edge e_check = edge_code(uc, vc);
						delete_edge(g,depth-1,i);
						delete_edge(g,i,depth-1);
						if (maing.edges.find(e_check) != maing.edges.end()) {
							edgetype ec = maing.edges.find(e_check)->second;
							if (ecolors)
							{ // colored egdes
								if (uc < vc) {
									switch (rm_colors(ec)) {
										case DIR_U_T_V : set_edge(g,i,depth-1,get_color_u_v(ec));	break;
										case DIR_V_T_U : set_edge(g,depth-1,i,get_color_v_u(ec));	break;
										case UNDIR_U_V : set_edge(g,i,depth-1,get_color_u_v(ec));
														 set_edge(g,depth-1,i,get_color_v_u(ec));	break;
									}
								} else {
									switch (rm_colors(ec)) {
										case DIR_U_T_V : set_edge(g,depth-1,i,get_color_u_v(ec));	break;
										case DIR_V_T_U : set_edge(g,i,depth-1,get_color_v_u(ec));	break;
										case UNDIR_U_V : set_edge(g,depth-1,i,get_color_u_v(ec));
														 set_edge(g,i,depth-1,get_color_v_u(ec));	break;
									}
								}
							} else { //no colored edges
								if (uc < vc) {
									switch (rm_colors(ec)) {
										case DIR_U_T_V : set_edge(g,i,depth-1);	break;
										case DIR_V_T_U : set_edge(g,depth-1,i);	break;
										case UNDIR_U_V : set_edge(g,i,depth-1);
														 set_edge(g,depth-1,i);	break;
									}
								} else {
									switch (rm_colors(ec)) {
										case DIR_U_T_V : set_edge(g,depth-1,i);	break;
										case DIR_V_T_U : set_edge(g,i,depth-1);	break;
										case UNDIR_U_V : set_edge(g,i,depth-1);
														 set_edge(g,depth-1,i);	break;
									}
								}                                
							}
								
						}
						
					}
					// end edit nauty-graph

					min_scope[depth+1] = max_scope[depth];
					scope_place[depth+1] = scope_place[depth];
					max_scope[depth+1] = min_scope[depth+1];
					++depth;
					counter = maing.num_neighbours[SOURCE_VERTEX] - 1;
					if (counter > 0)   {
						FILL_VERTEX = maing.neighbours[SOURCE_VERTEX][counter];
					}
					while (counter > -1 && maing.neighbours[SOURCE_VERTEX][counter] > v)	{
						FILL_VERTEX = maing.neighbours[SOURCE_VERTEX][counter];
						if (maing.v_util[FILL_VERTEX] == NILLVERTEX)	{
							v_extension[max_scope[depth]] = FILL_VERTEX;
							maing.v_util[FILL_VERTEX] = SOURCE_VERTEX;
							++max_scope[depth];
						}
						--counter;
					}
					}//!!
					
					//GENERATE SELECTION AND RESET SCOPE PLACE POINT
					if (! fullenumeration)
					{
						scope_place_point[depth] = 0;
						gen_selection(scope_place[depth], max_scope[depth],
						              prob[depth-1],scope_place_loc[depth],rand);
						scope_place[depth] = scope_place_loc[depth][0];
					}
				}
			}
		}
		maing.v_util[v] = NILLVERTEX;
		}
	}

	}

    //Clean up
    //delete[] v_subgraph;
    //delete[] min_scope;
    //delete[] max_scope;
    //delete[] scope_place;
    delete[] start_vertices;
    //delete[] scope_place_point;
	for (int i = 0; i != G_N+1; ++i) {
		delete[] scope_place_loc[i];
	}
    //delete[] scope_place_loc;

   	// Done with sampling / enumeration: stop the clock, return the time
    return double (clock() - start_time) / CLOCKS_PER_SEC;
}


/********************* RANDOMIZATION **********************************/


double randomize_graph(maingraph & maing, short random_type,
                     int num_exchanges, int num_tries, 
                     const bool vertex_color_matters, const bool edge_color_matters,
                     EdgeContainer & EC, long & total_tries, long & total_success,
                     randlib::rand &rand) {
    
    clock_t start_time(clock()); // Start the time measurement                     
    long rule4 = 0, rule5 = 0;
    
    switch (random_type) {
    case NO_REGARD :  

         for (int exchanges = 0; exchanges != num_exchanges; ++exchanges) {
    	     EC.resetIteration();
    	     while (EC.hasNextBag()) {
                   EC.goNextBag();
                   bagID bag_id = EC.getCurrentBag();
                   bag_index bag_size = EC.getBagSize(bag_id);
    		       while (EC.hasNextElement()) {
    			         EC.goNextElement();
    			         bag_index current_index = EC.getCurrentIndex();
    			         //Begin switch attempt
    			         for (int tries = 0; tries != num_tries; tries++) {
                             ++total_tries;
                             if (bag_size == 1)
                             { total_tries+= (num_tries-1);  break;}

							 //choose random exchange partner
                             bag_index random_index = (++rand) % bag_size;
                             while (random_index == current_index) {
                                   random_index = (++rand) % bag_size;
                             }
                             
                             //Check if switch can be performed

							 //Uses invariant that only unidirectional edges 
                             //may be present in the EdgeContainer
                             
                             //May be noncanonical edgecode!! (Direction
                             //is implicitly stored here)
                             const edge e = EC.getCurrentElement();
                             const edge f = EC.getIndexedElement(random_index);
							 
                             const vertex eu = edge_get_u(e);
                             const vertex ev = edge_get_v(e);
                             const vertex fu = edge_get_u(f);
                             const vertex fv = edge_get_v(f);
                             
                             if (eu != fu && eu != fv && ev != fu && ev != fv) { //check vertices
                             
                             const edge e_code = edge_code(eu,ev);
                             const edge f_code = edge_code(fu,fv);
                             const edge eufv_code = edge_code(eu,fv);
                             const edge fuev_code = edge_code(fu,ev);

                             const edgetype et = maing.edges[e_code],
                                            ft = maing.edges[f_code];
                             
                             const short color_euev = (eu < ev) ? get_color_u_v(et) : get_color_v_u(et);
                             const short color_fufv = (fu < fv) ? get_color_u_v(ft) : get_color_v_u(ft);
                             
                             const bool evfu_exists = maing.edges.find(fuev_code) != maing.edges.end();
                             const bool fveu_exists = maing.edges.find(eufv_code) != maing.edges.end();
                             
                             bool can_switch = true;
                             
                            // if (rm_colors(et) == UNDIR_U_V || rm_colors(ft) == UNDIR_U_V)
                            //    can_switch = false;
                             
                             if (fveu_exists) {

                                const edgetype n_et = maing.edges[eufv_code];
                                if (eu < fv && ((n_et & DIR_U_T_V) == DIR_U_T_V))
                                   can_switch = false;
                                else if (eu > fv && ((n_et & DIR_V_T_U) == DIR_V_T_U))
                                   can_switch = false;                                
                             }
                             if (evfu_exists) { 

                                const edgetype n_ft = maing.edges[fuev_code];
                                if (fu < ev && ((n_ft & DIR_U_T_V) == DIR_U_T_V))
                                   can_switch = false;
                                else if (fu > ev && ((n_ft & DIR_V_T_U) == DIR_V_T_U))
                                   can_switch = false;                                
                             }
                             
                             if (can_switch) {
                                //Perform
                                //perform on EC (bagID remains here)
                                
                                EC.setIndexedElement(current_index,new_edge(eu,fv));
                                EC.setIndexedElement(random_index,new_edge(fu,ev));                                
                                ///*
                                //4 Things to consider when manipulating HashMap: 
                                // 1. del euev: may be bidirected
                                const edgetype net = (eu < ev) ? get_vu(et) : get_uv(et);
                                if (net == 0) {
								    maing.edges.erase(e_code);
                                } else {
                                    maing.edges[e_code] = net;
                                }
                                
                                // 2. del fufv: may be bidirected

                                const edgetype nft = (fu < fv) ? get_vu(ft) : get_uv(ft);
                                if (nft == 0) {
								    maing.edges.erase(f_code);
                                } else {
                                    maing.edges[f_code] = nft;
                                }
                                
                                // 3. set eu->fv: eu<-fv may already exist
                                edgetype type_eufv = (eu < fv) ? (colored_u_v(DIR_U_T_V,color_euev))
                                                               : (colored_v_u(DIR_V_T_U,color_euev));
                                if (fveu_exists) {
                                   maing.edges[eufv_code] |= type_eufv;
                                } else {
                                   maing.edges[eufv_code] = type_eufv;  
                                } 
                                
                                // 4. set fu->ev: fu<-ev may already exist
                                edgetype type_fuev = (fu < ev) ? (colored_u_v(DIR_U_T_V,color_fufv))
                                                               : (colored_v_u(DIR_V_T_U,color_fufv));
                                if (evfu_exists) {
                                   maing.edges[fuev_code] |= type_fuev;
                                } else {
                                   maing.edges[fuev_code] = type_fuev;  
                                } 
                                                             
                                //*/
                                
                                ++total_success;
                                break;//Break try-loop                                                                                                          
                             }//end exchange
                             }//end check vertices
                         }//end exchange attempt
                         
                   }
             }
         }
         
         //Update Maingraph
         build_graph(maing);
         
    break; // break case NO_REGARD
    
    case GLOBAL_CONST : 


         for (int exchanges = 0; exchanges != num_exchanges; ++exchanges) {
    	     EC.resetIteration();
    	     while (EC.hasNextBag()) {
                   EC.goNextBag();
                   bagID bag_id = EC.getCurrentBag();
                   bag_index bag_size = EC.getBagSize(bag_id);
    		       while (EC.hasNextElement()) {
    			         EC.goNextElement();
    			         bag_index current_index = EC.getCurrentIndex();
    			         //Begin switch attempt
    			         for (int tries = 0; tries != num_tries; tries++) {

                             ++total_tries;
                             const edge e = EC.getCurrentElement();
                             vertex eu = edge_get_u(e);
                             vertex ev = edge_get_v(e);
                             const edgetype et = maing.edges[e];   
							 short color_euev = get_color_u_v(et);
							 short color_eveu = get_color_v_u(et);                                                       
                             
//-------------------------------------------------
if (rm_colors(et) == UNDIR_U_V) {
//---------------------Start new global BIDIR rules

    //Try to apply rule 5
    vertex cand_w = NILLVERTEX, cand_x = NILLVERTEX, cand_y = NILLVERTEX;
    short cand_color_wx = 0, cand_color_xy = 0;
    double num_rule5_candiates = 0.0;

    //Determine the direction in which we are searching
    if (rand.trueWithProb(0.5)) {
       const vertex tmp = eu;
       eu = ev;
       ev = tmp;
       const short tmpcolor = color_euev;
       color_euev = color_eveu;
       color_eveu = tmpcolor;                           
    }
    
    //find appropriate edge in bag
    //retrieve bagID
    short colr_u = vertex_color_matters ? maing.vertex_colors[eu] : 1;
    short colr_v = vertex_color_matters ? maing.vertex_colors[ev] : 1;
    short colr_euev = edge_color_matters ? color_euev : 1;
//    short colr_eveu = edge_color_matters ? color_eveu : 1;    
    const bagID cand_ebag = getBagID(colr_euev,0,colr_u,colr_v);
    const bag_index ebag_size = EC.getBagSize(cand_ebag);    
        
    if (ebag_size > 0) {
    const edge edge_wx = EC.getElement( (++rand)%ebag_size ,cand_ebag);
    const edgetype etype_wx = maing.edges[edge_wx];
    const vertex vert_w = (rm_colors(etype_wx) == DIR_U_T_V) ? edge_get_u(edge_wx) : edge_get_v(edge_wx);    
    const vertex vert_x = (rm_colors(etype_wx) == DIR_U_T_V) ? edge_get_v(edge_wx) : edge_get_u(edge_wx);    
    const short color_wx = (vert_w<vert_x) ? get_color_u_v(etype_wx) : get_color_v_u(etype_wx);


    //go over neighbours and find candidates
    if (vert_w != eu && vert_w != ev && vert_x != eu && vert_x != ev
        && maing.edges.find(edge_code(eu,vert_x)) == maing.edges.end()
        && maing.edges.find(edge_code(ev,vert_w)) == maing.edges.end()) { //vertices ok?
        for (unsigned long i = 0; i != maing.num_neighbours[vert_x]; ++i) {
            vertex y = maing.neighbours[vert_x][i];
            if (y!= eu && y != ev && maing.edges.find(edge_code(ev,y)) == maing.edges.end()) {
                const edgetype et_xy = maing.edges[edge_code(vert_x,y)];   
                const bool direction_xy_ok = (vert_x<y) ? (rm_colors(et_xy) == DIR_U_T_V) : (rm_colors(et_xy) == DIR_V_T_U);
          	    const short color_xy = (vert_x<y) ? get_color_u_v(et_xy) : get_color_v_u(et_xy);
          	    const bool xy_vcolor_ok = vertex_color_matters ? (maing.vertex_colors[eu] == maing.vertex_colors[y]) : true;
          	    const bool xy_ecolor_ok = edge_color_matters ? (color_eveu == color_xy) : true;	    
                if (direction_xy_ok && xy_vcolor_ok && xy_ecolor_ok) {
                    if (rand.trueWithProb(1.0 / (++num_rule5_candiates))) {
                       cand_w = vert_w;
                       cand_x = vert_x; 
                       cand_y = y;   
                       cand_color_xy = color_xy;                       
                       cand_color_wx = color_wx;                        
                    }                             
                }      
            }
        } 
    }//end vertices ok check
    //if candidates were found, perform the switch
    if (cand_w != NILLVERTEX && cand_x != NILLVERTEX && cand_y != NILLVERTEX) {
 /*          wxMessageBox(wxString("ex: ") << eu << " " << ev << " " << cand_w 
                                  << " " << cand_x << " " << cand_y 
                                  << " " << color_euev << " " << cand_color_wx
                                  << "  " << color_eveu << " " << cand_color_xy, "Error reading file", wxICON_ERROR | wxOK);
*/                                  
           ++rule5;
           const edge e_ux = edge_code(eu,cand_x),
                      e_wx = edge_code(cand_w,cand_x),
                      e_wv = edge_code(cand_w,ev),
                      e_xy = edge_code(cand_x,cand_y),
                      e_vy = edge_code(ev,cand_y);                                            
                                            
           //adjust edges in maingraph
           
           const edgetype et_ux = colored_edge(eu,cand_x,color_euev,color_eveu),
                          et_wv = colored_edge(ev,cand_w,0,cand_color_wx),
                          et_vy = colored_edge(ev,cand_y,cand_color_xy,0);
           
           maing.edges[e_ux] = et_ux;               
           maing.edges[e_wv] = et_wv;
           maing.edges[e_vy] = et_vy;                      

           maing.edges.erase(e);
           maing.edges.erase(e_wx);
           maing.edges.erase(e_xy);
                                                                  
           //adjust neighbors in maingraph
           replace_neighbour(maing,eu,ev,cand_x);
           replace_neighbour(maing,ev,eu,cand_y);
           add_neighbour(maing,ev,cand_w); 
           replace_neighbour(maing,cand_w,cand_x,ev);
           replace_neighbour(maing,cand_x,cand_y,eu);
           del_neighbour(maing,cand_x,cand_w);
           replace_neighbour(maing,cand_y,cand_x,ev);                      
           
           //adjust edgecontainer               
           EC.loggedReplaceElement(e,e_ux);           
           EC.loggedReplaceElement(e_wx,e_wv);
           EC.loggedReplaceElement(e_xy,e_vy);
       break;
    }
    
    }//end if ebagsize > 0
    
    
    
    //Try to apply rule 4
    //find candidates, choose random one
    //Perform if success and break edge exchange try
    vertex candidate_w = NILLVERTEX, candidate_x = NILLVERTEX;
    short candidate_ecolor_tow = 0, candidate_ecolor_wx = 0;
    double number_of_candidates = 0.0;
    
    //iterate over v neighbors
    for (unsigned long i = 0; i != maing.num_neighbours[ev]; ++i) {
        vertex w = maing.neighbours[ev][i];
        const edgetype et_vw = maing.edges[edge_code(ev,w)];   
        const bool direction_vw_ok = (ev<w) ? (rm_colors(et_vw) == DIR_U_T_V) : (rm_colors(et_vw) == DIR_V_T_U);
  	    const short color_vw = (ev<w) ? get_color_u_v(et_vw) : get_color_v_u(et_vw);
  	    const bool uw_vcolor_ok = vertex_color_matters ? (maing.vertex_colors[eu] == maing.vertex_colors[w]) : true;
  	    const bool uw_ecolor_ok = edge_color_matters ? (color_eveu == color_vw) : true;	    

        if ( direction_vw_ok  //check if single edge
             && uw_vcolor_ok //vertex colors ok?
             && uw_ecolor_ok ) {
             //neighbor w would be ok      
             for (unsigned long j = 0; j != maing.num_neighbours[w]; ++j) {
                 vertex x = maing.neighbours[w][j];
                 const edgetype et_wx = maing.edges[edge_code(w,x)];   
                 const bool direction_wx_ok = (w<x) ? (rm_colors(et_wx) == DIR_U_T_V) : (rm_colors(et_wx) == DIR_V_T_U);
  	             const short color_wx = (w<x) ? get_color_u_v(et_wx) : get_color_v_u(et_wx);
  	             const bool wx_vcolor = vertex_color_matters ? (maing.vertex_colors[ev] == maing.vertex_colors[x]) : true;
  	             const bool wx_ecolor = edge_color_matters ? (color_euev == color_wx) : true;	    

                 if (direction_wx_ok && x!=eu && x!=ev              
                     && wx_vcolor //vertex colors ok?
                     && wx_ecolor 
                     && (maing.edges.find(edge_code(eu,x)) == maing.edges.end() )  ) {
                    //We have a candidate
                    //This candidate is set with a certain probability
                    if (rand.trueWithProb(1.0 / (++number_of_candidates))) {
                       candidate_w = w; 
                       candidate_x = x;  
                       candidate_ecolor_tow = color_vw; 
                       candidate_ecolor_wx = color_wx;                                               
                    }                    
                 }
                   
             }      
        }
    }
    bool w_adjacent_v = true;
        
    //iterate over u neighbors
   for (unsigned long i = 0; i != maing.num_neighbours[eu]; ++i) {
        vertex w = maing.neighbours[eu][i];
        const edgetype et_uw = maing.edges[edge_code(eu,w)];   
        const bool direction_uw_ok = (eu<w) ? (rm_colors(et_uw) == DIR_U_T_V) : (rm_colors(et_uw) == DIR_V_T_U);
  	    const short color_uw = (eu<w) ? get_color_u_v(et_uw) : get_color_v_u(et_uw);
  	    const bool vw_vcolor = vertex_color_matters ? (maing.vertex_colors[ev] == maing.vertex_colors[w]) : true;
  	    const bool vw_ecolor = edge_color_matters ? (color_euev == color_uw) : true;	    
  	    
        if ( direction_uw_ok  //check if single edge
             && vw_vcolor //vertex colors ok?
             && vw_ecolor  ) {
             //neighbor w would be ok      
             for (unsigned long j = 0; j != maing.num_neighbours[w]; ++j) {
                 vertex x = maing.neighbours[w][j];
                 const edgetype et_wx = maing.edges[edge_code(w,x)];   
                 const bool direction_wx_ok = (w<x) ? (rm_colors(et_wx) == DIR_U_T_V) : (rm_colors(et_wx) == DIR_V_T_U);
  	             const short color_wx = (w<x) ? get_color_u_v(et_wx) : get_color_v_u(et_wx);
  	             const bool wx_vcolor = vertex_color_matters ? (maing.vertex_colors[eu] == maing.vertex_colors[x]) : true;
  	             const bool wx_ecolor = edge_color_matters ? (color_eveu == color_wx) : true;	    

                 if (direction_wx_ok && x!=eu && x!=ev              
                     && wx_vcolor //vertex colors ok?
                     && wx_ecolor 
                     && (maing.edges.find(edge_code(ev,x)) == maing.edges.end() )  ) {
                    //We have a candidate
                    //This candidate is set with a certain probability
                    if (rand.trueWithProb(1.0 / (++number_of_candidates))) {
                       w_adjacent_v = false;
                       candidate_w = w; 
                       candidate_x = x;     
                       candidate_ecolor_tow = color_uw; 
                       candidate_ecolor_wx = color_wx;                                            
                    }                    
                 }
                   
             }      
        }
    }    
    
    if (candidate_w != NILLVERTEX && candidate_x != NILLVERTEX) { //candiates found
       ++rule4;
       if (w_adjacent_v) {
           //adjust edges in maingraph
           maing.edges.erase(edge_code(candidate_w,candidate_x));  
           const edgetype new_type_e = (eu < ev) ? colored_v_u(DIR_V_T_U,candidate_ecolor_tow) 
                                                 : colored_u_v(DIR_U_T_V,candidate_ecolor_tow);         
           maing.edges[e] = new_type_e;
           const edgetype new_type_vw = ((eu<ev)&&(candidate_w<ev)||(eu>ev)&&(candidate_w>ev)) ? et : reverse(et);
           maing.edges[edge_code(ev,candidate_w)] = new_type_vw;
           const edgetype new_type_ux = (eu < candidate_x) ? colored_u_v(DIR_U_T_V,candidate_ecolor_wx) 
                                                           : colored_v_u(DIR_V_T_U,candidate_ecolor_wx);         
           maing.edges[edge_code(eu,candidate_x)] = new_type_ux;
                                                                  
           //adjust neighbors in maingraph
           add_neighbour(maing,eu,candidate_x);
           del_neighbour(maing,candidate_w,candidate_x);     
           replace_neighbour(maing,candidate_x,candidate_w,eu);
           
           //adjust edgecontainer
           EC.loggedSwapElements(e,edge_code(ev,candidate_w));
           EC.loggedReplaceElement(edge_code(candidate_w,candidate_x),edge_code(eu,candidate_x));           
           
       } else {
           maing.edges.erase(edge_code(candidate_w,candidate_x));  
           const edgetype new_type_e = (eu < ev) ? colored_u_v(DIR_U_T_V,candidate_ecolor_tow)
                                                 : colored_v_u(DIR_V_T_U,candidate_ecolor_tow);         
           maing.edges[e] = new_type_e;
           const edgetype new_type_uw = ((ev<eu)&&(candidate_w<eu)||(ev>eu)&&(candidate_w>eu)) ? et : reverse(et);
           maing.edges[edge_code(eu,candidate_w)] = new_type_uw;
           const edgetype new_type_vx = (ev < candidate_x) ? colored_u_v(DIR_U_T_V,candidate_ecolor_wx) 
                                                           : colored_v_u(DIR_V_T_U,candidate_ecolor_wx);         
           maing.edges[edge_code(ev,candidate_x)] = new_type_vx;
                                                                  
           //adjust neighbors in maingraph
           add_neighbour(maing,ev,candidate_x);
           del_neighbour(maing,candidate_w,candidate_x);     
           replace_neighbour(maing,candidate_x,candidate_w,ev);
           
           //adjust edgecontainer
           EC.loggedSwapElements(e,edge_code(eu,candidate_w));
           EC.loggedReplaceElement(edge_code(candidate_w,candidate_x),edge_code(ev,candidate_x));           
    
       }
       break;
    }

//-----------------------End new global BIDIR rules                 
}
//-------------------------------------------------            
                            
                             if (bag_size == 1)
                             { total_tries+= (num_tries-1);  break;}

							 //choose random exchange partner
                             bag_index random_index = (++rand) % bag_size;
                             while (random_index == current_index) {
                                   random_index = (++rand) % bag_size;
                             }

                             //check if switch can be performed
                             
                             const edge f = EC.getIndexedElement(random_index);
							 
                             vertex fu = edge_get_u(f);
                             vertex fv = edge_get_v(f);

                             const edgetype ft = maing.edges[f];

							 short color_fufv = get_color_u_v(ft);
							 short color_fvfu = get_color_v_u(ft);

							 bool reverse_edge_f = (rm_colors(et) == rm_colors(ft)) ? false : true;
							 if (edge_color_matters &&  color_euev != color_fufv)
								  reverse_edge_f = true;	
							 if (reverse_edge_f)
								{ vertex tmp = fu; fu = fv; fv = tmp;
							      tmp = color_fufv; color_fufv = color_fvfu; color_fvfu = tmp; }

							 const edge new_e = edge_code(eu,fv);
							 const edge new_f = edge_code(fu,ev);

							 if ( maing.edges.find(new_e) == maing.edges.end() 
								 && maing.edges.find(new_f) == maing.edges.end()
                                 && eu != fu && eu != fv && ev != fu && ev != fv) {

								 //perform on EC (bagID remains here)
								 EC.loggedSetIndexedElement(current_index,new_e);
								 EC.loggedSetIndexedElement(random_index,new_f);

								 //perform on HashMap
								 edgetype new_et = 0; 
								 if (eu < fv) {
                                    if (color_euev > 0) {
                                       new_et |= DIR_U_T_V;
                                       set_color_u_v(new_et, color_euev);               
                                    }
                                    if (color_eveu > 0) {
                                       new_et |= DIR_V_T_U;
                                       set_color_v_u(new_et, color_eveu);               
                                    }
								 } else {
                                    if (color_euev > 0) {
                                       new_et |= DIR_V_T_U;
                                       set_color_v_u(new_et, color_euev);               
                                    }
                                    if (color_eveu > 0) {
                                       new_et |= DIR_U_T_V;
                                       set_color_u_v(new_et, color_eveu);               
                                    }
								 }

								 edgetype new_ft = 0; 
								 if (fu < ev) {
                                    if (color_fufv > 0) {
                                       new_ft |= DIR_U_T_V;
                                       set_color_u_v(new_ft, color_fufv);               
                                    }
                                    if (color_fvfu > 0) {
                                       new_ft |= DIR_V_T_U;
                                       set_color_v_u(new_ft, color_fvfu);               
                                    }
								 } else {
                                    if (color_fufv > 0) {
                                       new_ft |= DIR_V_T_U;
                                       set_color_v_u(new_ft, color_fufv);               
                                    }
                                    if (color_fvfu > 0) {
                                       new_ft |= DIR_U_T_V;
                                       set_color_u_v(new_ft, color_fvfu);               
                                    }
								 }

								 maing.edges.erase(e);
								 maing.edges.erase(f);
								 maing.edges[new_e] = new_et;
								 maing.edges[new_f] = new_ft;
								 replace_neighbour(maing,eu,ev,fv);
								 replace_neighbour(maing,ev,eu,fu);
								 replace_neighbour(maing,fu,fv,ev);
								 replace_neighbour(maing,fv,fu,eu);
                                                                  								 
								 ++total_success;
								 break;//Break try-loop

							 }//if exchange can be performed
                         }
    			         //End switch attempt
                   }
             }
         }
         
         //Update Maingraph
         build_graph(maing);
         
//wxMessageBox(wxString("  ex: ") << rule4 << " " << rule5, "Error reading file", wxICON_ERROR | wxOK);
         
         
    break; // break case GLOBAL_CONST

    case LOCAL_CONST: 

         //Perform Switching
         for (int exchanges = 0; exchanges != num_exchanges; ++exchanges) {
    	     EC.resetIteration();
//wxMessageBox(wxString("  ex: ") << exchanges , "Error reading file", wxICON_ERROR | wxOK);

    	     while (EC.hasNextBag()) {
                   EC.goNextBag();
                   bagID bag_id = EC.getCurrentBag();
                   bag_index bag_size = EC.getBagSize(bag_id);
    		       while (EC.hasNextElement()) {
    			         EC.goNextElement();
    			         bag_index current_index = EC.getCurrentIndex();
    			         //Begin switch attempt
    			         for (int tries = 0; tries != num_tries; tries++) {
                             ++total_tries;
                             if (bag_size == 1)
                             { total_tries+= (num_tries-1);  break;}

							 //choose random exchange partner
                             bag_index random_index = (++rand) % bag_size;
                             while (random_index == current_index) {
                                   random_index = (++rand) % bag_size;
                             }

                             //check if switch can be performed

                             const edge e = EC.getCurrentElement();
                             const edge f = EC.getIndexedElement(random_index);
							 
                             const vertex eu = edge_get_u(e);
                             const vertex ev = edge_get_v(e);
                             vertex fu = edge_get_u(f);
                             vertex fv = edge_get_v(f);

                             const edgetype et = maing.edges[e],
                                            ft = maing.edges[f];

							 const short color_euev = get_color_u_v(et);
							 const short color_eveu = get_color_v_u(et);
							 short color_fufv = get_color_u_v(ft);
							 short color_fvfu = get_color_v_u(ft);

							 bool reverse_edge_f = (rm_colors(et) == rm_colors(ft)) ? false : true;
							 if (edge_color_matters &&  color_euev != color_fufv)
								  reverse_edge_f = true;	
							 if (reverse_edge_f)
								{ vertex tmp = fu; fu = fv; fv = tmp;
							      tmp = color_fufv; color_fufv = color_fvfu; color_fvfu = tmp; }

							 const edge new_e = edge_code(eu,fv);
							 const edge new_f = edge_code(fu,ev);

							 if ( maing.edges.find(new_e) == maing.edges.end() 
								 && maing.edges.find(new_f) == maing.edges.end()
                                 && eu != fu && eu != fv && ev != fu && ev != fv) {

								 //perform on EC (bagID remains here)
								 EC.setIndexedElement(current_index,new_e);
								 EC.setIndexedElement(random_index,new_f);

								 //perform on HashMap
								 edgetype new_et = 0; 
								 if (eu < fv) {
                                    if (color_euev > 0) {
                                       new_et |= DIR_U_T_V;
                                       set_color_u_v(new_et, color_euev);               
                                    }
                                    if (color_eveu > 0) {
                                       new_et |= DIR_V_T_U;
                                       set_color_v_u(new_et, color_eveu);               
                                    }
								 } else {
                                    if (color_euev > 0) {
                                       new_et |= DIR_V_T_U;
                                       set_color_v_u(new_et, color_euev);               
                                    }
                                    if (color_eveu > 0) {
                                       new_et |= DIR_U_T_V;
                                       set_color_u_v(new_et, color_eveu);               
                                    }
								 }

								 edgetype new_ft = 0; 
								 if (fu < ev) {
                                    if (color_fufv > 0) {
                                       new_ft |= DIR_U_T_V;
                                       set_color_u_v(new_ft, color_fufv);               
                                    }
                                    if (color_fvfu > 0) {
                                       new_ft |= DIR_V_T_U;
                                       set_color_v_u(new_ft, color_fvfu);               
                                    }
								 } else {
                                    if (color_fufv > 0) {
                                       new_ft |= DIR_V_T_U;
                                       set_color_v_u(new_ft, color_fufv);               
                                    }
                                    if (color_fvfu > 0) {
                                       new_ft |= DIR_U_T_V;
                                       set_color_u_v(new_ft, color_fvfu);               
                                    }
								 }

								 maing.edges.erase(e);
								 maing.edges.erase(f);
								 maing.edges[new_e] = new_et;
								 maing.edges[new_f] = new_ft;

								 ++total_success;
								 break;//Break try-loop

							 }//if exchange can be performed
                         }
    			         //End switch attempt
                   }
             }
         }
         
         //Update Maingraph
         
         build_graph(maing);
         
         
    break; // break case LOCAL_CONST
    }
           
    


    return double (clock() - start_time) / CLOCKS_PER_SEC;                     
} 





/*
double randomize_graph(maingraph & maing, short random_type,
                     int num_exchanges, int num_tries, long & total_tries,
                     long & total_success)
{
    randlib::rand rand(time(NULL)); // Init Random-Generator
    edge* edgelist; // The edges are written into this array
    long* open_positions = NULL; // Stores the free places in the edgelist array.
    long op_index = 0;    // The index one behind the last open position.

    edge e, p, ne1, ne2, ne3, single, undir, third;
    // e is the edge, which needs an exchange partner, p is this (random) found partner,
    // ne1-3 are newly created edges (often it is first checked, whether those edges already exist)
    // single, undir and third are used for the global complex rules.
    // you start with an single/directed edge, pick an undirected edge, and a third directed edge.
    vertex eu, ev, pu, pv, unu, unv, tu, tv, begin_path, middle, end_path, other;
    // the first are the u and v of the above defined vertices.
    // begin_path, middle, end_path and other are again used for the complex rules.
    // there is an path between three vertices, and the fourth is the "other" vertex.
    edgetype et, pt, est, pst, ne1t, ne2t, st, tt;
    // Types of the edges above: st is type of edge "single", tt type of the edge "third".
    // est and pst are used for "no regard"-graphs: edge split type and partner split type
    hash_map < edge, edgetype >::const_iterator find_ne1, find_ne2;
    // those store the values of finding ne1 and ne2 in the hashmap.
    bool success, found_path;
    // success is set if an exchange was made, found_path, if an path between the vertices is found
    // These longs store the numbers the edges have in the edgelist array.
    long p_number, undir_number, single_number;
    // These vars do some statistics
    long complex4used = 0, complex4tried = 0,
         complex5used = 0, complex5tried = 0, failedbidir = 0;

    // In NO_REGARD-Graphs the number of edges can increase due to an edge-split
    if (random_type == NO_REGARD){
        edgelist = new edge[maing.m+maing.num_undir_edges];
        open_positions = new long[maing.m];
    } else
        edgelist = new edge[maing.m];

	clock_t start_time(clock()); // Start the time measurement

    // Turn edges hashmap into an array. It is sorted: First dir edges, then undir edges
    { long i=0, j=maing.num_dir_edges;
    for (hash_map < edge, edgetype >::const_iterator iter = maing.edges.begin();
         iter != maing.edges.end(); ++iter) {
        if (iter->second == UNDIR_U_V){
            edgelist[j] = iter->first;
            ++j;
        } else {
            edgelist[i] = iter->first;
            ++i;
        }
    }
    }

    for (int exchange_ctr = 0; exchange_ctr < num_exchanges; ++exchange_ctr)
    {
        // for each edge: try to find a partner
        for (unsigned long e_number = 0; e_number < maing.m; ++e_number){
            e = edgelist[e_number];  // select edge according to counter
            if (e == 0xFFFFFFFFFFFFFFFFULL) { // if the edge has been erased (which can happen in NO_REGARD-Graphs only)
                continue;                     // continue the for-loop at the next edge.
            }
            eu = edge_get_u(e);  // do the lookups for e
            ev = edge_get_v(e);
            et = maing.edges[e];
            success = false;       // reset success flag
            // Try to use a rule until one is successful or num_tries is exceeded
            for (int tries_ctr = 0; (!success) && (tries_ctr < num_tries); ++tries_ctr){
                found_path = false;  // reset found_path flag
                ++total_tries;       // count tries for the statistic
                switch (random_type){
// NO_REGARD ------------
                case NO_REGARD:
                    p_number = ++rand % maing.m;
                    p = edgelist[p_number]; // find partner edge
                    while (p == 0xFFFFFFFFFFFFFFFFULL) { // if the edge has been erased
                        p_number = ++rand % maing.m;     // chose a new one.
                        p = edgelist[p_number];
                    }
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];
                    if (pu != eu && pv != ev && pu != ev && pv != eu){ // ensure that vertices are disjunct
                        // In case e or p are bidir. choose only one direction.
                        est = (et == UNDIR_U_V) ? (++rand % 2)+1 : et; // split undirected edge
                        pst = (pt == UNDIR_U_V) ? (++rand % 2)+1 : pt;
                        if (est == pst) { // Determine directions and type of new edges
                            if (eu < pv){ // This creates ne1 as an new edge from eu to pv
                                ne1 = new_edge(eu, pv);
                                ne1t = est;
                            } else {      // In case pv is smaller than eu the direction has to be reversed.
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(est);
                            }
                            if (pu < ev){ // ne2 and its type are defined going from pu to ev.
                                ne2 = new_edge(pu, ev);
                                ne2t = pst;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pst);
                            }
                        } else { // est != pst In case the types of e and p are unequal,
                                 // p is mentally turned, so that the new edges now
                                 // lead from eu to pu and from pv to ev.
                            if (eu < pu){
                                ne1 = new_edge(eu, pu);
                                ne1t = est;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(est);
                            }
                            if (pv < ev){
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pst); // because p is handled as if reversed, we have to reverse the type of the new edge.
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pst; // This is really reverse(reverse(pst))!
                            }
                        } // end if types equal
                        find_ne1 = maing.edges.find(ne1); // Look up the edges we want to create.
                        find_ne2 = maing.edges.find(ne2);
                        // Check if they are not already existing.
                        if ((find_ne1 == maing.edges.end() || find_ne1->second == reverse(ne1t)) &&
                            (find_ne2 == maing.edges.end() || find_ne2->second == reverse(ne2t))){
                            // Add to hashmap and update edgelist.
                            // The types are not regarded, as the list does not need to be sorted for this random-type
                            if (maing.edges.find(ne1) == maing.edges.end()){ // ne1 does not exist:
                                maing.edges[ne1] = ne1t;         // Create entry for ne1
                                if (et == est){                  // In this case, e has to be erased:
                                    edgelist[e_number] = ne1;    // Overwrite e with ne1
                                } else {                         // in this case, e should stay in the list:
                                    if (op_index == 0){          // There are no open positions:
                                        edgelist[maing.m] = ne1; // Add ne1 at the end of list.
                                        ++maing.m;
                                    } else {                     // There are open positions:
                                        --op_index;              // Write ne2 in one of them.
                                        edgelist[open_positions[op_index]] = ne1;
                                    }
                                }
                            }
                            else {                               // ne1 exists:
                                maing.edges[ne1] = UNDIR_U_V;    // ne1 is now undirected.
                                if (et == est){                  // In this case, e has to be erased:
                                    edgelist[e_number] = 0xFFFFFFFFFFFFFFFFULL; // This is an kind of erase...
                                    open_positions[op_index] = e_number; // The position of e is now open.
                                    ++op_index;
                                }
                                // If (et != est) e should remain in the list and ne1 is already in the list:
                                // There is no work at all.
                            }

                            if (maing.edges.find(ne2) == maing.edges.end()){ // ne2 does not exist:
                                maing.edges[ne2] = ne2t;         // Create entry for ne2
                                if (pt == pst){                  // In this case, p has to be erased:
                                    edgelist[p_number] = ne2;    // Overwrite p with ne2
                                } else {                         // in this case, p should stay in the list:
                                    if (op_index == 0){          // There are no open positions:
                                        edgelist[maing.m] = ne2; // Add ne2 at the end of list.
                                        ++maing.m;
                                    } else {                     // There are open positions:
                                        --op_index;              // Write ne2 in one of them.
                                        edgelist[open_positions[op_index]] = ne2;
                                    }
                                }
                            }
                            else {                               // ne2 exists:
                                maing.edges[ne2] = UNDIR_U_V;    // ne2 is now undirected.
                                if (pt == pst){                  // In this case, p has to be erased:
                                    edgelist[p_number] = 0xFFFFFFFFFFFFFFFFULL; // This is an kind of erase...
                                    open_positions[op_index] = p_number; // The position of p is now open.
                                    ++op_index;
                                }
                                // If (pt != pst) p should remain in the list and ne2 is already in the list:
                                // There is no work at all.
                            }

                            // Erase from Hashmap:
                            if (et == est)                    // In this case e has to be erased.
                                maing.edges.erase(e);
                            else {                            // et != est: e has to stay in the map
                                maing.edges[e] = reverse(est);// Set e only to the not switched dir.
                            }

                            if (pt == pst)                    // In this case p has to be erased.
                                maing.edges.erase(p);         // pt != pst: p has to stay in the map
                            else {
                                maing.edges[p] = reverse(pst);// Set p only to the not switched dir.
                            }
                            success=true;
                        }  // end if edges do not exist
                    } // end if vertices unequal
                break;
                case GLOBAL_CONST:
// GLOBAL_CONST
// Complex 4 vertex rule
                    // If there are undir edges and directed edges, with an 50% probability:
                    // Use complex 4-vertex-rule (probability can be adjusted)
                    if (maing.num_undir_edges != 0 && maing.num_dir_edges != 0
                        && ++rand % 2 == 0){
                        ++complex4tried; // Count tries of this rule
                        // If edge e is undirected: randomly select directed edge
                        if (et == UNDIR_U_V){
                            single_number = ++rand % maing.num_dir_edges;
                            undir_number = e_number; // store the numbers to be able to update edgelist array
                            single = edgelist[single_number];
                            st = maing.edges[single]; // lookup the single edge
                            // define the middle and the end of the path according to the direction of "single"
                            middle = (st == DIR_U_T_V) ? edge_get_u(single) : edge_get_v(single);
                            end_path = (st == DIR_U_T_V) ? edge_get_v(single) : edge_get_u(single);
                            undir = e; // edge e is from now on called "undir"
                            unu = eu;  // its vertices are unu and unv
                            unv = ev;
                        } else {  // If edge e is directed randomly select undirected edge
                            undir_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                            single_number = e_number; // store the numbers to be able to update edgelist array
                            undir = edgelist[undir_number]; // the undirected edge is looked up
                            unu = edge_get_u(undir);        // and its vertices defined
                            unv = edge_get_v(undir);
                            single = e;                     // edge e is from now on called "single"
                            st = et;                        // its type is st (It has to be one of the dir. types)
                            middle = (st == DIR_U_T_V) ? eu : ev;    // middle and end_path are defined
                            end_path =  (st == DIR_U_T_V) ? ev : eu; // according to the direction of edge "single"
                        }
                        // Ensure that no two vertices are the same:
                        if (unu != middle && unv != end_path && unu != end_path && unv != middle){
                            begin_path = unu;   // first try: unu is the begin of the path
                            other = unv;        // unv is the "other" vertex, which is not in the path
                            // Define u and v of the "third" edge
                            if (unu < middle) {
                               tu = begin_path;
                               tv = middle;
                            } else {
                               tu = middle;
                               tv = begin_path;
                            }
                            third = new_edge(tu, tv); // Define the third edge.

                            if (maing.edges.find(third) != maing.edges.end()) { // third should exist
                                tt = maing.edges[third];
                                if (tt == DIR_U_T_V && tu == begin_path) found_path = true; // third is pointing into the right direction
                                if (tt == DIR_V_T_U && tv == begin_path) found_path = true; // a path is found.
                            }
                            if (!found_path){ // If no path could be found with unu as begin path:
                                begin_path = unv; // Use unv as begin
                                other = unu;      // unu then is the other vertex, not in the path
                                if (unv < middle) {  // same as above
                                    tu = begin_path;
                                    tv = middle;
                                } else {
                                    tu = middle;
                                    tv = begin_path;
                                }
                                third = new_edge(tu, tv); // third is looked up again
                                if (maing.edges.find(third) != maing.edges.end()) {
                                   tt = maing.edges[third];  // as above: Is third pointing into the right direction?
                                   if (tt == DIR_U_T_V && tu == begin_path) found_path = true;
                                   if (tt == DIR_V_T_U && tv == begin_path) found_path = true;
                               }
                            } // end if use unv as begin
                            if (found_path) { // if an path was found:
                                if (other < end_path){ // Define ne1 as pointing from other to end path.
                                    ne1 = new_edge(other, end_path);
                                    ne1t = DIR_U_T_V;
                                 } else {
                                    ne1 = new_edge(end_path, other);
                                    ne1t = DIR_V_T_U;
                                }
                                // if this new edge ne1 does not exist: begin the switching:
                                if (maing.edges.find(ne1) == maing.edges.end()){
                                    maing.edges[ne1] = ne1t;        // ne1 is created
                                    maing.edges[third] = UNDIR_U_V; // the third edge is nor undirected
                                    if (begin_path < other){  // The undirected edge is now pointing from begin_path to other
                                        maing.edges[undir] = DIR_U_T_V;
                                    } else {
                                        maing.edges[undir] = DIR_V_T_U;
                                    }

                                    maing.edges.erase(single); // the single edge is erased.
                                    // ne1 is an dir edge, therefore it can take the place of single in the edgelist
                                    edgelist[single_number]= ne1;
                                    // undir has become directed, therefore it
                                    // takes the place of third, which used to be directed
                                    find_and_replace(edgelist, third, undir);
                                    // third is now undirected and placed at "undir"s position.
                                    edgelist[undir_number] = third;
                                    success = true;
                                    ++complex4used; // The rule was used successfully
                                } // end if switching
                            } // end if found path
                        } // end if vertices unequal

// Easy rules -------------
                    } else { // Use easy rules or 5-vertex-rule

                    // If edge e is undirected: randomly select undirected edge
                    if (et == UNDIR_U_V){ // randomly select undirected edge
                       p_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                    } else { // et = DIR: select another directed edge
                       p_number = ++rand % maing.num_dir_edges;
                    } // end if et = UNDIR
                    p = edgelist[p_number]; // find partner edge
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];    // lookup type of partner edge
                    if (pu != eu && pv != ev && pu != ev && pv != eu){
                        if (et == pt) { // Determine directions and type of new edges
                            if (eu < pv){  // These are the same procedures as in "NO REGARD" above
                                ne1 = new_edge(eu, pv);  // New edge ne1 from eu to pv
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(et);
                            }
                            if (pu < ev){                // New edge ne2 from pu to ev
                                ne2 = new_edge(pu, ev);
                                ne2t = pt;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pt);
                            }
                        } else { // et != pt Two directed edges point in opposite directions.
                        // as above p is then mentally turned.
                            //if (et==UNDIR_U_V) {} // Assertion: This should never ever happen.
                            if (eu < pu){               // New edge ne1 from eu to pu
                                ne1 = new_edge(eu, pu);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(et);
                            }
                            if (pv < ev){               // New edge ne2 from pv to ev
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pt);
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pt;
                            }
                        } // end if types equal

                        // Look up the edges we want to create
                        find_ne1 = maing.edges.find(ne1);
                        find_ne2 = maing.edges.find(ne2);
                        // Check whether they are not already existing:
                        if (find_ne1 == maing.edges.end() &&
                            find_ne2 == maing.edges.end()) {
                            maing.edges.erase(e);   //Update Hashmap
                            maing.edges.erase(p);   // e and p are erased.
                            maing.edges[ne1] = ne1t; // ne1 and ne2 added with their types
                            maing.edges[ne2] = ne2t;
                            edgelist[e_number]= ne1;         // Update edgelist
                            edgelist[p_number] = ne2; // ne1 takes the position of e, ne2 the position of p
                            success=true;
// 5 Vertex Rule ---------------
                        // There is an connection between the vertices: try to use 5-vertex-rule
                        // Edge e has to be dir as it shall be part of the path.
                        // There has to be at least one undir edge and at least two dir edges, to use the rule
                        } else if (et != UNDIR_U_V && maing.num_undir_edges != 0 && maing.num_dir_edges > 1) {
                            ++complex5tried; // Count tries of this rule
                            if (find_ne1 != maing.edges.end()) { // If edge ne1 exists:
                                // try to use it as the third edge for the switching rule
                                // The name third comes from the 4 vertex rule, where this is really the third edge to be picked.
                                third = ne1;
                                // edges changed are: e, undir (defined later), and third
                                tt = find_ne1->second; // the type of third according to the hashtable.
                                if (et == DIR_U_T_V){
                                   // if edge e points from u to v, it is the second edge in the path,
                                   // because ne1 has eu as one vertex, so eu is the middle vertex of the path.
                                   if (edge_get_u(third) == eu){ // If third and e have the same "u"
                                       if (tt == DIR_V_T_U){ // and the directions are ok
                                           found_path = true;                // we have found an path
                                           begin_path = edge_get_v(third);   // begin, middle and end are set accordingly
                                           middle = eu;
                                           end_path = ev;

                                       }
                                   } else { // edge_get_v(third) == eu
                                       // third and e do not have the same u =>
                                       // third has to point in the opposite direction for the path to be correct.
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = edge_get_u(third);
                                           middle = eu;
                                           end_path = ev;
                                       }
                                   }
                                } else { // et == DIR_V_T_U
                                   // if edge e points from v to u, it is the first edge in the path.
                                   if (edge_get_u(ne1) == eu){
                                       if (find_ne1->second == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = ev;
                                           middle = eu;
                                           end_path = edge_get_v(third);
                                       }
                                   } else { // edge_get_v(third) == eu
                                       if (find_ne1->second == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = ev;
                                           middle = eu;
                                           end_path = edge_get_u(third);
                                       }
                                   }
                                } // end if lookup of et
                            } // end if "ne1 exists"
                            if (!found_path && find_ne2 != maing.edges.end()) {
                            // If we haven't already found a path and ne2 exists,
                            // we are trying to use it as "third" edge.
                                third = ne2; // edges changed are: e, undir (defined later), and third
                                // ne2 certainly has ev as one vertex.
                                tt = find_ne2->second; // Type of third according to the hashtable
                                if (et == DIR_U_T_V){  // These things are going the same way, as with ne1 above
                                   if (edge_get_u(third) == ev){
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = eu;
                                           middle = ev;
                                           end_path = edge_get_v(third);
                                       }
                                   } else { // edge_get_v(third) == eu
                                       if (tt == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = eu;
                                           middle = ev;
                                           end_path = edge_get_u(third);
                                       }
                                   }
                                } else { // et == DIR_V_T_U
                                   if (edge_get_u(third) == ev){
                                       if (tt == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = edge_get_v(third);
                                           middle = ev;
                                           end_path = eu;
                                       }
                                   } else { // edge_get_v(ne2) == ev
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = edge_get_u(third);
                                           middle = ev;
                                           end_path = eu;
                                       }
                                   }
                                } // end if lookup of et
                            } // end if ne2 exists
                        } // end if new edges do not exist
                     } // end if vertices unequal (These brackets finish the easy rule usage)
                     if (found_path){ // Path was found in section above: Using 5-vertex-rule
                     // Try to find an bidir edge fitting to the path num_tries times.
                     for (int tries_ctr2 = 0; (!success) && (tries_ctr2 < num_tries); ++tries_ctr2){
                         // randomly select an undir edge, number is stored once again for the edgelist update.
                         undir_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                         undir = edgelist[undir_number];
                         unu = edge_get_u(undir);
                         unv = edge_get_v(undir);
                         // if vertices are unequal
                         if (unu != begin_path && unu != middle && unu != end_path &&
                             unv != begin_path && unv != middle && unv != end_path) {
                             if (begin_path < unv) { // The new edge ne1 is pointing from begin_path to unv.
                                 ne1 = new_edge(begin_path, unv);
                                 ne1t = DIR_U_T_V;
                             } else {
                                 ne1 = new_edge(unv, begin_path);
                                 ne1t = DIR_V_T_U;
                             }
                             if (unv < end_path) {  // The new edge ne2 is pointing from unv to end_path.
                                 ne2 = new_edge(unv, end_path);
                                 ne2t = DIR_U_T_V;
                             } else {
                                 ne2 = new_edge(end_path, unv);
                                 ne2t = DIR_V_T_U;
                             }
                             if (unu < middle) {   // The new edge ne3 is undirected between unu and middle.
                                 ne3 = new_edge(unu, middle);
                             } else {
                                 ne3 = new_edge(middle, unu);
                             }
                             // Those three edges we want to create should not exist:
                             if (maing.edges.find(ne1) == maing.edges.end() &&
                                 maing.edges.find(ne2) == maing.edges.end() &&
                                 maing.edges.find(ne3) == maing.edges.end()) {
                                 // start the 5-vertice switching
                                 maing.edges[ne1] = ne1t;      // Insert ne1-3 with their types
                                 maing.edges[ne2] = ne2t;
                                 maing.edges[ne3] = UNDIR_U_V;
                                 maing.edges.erase(e);         // Erase the edges we worked with:
                                 maing.edges.erase(undir);     // e, undir and third
                                 maing.edges.erase(third);
                                 // third was undirected, ne1 is undir, so it can replace third.
                                 find_and_replace(edgelist, third, ne1);
                                 // e was undirected, ne2 ist undir, so it can replace e.
                                 edgelist[e_number] = ne2;
                                 // ne3 is undirected, so it can replace undir.
                                 edgelist[undir_number] = ne3;
                                 success = true;
                             } // end if new edges do not exist
                         } // end if vertices unequal
                         } // end for tries to find bidir edge
                         if (success) ++complex5used; // count the usage of this rule.
                         // if we were here, but did not succeed, we had an path, but no bidir edge for it.
                         else ++failedbidir;
                     } // end if found_path
                     } // end if which rule to use
                break;
// LOCAL_CONST
                case LOCAL_CONST: // The local const rules are exactly the same as the easy global const rules.
                                  // so comments become rare here.
                    if (et == UNDIR_U_V){ // randomly select undirected edge
                       p_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                    } else { // et = DIR
                       p_number = ++rand % maing.num_dir_edges;
                    } // end if et = UNDIR
                    p = edgelist[p_number]; // find partner edge
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];
                    if (pu != eu && pv != ev && pu != ev && pv != eu){
                        if (et == pt) { // Determine directions and type of new edges
                            if (eu < pv){
                                ne1 = new_edge(eu, pv);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(et);
                            }
                            if (pu < ev){
                                ne2 = new_edge(pu, ev);
                                ne2t = pt;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pt);
                            }
                        } else { // et != pt
                            //if (et==UNDIR_U_V)  {}// This should never ever happen.
                            if (eu < pu){
                                ne1 = new_edge(eu, pu);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(et);
                            }
                            if (pv < ev){
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pt);
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pt;
                            }
                        } // end if types equal
                        if (maing.edges.find(ne1)== maing.edges.end() &&
                            maing.edges.find(ne2)== maing.edges.end()) {
                            maing.edges.erase(e);   //Update Hashmap
                            maing.edges.erase(p);
                            maing.edges[ne1] = ne1t;
                            maing.edges[ne2] = ne2t;
                            edgelist[e_number]= ne1;         // Update edgelist
                            edgelist[p_number] = ne2;
                            success=true;
                        } // end if new edges do not exist
                     } // end if vertices unequal
                break;
                } // end switch
            } // end for tries
            if (success) total_success++; // count the amount of successful exchanges
        } // end for edges
    } // end for exchanges

    delete[] edgelist;
    if (random_type == NO_REGARD)
     delete[] open_positions;



    build_graph(maing); // Update neighbour arrays
    return double (clock() - start_time) / CLOCKS_PER_SEC;


}*/
