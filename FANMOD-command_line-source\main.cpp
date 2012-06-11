#include <stdlib.h>
#include <time.h>
#include <iostream>
using std::endl;
#include <fstream>
#include <sstream>
#include <vector>
using std::vector;
#include <stdexcept>
using std::domain_error;
#include <string>
using std::string;

using namespace std;


#include "hashmap.h"
#include "maingraph.hpp"
#include "output.hpp"
#include "random.hpp"


// input parameters

// algorithm options
// 1- subgraph (motif) size [default: 3]
// 2- # of samples used to determine approx. # of subgraphs [100000]
// 3- full enumeration? 1(yes)/0(no) [1]
// sampling probabilities will be input at end of command line, variable #!

// input file
// 4- infile name
// 5- directed? 1(yes)/0(no) [1]
// 6- colored vertices? 1(yes)/0(no) [0]
// 7- colored edges? 1(yes)/0(no) [0]

// random networks
// 8- random type: 0(no regard)/1(global const)/2(local const) [2]
// 9- regard vertex colors? 1(yes)/0(no) [0]
// 10- regard edge colors? 1(yes)/0(no) [0]
// 11- reestimate subgraph number? 1(yes)/0(no) [0]
// 12- # of random networks [1000]
// 13- # of exchanges per edge [3]
// 14- # of exchange attempts per edge [3]

// output file
// 15- outfile name
// 16- outfile format 1(ASCII - human readable)/0(CSV - for easy parsing) [1]
// 17- create dumpfile? 1(yes)/0(no) [0]

// sampling probabilities (if NOT full enumeration)
// 18->25 (depending on subgraph size) [.5]


// utility function
uint64 calc_expected_samples(double ntrees, const short & G_N, const double *prob);


int main(int argc, char** argv){
           
    char num_to_lett[] = {'0','1','2','3','4','5','6','7',
		'8','9','A','B','C','D','E','F', 
		'#','#','#','#','#','#','#','#',}; //last line for safety only
        
    // time counter
    time_t oldtime, newtime;
    struct tm * timeinfo;
    time (&oldtime); // begin time
    timeinfo = localtime (&oldtime);

    // command line arguments errors array (in order to save them and write to 
    // log file later)
    string *errors_array[100];
    for (int i=0; i<100; i++) {errors_array[i]=NULL;}
    
    if (argc<18) {
       cout << "Usage: fanmod <(at least) 17 basic input parameters as specified>\n\ninput parameters:\n-----------------\n\nalgorithm options\n1- subgraph (motif) size [default: 3]\n2- # of samples used to determine approx. # of subgraphs [100000]\n3- full enumeration? 1(yes)/0(no) [1]\nsampling probabilities will be input at end of command line, variable #!\n\ninput file\n4- infile name\n5- directed? 1(yes)/0(no) [1]\n6- colored vertices? 1(yes)/0(no) [0]\n7- colored edges? 1(yes)/0(no) [0]\n\nrandom networks\n8- random type: 0(no regard)/1(global const)/2(local const) [2]\n9- regard vertex colors? 1(yes)/0(no) [0]\n10- regard edge colors? 1(yes)/0(no) [0]\n11- reestimate subgraph number? 1(yes)/0(no) [0]\n12- # of random networks [1000]\n13- # of exchanges per edge [3]\n14- # of exchange attempts per edge [3]\n\noutput file\n15- outfile name\n16- outfile format 1(ASCII - human readable)/0(CSV - for easy parsing) [1]\n\n17- create dumpfile? 1(yes)/0(no) [0]\n\nsampling probabilities (if NOT full enumeration)\n18->25 (depending on subgraph size) [.5]\n";
       exit (-1);                                   
    }             
        
    // create "dump" file (Outputs a file named <outputfile>.dump that contains a list of all subgraphs found in the input graph as a list (including their adjacency matrix and the participating vertices)
    bool gen_dumpfile = atoi(argv[17]); 
           
	// Some vars which are necessary, but not params
	uint64 SMPLS = 1000000;
	long total_tries = 0, total_success = 0;
	int perc_number = 1, netevent_number = 1; // How often are percentage events send to the main program
	std::ostringstream info; // String which is sent to the main program

	// algorithm options
	short G_N = atoi(argv[1]); // subgraph size
	long TREESMPLS = atol(argv[2]);
	if (TREESMPLS>10000000 || TREESMPLS<0){
	    errors_array[0] = new string("Number of samples should be between 0 and 1e7\n");
	}
	bool fullenumeration = atoi(argv[3]);
	
	// input file
	string inputfile = argv[4];
	bool directed = atoi(argv[5]);
	bool has_vertex_colors = atoi(argv[6]);
	bool has_edge_colors = atoi(argv[7]);
    // random networks
	short random_type;
	// random type should be 'local const' if graph is undir
	if (directed){
	    random_type = atoi(argv[8]);
	}
	else{
	    random_type = 2;
	}
	// to ensure that vertex/edge color does not matter, if there is none
	bool vertex_color_matters = (atoi(argv[9]) && has_vertex_colors);
	bool edge_color_matters = (atoi(argv[10]) && has_edge_colors);
	bool reest_size = atoi(argv[11]);
	long num_r_nets = atol(argv[12]);
	if (num_r_nets<1){
	    errors_array[1] = new string("Number of random networks should be a positive integer\n");
    }
	long num_exchanges = atol(argv[13]);
	if (num_exchanges<1){
	    errors_array[2] = new string("Number of exchanges per edge should be a positive integer\n");
	}
	long num_tries = atol(argv[14]);
	if (num_tries<1){
	    errors_array[3] = new string("Number of tries per edge should be a positive integer\n");
	}
	
    // output file
	string outputfile = argv[15];
	bool text_output = atoi(argv[16]);
	
    // sampling probabilities (if NOT full enumeration)
	double prob[G_N];
    if (!fullenumeration){
	    for (int i=0; i<G_N; i++){
		    prob[i] = atof(argv[i+18]);
		    if (prob[i]<0 || prob[i]>1){
			   errors_array[4] = new string("Sampling probabilities should be between 0 and 1\n");	
		    }
	    }
	}
        
    // open log file
    std::ofstream logfile((outputfile+".log").c_str(), ios::out);
	logfile<<"FANMOD log file"<<endl;
	logfile<<"Run began at "<<asctime(timeinfo)<<endl;
	
	// write all errors to the log file and exit (if there are any errors)
	bool errors_exist = false;
	for (int i=0 ; i< 100; i++){
        if (errors_array[i]!=NULL) {
           if (errors_exist) {
              logfile<<*errors_array[i]<<endl;            
           }                  
           else {
              logfile<<"Errors in command line arguments:"<<endl;  
              logfile<<*errors_array[i]<<endl;  
              errors_exist = true;
           }  
           delete errors_array[i]; 
        }
     }

     logfile.close();
     if (errors_exist) exit(-1);

    logfile.open((outputfile+".log").c_str(),ios::app);
    logfile<<endl<<"Runtime messages:"<<endl<<endl;
    logfile.close();
    
    // got all params, starting calculations
    
	// Initialize Graph structure from file
	maingraph maing;
	if (!read_graph(maing, inputfile, directed, has_vertex_colors, 
		has_edge_colors, G_N, gen_dumpfile, outputfile)) {           // If reading of the graph failed
        logfile.open((outputfile+".log").c_str(),ios::app);
		logfile << "ERROR: Reading of the graph failed" << endl;
		logfile.close();
        exit(-1);
	}

	// Fill the arrays which store neighbour information
	build_graph(maing);
	
	// v_extension is used in the sampling and est_tree_size-functions
	register long* v_extension = new long[maing.maxnumneighbours*G_N];

	// Write general info into the outputfile
	std::ofstream outfile (outputfile.c_str());
	outfile << "FANMOD 1.1 subgraph "
		<< (fullenumeration ? "enumeration" : "sampling") << endl;
	outfile << "-------------------------------\n" << endl;
	outfile << "Network name: " << inputfile << endl;
	outfile << "Network type: " << (maing.directed? "Directed" : "Undirected") << endl;
	outfile << "Number of nodes: " << maing.n << endl
		<< "Number of edges: " << maing.m << endl;
	outfile << "Number of single edges: " << maing.num_dir_edges << endl
		<< "Number of mutual edges: " << maing.num_undir_edges << endl;
	if (has_vertex_colors) outfile << "Number of vertex colors: " << maing.num_vertex_colors << endl;    
	if (has_edge_colors) outfile << "Number of edge colors: " << maing.num_edge_colors << endl;
	outfile << endl << "Algorithm: " << (fullenumeration ? "enumeration" : "sampling") << endl;
	if (!fullenumeration) {
		outfile << "Sampling parameters = {";
		for (int i = 0; i!= G_N; ++i)
			outfile << " " << prob[i];
		outfile << " }" << endl;
	}
	outfile << "Subgraph size: " << G_N << endl << endl;
	string type_str = "";
	if (num_r_nets == 0)
	{
		outfile << "No random graphs were generated " << endl;
	} else {
		switch (random_type) {
			case 0: type_str = "no respect to bidirectional edges,"; break;
			case 1: type_str = "globally constant number of bidirectional edges,";
				break;
			case 2: type_str = "locally constant number of bidirectional edges,";
				break;
		}
		if (has_edge_colors && edge_color_matters)
			type_str += "\n   exchanging only edges with the same color,";
		if (has_vertex_colors && vertex_color_matters)
			type_str += "\n   exchanging only edges with the same vertex colors,";
		outfile << "Generated " << num_r_nets << " random networks" << endl
			<< "   with " << type_str << endl
			<< "   " << num_exchanges << " exchanges per edge and "
			<<  num_tries << " tries per edge." << endl << endl;
	}
	outfile.close(); // Close the file so it is not open during main algorithm.

	// Initializers for main algorithm 
	randlib::rand rand(time(NULL));

	// Estimate number of subgraphs by random tree sampling
	uint64 numtrees = 0x28F5C28F5C28F5CULL;

	if (TREESMPLS > 0) {
		numtrees = est_tree_size(maing, v_extension, TREESMPLS, G_N, rand) / TREESMPLS * maing.n;
	    // print to log file
        logfile.open((outputfile+".log").c_str(),ios::app);
		logfile << "Approximate number of subgraphs: " << numtrees
			<< " (based on " << TREESMPLS << " samples)\n"; 
	
		
		if (!fullenumeration){
			SMPLS = calc_expected_samples(double(numtrees), G_N, prob);
			logfile << "Expected number of sampled subgraphs: " << SMPLS << endl;
		}
		logfile.close();
	}

	// For the Status bar: pass the number of expected subgraphs to sampling
	uint64 equiv100p = fullenumeration ? numtrees : SMPLS;

	// Start the sampling / enumeration of subgraphs

	// Main enumeration / sampling loop 
	int total_num_nets = num_r_nets + 1; // We sample num_r_nets and the original graph
	// In this hashmap, the intermediate result is stored before entered into the result - hashmap.
	hash_map < graphcode64, uint64 > inter_result;
	// In this hashmap, the results of the sampling and the randomization are stored.
	hash_map < graphcode64, uint64* > result_graphs;
	uint64 *count_subgr = new uint64[total_num_nets];
	uint64 total_subgr = 0;
	double sampling_time = 0.0, random_time = 0.0;
	uint64 *current_array; // Abbrevation for the array currently updated

	// Sample the original graph and random graphs, if necessary.
	EdgeContainer EC;    

	for (int nets_ctr = 0; nets_ctr < total_num_nets; ++nets_ctr){

        // Reestimates tree-size if the user wishes.
		// and it is not the original network being sampled
		if (reest_size && nets_ctr != 0 && TREESMPLS > 0) {
			numtrees = est_tree_size(maing, v_extension, TREESMPLS, G_N, rand) / TREESMPLS * maing.n;
			if (!fullenumeration) SMPLS = calc_expected_samples(double(numtrees), G_N, prob);
			equiv100p = fullenumeration ? numtrees : SMPLS;
		}

		//Create vector where subgraphs can be dumped to
		vector<subgraph> subgraphdump;
		const bool subgraph_dumpfile = (gen_dumpfile && (nets_ctr == 0));

		// Sample the graph
		// sampling is called with its options,
		// true to show the status-bar and inter_result and count_subgr as result-parameters
		// The frame and this thread are passed to send events and "TestDestroy()"
		sampling_time += sampling(maing, v_extension, G_N, fullenumeration, prob, 
			equiv100p, perc_number, inter_result, 
			count_subgr[nets_ctr], rand,
			subgraph_dumpfile, subgraphdump);
		total_subgr += count_subgr[nets_ctr];

		//Write the subgraph dump
		if (subgraph_dumpfile) {
			//TO DO  
			//string outputfile 
			std::ofstream dumpfile ((outputfile+".dump").c_str(), ios::out);
			dumpfile << "Number of subgraphs: " << subgraphdump.size() << endl;
			dumpfile << "Format: adjacency matrix, <participating vertices>" << endl;
			graph64 g;
			init_graph(g,G_N,maing.num_vertex_colors,maing.num_edge_colors,maing.directed); 
			for (vector<subgraph>::const_iterator iter = subgraphdump.begin();
				iter !=subgraphdump.end(); ++iter) {
					readHashCode(g, iter->gc);
					for (short i = 0; i != G_N; ++i) {
						for (int j = 0; j != G_N; ++j) {
							dumpfile << num_to_lett[get_element(g,i,j)];
						}
					}
					for (short i = 0; i != G_N; ++i) {
						dumpfile << "," << getvertex(*iter,i);
					}
					dumpfile << endl;
				}
				dumpfile.close();
				subgraphdump.clear();
		}

		if (nets_ctr == 0){
	        logfile.open((outputfile+".log").c_str(),ios::app);
			logfile << count_subgr[0] << " subgraphs were "
				<< (fullenumeration ? "enumerated" : "sampled")
				<< " in the original network." << endl << endl;
		    logfile.close();
		}
		
		// Write results into result_graphs-hashmap
		for (hash_map < graphcode64, uint64 >::const_iterator iter = inter_result.begin();
			iter !=inter_result.end(); ++iter){
				if (result_graphs.find(iter->first) == result_graphs.end()){
					// create the array at the graph's hashmap position.
					current_array = (result_graphs[iter->first] = new uint64[total_num_nets]);
					// Initialize the array
					for (int i=0; i < total_num_nets; ++i)
						current_array[i] = 0;
				} else {
					// set current_array as the array at the graphs map position
					current_array = result_graphs[iter->first];
				}
				// Write the number of subgraphs in the array
				current_array[nets_ctr] = iter->second;
			} // end for
			inter_result.clear();       // Empty the intermediate result hashmap.

			// Randomize the graph
			if (num_r_nets != 0) {

				if (nets_ctr == 0) {
					//build EC, needs only be done once because sampling
					//does not change the edges and hence the content remains
					//correct
					short color_u = 1, color_v = 1, color_uv, color_vu;
					for (hash_map < edge, edgetype >::const_iterator iter = maing.edges.begin();
						iter != maing.edges.end(); ++iter) {
							const edge e = iter->first;
							if (vertex_color_matters) {
								color_u = maing.vertex_colors[edge_get_u(e)];
								color_v = maing.vertex_colors[edge_get_v(e)];
							}
							const edgetype et = iter->second;
							if (edge_color_matters) {
								color_uv = get_color_u_v(et);
								color_vu = get_color_v_u(et);              
							} else {
								color_uv = et & DIR_U_T_V;
								color_vu = (et & DIR_V_T_U) >> 1;              
							}
							//Put edge into EC. For NO_REGARD, the bidirectional
							//edge is split up, additionally, the edge code shows
							//an edge's direction
							switch (random_type) {
					 case NO_REGARD    : if ((et & DIR_U_T_V) == DIR_U_T_V)
											 EC.put(new_edge(edge_get_u(e),edge_get_v(e)), getBagID(color_uv, 0 ,color_u,color_v) );
						 if ((et & DIR_V_T_U) == DIR_V_T_U)
							 EC.put(new_edge(edge_get_v(e),edge_get_u(e)), getBagID( 0 , color_vu ,color_u,color_v) );
						 break;       
					 case GLOBAL_CONST : EC.loggedPut(e, getBagID(color_uv,color_vu,color_u,color_v) );
						 break;
					 case LOCAL_CONST  : EC.put(e, getBagID(color_uv,color_vu,color_u,color_v) );
						 break;
							}
						}

				}
				random_time += randomize_graph(maing, random_type, num_exchanges, num_tries, 
					vertex_color_matters, edge_color_matters, 
					EC, total_tries, total_success, rand);
			}  
			
      // time counter - calculates time after each sampled random network
      time (&newtime);
      timeinfo = localtime (&newtime);
      double dif = difftime (newtime,oldtime);
      oldtime = newtime;
      logfile.open((outputfile+".log").c_str(),ios::app);
      if (nets_ctr==0){
         logfile<<"Original network: finished processing at "<<asctime(timeinfo)<<" after "<<(dif/60)<<" minutes"<<endl;              
      }
      else{
         logfile<<"Random network #"<<nets_ctr<<" finished processing at "<<asctime(timeinfo)<<" after "<<(dif/60)<<" minutes"<<endl;
      }
      logfile.close();                     

	} // end for

	// Free memory
	maing.edges.clear();
	for (vertex i = 0; i != maing.n; ++i) {
		delete[] maing.neighbours[i] ;
	}
	delete[] maing.num_neighbours;
	delete[] maing.neighbours;
	//delete[] maing.num_larger_neighbours;
	delete[] maing.v_util;
	delete[] v_extension;

	outfile.open(outputfile.c_str(), std::ofstream::out | std::ofstream::app);
	outfile << count_subgr[0] << " subgraphs were "
		<< (fullenumeration ? "enumerated" : "sampled")
		<< " in the original network." << endl;
	if (num_r_nets != 0) {
		outfile << (total_subgr-count_subgr[0]) << " subgraphs were "
			<< (fullenumeration ? "enumerated" : "sampled")
			<< " in the random networks." << endl;
		outfile << total_subgr << " subgraphs were "
			<< (fullenumeration ? "enumerated" : "sampled")
			<< " in all networks." << endl << endl;

		outfile << "For the random networks: " << total_tries << " tries were made, "
			<< total_success << " were successful." << endl;
		outfile << "Randomization took " << random_time << " seconds." << endl;
	} else {
		outfile << endl;
	}
	outfile << (fullenumeration ? "Enumeration took " : "Sampling took ")
		<< sampling_time << " seconds.\n\n" << endl;

	// Output the results table
	// if text_output is false, CSV output is given.
	pretty_output(text_output, result_graphs, G_N, maing.num_vertex_colors,
		maing.num_edge_colors, maing.directed, count_subgr, 
		num_r_nets, outfile);
	outfile.close();

	// Free Memory of result storage
	for (hash_map < graphcode64, uint64* >::const_iterator iter = result_graphs.begin();
		iter != result_graphs.end(); ++iter){
			delete[] iter->second;
		}
	result_graphs.clear();
	delete[] count_subgr;

	return 0;
}

// Calculates the expected number of sampled subgraph
uint64 calc_expected_samples(double ntrees, const short & G_N, const double *prob)
{
	for (int i = 0; i != G_N; ++i)
	{
		ntrees *= prob[i];
	}
	return (uint64)(ntrees+0.5);
}
