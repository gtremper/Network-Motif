#include "output.hpp"
using std::setw;

double calc_mean(double* array, int num){

    if (num == 0) return NaN;
    double ret = 0;
    for (int i=1; i <= num; ++i){
        ret += array[i];
    }
    return ret / num;
}

double calc_deviation(double* array, int num, double mean){

    if (num <= 1) return NaN;
    double ret = 0;
    for (int i=1; i <= num; ++i){
        ret += std::pow(array[i] - mean, 2);
    }
    ret /= num-1;
    return std::sqrt(ret);
}

// Converts an int into string
string int_to_str(const int & i)
{
   std::stringstream s;
   s<<i;
   return s.str();
}

// Get the number of digits an uint64 has.
int num_digits(uint64 i)
{
    int ret=1;
    while (i > 9) {
        i /= 10;
        ++ret;
    }
    return ret;
}

char numberstoletters[] = {'0','1','2','3','4','5','6','7',
                           '8','9','A','B','C','D','E','F', 
                           '#','#','#','#','#','#','#','#',}; //last line for safety only

void
pretty_output(const bool textout, hash_map < graphcode64, uint64* > & res_graphs, 
              short G_N, unsigned short num_v_colors, unsigned short num_e_colors,
              bool directed, uint64* count_subgr, int num_r_nets, std::ofstream & outfile)
{
    vector < subgr_res > result_vec(res_graphs.size());
    subgr_res gr;
    int total_num_nets = num_r_nets + 1; // Random nets plus the original
    double *concentration = new double[total_num_nets];
    uint64 max_ID = 0; //For the result table: max_ID has to fit in first column

    int idx = 0;
    for (hash_map < graphcode64, uint64* >::const_iterator iter =
	 res_graphs.begin(); iter != res_graphs.end(); ++iter) {
		gr.id = iter->first;
        // Get the largest ID
        if (gr.id > max_ID) max_ID = gr.id;

        // Calculate frequency in the original graph
		gr.freq = ((double) (iter->second)[0] / (double) count_subgr[0]);
        if (gr.freq > 0) { // NEW: Only if original freq > 0, write the graph to results!
            concentration[0] = gr.freq;
            if (num_r_nets > 0) { // at least one random network was sampled
                gr.p_value = 0;
                for (int i=1; i < total_num_nets; ++i){
                    // Build the concentration array for this subgraph
                    concentration[i] = ((double) (iter->second)[i] / (double) count_subgr[i]);
                    // Calculate the p-value of the subgraph
                    if (concentration[i] > gr.freq) ++gr.p_value;
                }
                gr.p_value /= double(num_r_nets);
                // Calculate the Z-Score via mean and standard deviation
                gr.rand_mean = calc_mean(concentration, num_r_nets);
                gr.rand_sd = calc_deviation(concentration, num_r_nets, gr.rand_mean);
                if (gr.rand_sd == 0) gr.z_score= NaN;
                else
                   gr.z_score = (gr.freq - gr.rand_mean) / gr.rand_sd; 
            }
    		result_vec[idx++] = gr;
        }
    }
    sort(result_vec.begin(), result_vec.end(), compare);
    //const uint64 bitmask = 1ULL << (G_N * G_N - 1);

    // Width constants
    const int id_width = num_digits(max_ID), adj_width = G_N+1, freq_w = 10,
              mfreq_w = 13, sd_width = 14, z_width = 11, p_width = 9;
    outfile.precision(5);
    outfile.flags(std::ofstream::right);

    outfile << "Result overview:\n" << endl;
    // First header line of the result table
    if (textout)
        outfile << setw(id_width) << "ID" << setw(adj_width) << "Adj"
                << setw(freq_w+2) << "Frequency";
    else
        outfile << "ID,Adj-Matrix,Frequency";
    if (num_r_nets > 0)
        if (textout)
            outfile << setw(mfreq_w+1) << "Mean-Freq" << setw(sd_width) << "Standard-Dev"
                    << setw(z_width) << "Z-Score" << setw(p_width) << "p-Value";
        else
            outfile << ",Mean-Freq,Standard-Dev,Z-Score,p-Value";
    outfile << endl;

    // Second header line of the result table
    if (textout)
        outfile << setw(id_width) << ' ' << setw(adj_width) << ' '
                << setw(freq_w+2) << "[Original]";
    else
        outfile << ",,[Original]";
        
    if (num_r_nets > 0)
        if (textout)
            outfile << setw(mfreq_w+1) << "[Random]" << setw(sd_width) << "[Random]";
        else
            outfile << ",[Random],[Random]";
    outfile << endl << endl;

    graph64 g;
    init_graph(g,G_N,num_v_colors,num_e_colors,directed);
     
    // Datalines of the result table
    for (vector < subgr_res >::const_iterator iter = result_vec.begin();
         iter != result_vec.end(); ++iter) {
        if (iter->freq > 0) { // NEW: Output only if freq > 0
            readHashCode(g, iter->id );
            if (textout)
      		    outfile << setw(id_width) << getGraphID(g,iter->id) << ' ';
            else
                outfile << getGraphID(g,iter->id) << ',';
                
            // Output the first line of the adj-matrix
            for (int j = 0; j != G_N; ++j) {
                outfile << numberstoletters[get_element(g,0,j)];
            }
            if (textout) { // write the rest of the line, either as text or comma separated
                outfile << ' ';
        		outfile	<< setw(freq_w) << iter->freq*100 << '%';
                if (num_r_nets > 0){ // Only if random nets have been sampled
                    outfile << setw(mfreq_w) << iter->rand_mean*100 << '%'
                            << setw(sd_width);
                    if (std::isinf(iter->rand_sd) || std::isnan(iter->rand_sd))
                        outfile << "undefined";
                    else outfile << iter->rand_sd;
                    outfile << setw(z_width);
                    if (std::isinf(iter->z_score) || std::isnan(iter->z_score))
                        outfile << "undefined";
                    else outfile << iter->z_score;
                    outfile << setw(p_width) << iter->p_value;
                }
                outfile << endl;
            } else {
                outfile << ',' << iter->freq*100 << '%';
                if (num_r_nets > 0){ // Only if random nets have been sampled
                    outfile << ',' << iter->rand_mean*100 << '%' << ',';
                    if (std::isnan(iter->rand_sd))
                        outfile << "undefined,";
                    else outfile << iter->rand_sd << ',';
                    if (std::isnan(iter->z_score))
                        outfile << "undefined";
                    else outfile << iter->z_score;
                    outfile << ',' << iter->p_value;
                }
                outfile << endl;
            }            
            // Output the remaining lines of the adj-matrix
    		for (int i = 1; i != G_N; ++i) {
                if (textout)
                    outfile << setw(id_width+1) << ' ';
                else
                    outfile << ',';
    			for (int j = 0; j != G_N; ++j) { // give one line of the matrix
        			outfile << numberstoletters[get_element(g,i,j)];
    			}
    			outfile << endl;
    		}
    		outfile << endl;
        } // end if freq > 0
    } // end for iterator

    delete[] concentration;
}
