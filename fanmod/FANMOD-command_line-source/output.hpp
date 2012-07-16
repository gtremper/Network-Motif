#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "graph64.hpp"
#include "maingraph.hpp"
#include "hashmap.h"


struct subgr_res {
    graphcode64 id;
    double freq;
    double rand_mean;
    double rand_sd;
    double p_value;
    double z_score;
};

const unsigned long nan_field[2]={0xffffffff, 0x7fffffff};
const double NaN = *( double* )nan_field;

inline bool compare(const subgr_res & a, const subgr_res & b)
{
    return a.freq > b.freq;
}

// Converts an int into string
std::string int_to_str(const int & i);

void
pretty_output(const bool textout, hash_map < graphcode64, uint64* > & res_graphs,
              short G_N, unsigned short num_v_colors, unsigned short num_e_colors, 
              bool directed, uint64* count_subgr, int num_r_nets, std::ofstream & outfile);

#endif
