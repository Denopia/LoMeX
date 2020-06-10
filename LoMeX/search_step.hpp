#ifndef SEARCH_STEP_H
#define SEARCH_STEP_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>

using namespace std;



void run_search_step(std::string work_dir, std::string kmers_path, std::string reads_path, int buffer_size, int total_length, int fixed_length, vector<bool> & character_status, int kmer_min, int curite, int maxite);


#endif