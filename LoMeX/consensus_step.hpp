#ifndef CONSENSUS_STEP_H
#define CONSENSUS_STEP_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>

using namespace std;


tuple<int, int, int> run_consensus_step_multirhread(std::string work_dir, std::string output_path, int total_length, int fixed_length, double relative_nuc_threshold, int abs_min_nuc, int delete_files, vector<bool> & character_status, int iterations, int iteration, int threads);

tuple<int, int, int, int> run_consensus_step(std::string work_dir, std::string output_path, std::string output_file_name, int total_length, int fixed_length, double relative_nuc_threshold, int abs_min_nuc, int delete_files, vector<bool> & character_status);


#endif