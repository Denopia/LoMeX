#ifndef FUN_CONSENSUS_H
#define FUN_CONSENSUS_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>

using namespace std;



/*
	Determine consensus k-mer in the unambiguous consensus case
*/
int determine_unambiguous_consensus(std::ofstream & output_file, int pattern_length, vector<vector<int> > & consensus_nucleotides, vector<vector<int> > & kmer_nucleotide_occurrences, int unambiguous_counter);


/*
	Determine consensus k-mer in the simple ambiguous consensus case
*/
int determine_simple_ambiguous_consensus(std::ofstream & output_file, int pattern_length, vector<vector<int> > & consensus_nucleotides, vector<vector<int> > & kmer_nucleotide_occurrences, int ambiguous_positions[], int simple_ambiguous_counter);


/*
	Determine consensus k-mer in the complex ambiguous consensus case
*/
int determine_complex_ambiguous_consensus(std::ofstream & output_file, int pattern_length, vector<vector<int> > & consensus_nucleotides, vector<vector<int> > & kmer_nucleotide_occurrences, int occurrences,  int ambiguous_positions[],
										  vector<vector<char> > & ambiguous_patterns, int ambiguous_count, int complex_ambiguous_counter);


/*
	Extract ambiguous character position sequences
*/
vector<vector<char> > extract_ambiguous_position_matrix(int pattocc_count, int ambipos_count, vector<std::string> & kmer_occurrences, 
										 			   int ambiguous_positions[], int total_length);

#endif