/*

	Consensus functionality

*/

#include "fun_consensus.hpp"
#include "fun_kmers.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>

//namespace po = boost::program_options;

using namespace std;


/*
	Array to map nucleotide integers to characters 
*/
char NUC_ARRAY[4] = {'C', 'A', 'T', 'G'};



int determine_unambiguous_consensus(std::ofstream& output_file, int pattern_length, vector<vector<int> > &consensus_nucleotides, int unambiguous_counter)
{
	std::string ua_kmer = "";
	// Loop through all pattern positions
	for (int i = 0; i < pattern_length; i++)
	{
		// Loop through the four possible nucleotides (one of which is marked with 1)
		for (int j = 0; j < 4; j++)
		{
			if (consensus_nucleotides[j][i] == 1)
			{
				//std::cout << NUC_ARRAY[j];
				//output_file << NUC_ARRAY[j];
				ua_kmer = ua_kmer + NUC_ARRAY[j];
				break;
			}
		}
	}
	//std::cout << std::endl;

	//std::string ua_kmer_rev = reverse_complement_seqstr(ua_kmer);
	//if (compare_seqs_string(ua_kmer, ua_kmer_rev)){output_file << ua_kmer;}
	//else{output_file << ua_kmer_rev;}

	output_file << ua_kmer;
	output_file << " U";
	output_file << '\n';

	return unambiguous_counter+1;
}

// Called multiple times, make it work as a single-callable function?
tuple<int, bool> determine_simple_ambiguous_consensus(std::ofstream& output_file, int pattern_length, vector<vector<int> > &consensus_nucleotides, int simple_ambiguous_counter)
{
	std::string sa_kmer = "";
	bool rerun_required = false; // Flag to show if another consensus k-mer can be found after this one
	// Loop through all pattern positions
	for (int i = 0; i < pattern_length; i++)
			{
				int printed_nuc = -1; // Store the printed nuclotide integer here
				// Loop through the four possible nucleotides (one of which is marked with 1)
				for (int j = 0; j < 4; j++)
				{
					if (consensus_nucleotides[j][i] == 1)
					{
						// If a nucleotide has already been printed for this position
						if (printed_nuc > -1)
						{
							// Rerun is needed
							rerun_required = true;
							// The printed nucleotide in an ambiguous position can be ignored in the future
							consensus_nucleotides[printed_nuc][i] = 0;
						}
						// If a nucleotide has not been printed already for this position
						else
						{
							//std::cout << NUC_ARRAY[j];
							//output_file << NUC_ARRAY[j];
							sa_kmer = sa_kmer + NUC_ARRAY[j];
							printed_nuc = j;
						}
					}
				}
			}
			//std::cout << std::endl;
			
			//std::string sa_kmer_rev = reverse_complement_seqstr(sa_kmer);
			//if (compare_seqs_string(sa_kmer, sa_kmer_rev)){output_file << sa_kmer;}
			//else{output_file << sa_kmer_rev;}
			
			output_file << sa_kmer;
			output_file << " S";
			output_file << '\n';

	return make_tuple(simple_ambiguous_counter+1, rerun_required);
}


int determine_complex_ambiguous_consensus(std::ofstream& output_file, int pattern_length, vector<vector<int> > &consensus_nucleotides, int occurrences,  int ambiguous_positions[],
										  vector<vector<char> > &ambiguous_patterns, int ambiguous_count, int complex_ambiguous_counter)
{
	std::string ca_kmer_rev;
	std::string ca_kmer;
	char fill_char;
	uint8_t fill_char_int;
	int ambipos;
	int ambiguous_position_matches;
	bool is_a_match;
	bool usable_kmer;

	// Array to store if an occurrence has already been used to fill in the ambiguous character positions 
	int used[occurrences];
	for (int i = 0; i < occurrences; i++)
	{
		used[i] = 0;
	}

	// Go through all k-mer matching occurrences and determine if it appears at least twice, and print the resulting k-mer
	for (int i = 0; i < occurrences; i++)
	{
		ambiguous_position_matches = 0; // How many times the ambiguous character sequence appears
		
		// If occurrence already used, skip it
		if (used[i] == 1){continue;}

		// Now go through all other later k-mer matching occurrences
		for (int j = i+1; j < occurrences; j++)
		{
			// If occurrence already used, skip it
			if (used[j] == 1){continue;}

			is_a_match = true; // Flag to tell if occurrences are the same
			// Loop through all ambiguous positions
			for (int k = 0; k < ambiguous_count; k++)
			{
				// If at any position the characters are different, sequences do not match
				if ((ambiguous_patterns[i][k] != ambiguous_patterns[j][k]))
				{
					is_a_match = false;
					break;
				}
			}
			// If match is found, mark sequence as used
			if (is_a_match)
			{
				ambiguous_position_matches += 1;
				used[j] = 1;
			}
		}

		// Now, if at least one another occurrence has the same ambiguous characters, we print the resulting long k-mer
		if (ambiguous_position_matches > 0)
		{
			ca_kmer = "";
			usable_kmer = true;
			ambipos = 0; // position in the ambigous characters sequence

			// For every position in the whole long k-mer
			for (int patt_pos = 0; patt_pos < pattern_length; patt_pos++)
			{
				// If position was determined ambiguous, print the character from the ambiguous character sequence
				if (ambiguous_positions[patt_pos] == 1)
				{

					fill_char = ambiguous_patterns[i][ambipos];
					fill_char_int = map_nuc2int(fill_char);
					if (consensus_nucleotides[fill_char_int][patt_pos] == 0){usable_kmer = false;}
					//std::cout << ambiguous_patterns[i][ambipos];
					//output_file << ambiguous_patterns[i][ambipos];
					ca_kmer = ca_kmer + fill_char;
					ambipos += 1;
					continue;
				}
				
				// If position is unambiguous, print the unambiguous character from consensus nucleotide matrix
				for (int nuc = 0; nuc < 4; nuc++)
				{
					if (consensus_nucleotides[nuc][patt_pos] == 1)
					{
						//std::cout << NUC_ARRAY[nuc];
						//output_file << NUC_ARRAY[nuc];
						ca_kmer = ca_kmer + NUC_ARRAY[nuc];
					}
				}
			}

			if(usable_kmer)
			{
				//ca_kmer_rev = reverse_complement_seqstr(ca_kmer);
				//if (compare_seqs_string(ca_kmer, ca_kmer_rev)){output_file << ca_kmer;}
				//else{output_file << ca_kmer_rev;}
				output_file << ca_kmer;
				output_file << " C";
				output_file << '\n';
				complex_ambiguous_counter += 1;
			}
		}
	}
	return complex_ambiguous_counter;
}

vector<vector<char> > extract_ambiguous_position_matrix(int pattocc_count, int ambipos_count, vector<std::string> &kmer_occurrences, 
										 				int ambiguous_positions[], int total_length)
{
	// Initialize ambiguous character matrix
	vector<char> ambi_vector_init(ambipos_count, 'L');
	vector<vector<char> > ambiguous_patterns(pattocc_count, ambi_vector_init);
	int ambipos = 0;

	for (int api = 0; api < total_length; api++)
	{
		if (ambiguous_positions[api] == 0){continue;}

		for (int koi = 0; koi < pattocc_count; koi++)
		{
			ambiguous_patterns[koi][ambipos] = kmer_occurrences[koi][api];
		}
		ambipos+=1;
	}
	return ambiguous_patterns;
}
