/*
 * This is a modified version of the original file appearing in Squeakr. 
 * Modification by Miika Leinonen
*/

/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include <fstream>

#include <iostream>
#include <regex>
#include <string>
#include <tuple>

#include "kmer.h"


#include "spdlog/spdlog.h" // BONUS


/*return the integer representation of the base */
char Kmer::map_int(uint8_t base)
{
	switch(base) {
		case DNA_MAP::A: { return 'A'; }
		case DNA_MAP::T: { return 'T'; }
		case DNA_MAP::C: { return 'C'; }
		case DNA_MAP::G: { return 'G'; }
		default:  { return DNA_MAP::G+1; }
	}
}

/*return the integer representation of the base */
uint8_t Kmer::map_base(char base)
{
	switch(base) {
		case 'A': { return DNA_MAP::A; }
		case 'T': { return DNA_MAP::T; }
		case 'C': { return DNA_MAP::C; }
		case 'G': { return DNA_MAP::G; }
		default:  { return DNA_MAP::G+1; }
	}
}

/**
 * Converts a string of "ATCG" to a uint64_t
 * where each character is represented by using only two bits
 */
__int128_t Kmer::str_to_int(std::string str)
{
	__int128_t strint = 0;
	for (auto it = str.begin(); it != str.end(); it++) {
		uint8_t curr = 0;
		switch (*it) {
			case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
		}
		strint = strint | curr;
		strint = strint << 2;
	}
	return strint >> 2;
}

/**
 * Converts a uint64_t to a string of "ACTG"
 * where each character is represented by using only two bits
 */
std::string Kmer::int_to_str(__int128_t kmer, uint64_t kmer_size)
{
	uint8_t base;
	std::string str;
	for (uint32_t i = kmer_size; i > 0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = Kmer::map_int(base);
		str.push_back(chr);
	}
	return str;
}

/* Return the reverse complement of a base */
int Kmer::reverse_complement_base(int x) { return 3 - x; }

/* Calculate the revserse complement of a kmer */
__int128_t Kmer::reverse_complement(__int128_t kmer, uint64_t kmer_size)
{
	__int128_t rc = 0;
	uint8_t base = 0;
	for (uint32_t i = 0; i < kmer_size; i++) {
		base = kmer & 3ULL; // 3 unsigned long long ?
		base = reverse_complement_base(base);
		kmer >>= 2;
		rc |= base;
		rc <<= 2;
	}
	rc >>=2;
	return rc;
}

/* Compare the kmer and its reverse complement and return the result 
 * Return true if the kmer is greater than or equal to its
 * reverse complement. 
 * */
bool Kmer::compare_kmers(__int128_t kmer, __int128_t kmer_rev)
{
	return kmer >= kmer_rev;
}


/* * [For the slow implementations]
 *
 * Take in a sequence of length f+g where
 * f is the number of fixed characters and
 * g is the number of characters falling into gaps
 * AND the gapped k-mer pattern.
 *
 * Returns the concatenation of the fixed characters
 * in the 2-bit representation AND a boolean flag
 * telling if the gapped k-mer is valid (i.e. contains
 * only characters ACTG).
 * */
std::tuple<uint64_t, bool> Kmer::read_gapped_kmer(std::string  & full_read, std::vector<bool> & character_status, int start, int length){
	uint64_t gapped_kmer = 0;

	for(int i = 0; i < length; i++)
	{
		if(character_status[i])
		{
			uint8_t current_character = Kmer::map_base(full_read[start+i]);
			if (current_character > DNA_MAP::G)
			{
				return std::make_tuple(gapped_kmer, false);
			}
			else
			{
				gapped_kmer <<= 2;
				gapped_kmer |= current_character; // bitwise OR operation
				

			}
		}
	}
	//gapped_kmer = gapped_kmer >> 2;
	return std::make_tuple(gapped_kmer, true);
}

/*
* Take in a string representing a gapped kmer structure
* and produce the corresponding boolean vector. 
*
*
*/
std::tuple<std::vector<bool>, int, int> Kmer::get_gapped_kmer_shape(std::string kmer_shape)
{
	std::regex rgx("-+");
    std::sregex_token_iterator iter(kmer_shape.begin(), kmer_shape.end(), rgx, -1);
    std::sregex_token_iterator last{};
    std::vector<int> segments;
    int number_of_segments = 0;
    int total_length = 0;
    int fixed_length = 0;

    for ( ; iter != last; ++iter)
    {
    	number_of_segments += 1;
    	int segment_length = stoi(*iter);
    	segments.insert(segments.end(), 1, segment_length);
    }

	std::vector<bool> character_status;
	bool current_status = true;
	for (int i = 0; i < number_of_segments; i++)
	{
		int segment_size = segments[i];
		character_status.insert(character_status.end(), segment_size, current_status);
		total_length += segment_size;
		if(current_status)
		{
			fixed_length += segment_size;	
		}
		current_status = !current_status;
	}


	return std::make_tuple(character_status, total_length, fixed_length);
}


/*
* Read kmers from a file 
* Requires file name, kmer size and reference to a set of kmers
* 
*
**/
//void Kmer::parse_kmers(const char *filename, uint64_t kmer_size, std::unordered_set<uint64_t>& kmerset, std::string kmer_shape) {
void Kmer::parse_kmers(const char *filename, uint64_t kmer_size, std::unordered_set<uint64_t>& kmerset, spdlog::logger* console) {
	
	console->info("Parsing k-mers.");

	std::ifstream ipfile(filename);
	std::string read;

	std::string kmer_shape = "6-4-6-4-6-4-6-4-6"; // 6s are fixed characters and 4s are don't care

	std::vector<bool> character_status;
	uint total_length, fixed_length;

	std::tie(character_status, total_length, fixed_length) = Kmer::get_gapped_kmer_shape(kmer_shape);

	// do while there are reads available
	while (ipfile >> read) {
		uint64_t item = 0;
		uint64_t reverse_gap_kmer = 0;


		for (uint i = 0; i < read.length()-total_length; i++){
			uint64_t gap_kmer;
			bool valid;
			std::tie(gap_kmer, valid) = Kmer::read_gapped_kmer(read, character_status, i, total_length);
			if(valid)
			{
				reverse_gap_kmer = Kmer::reverse_complement(gap_kmer, kmer_size);
				if(Kmer::compare_kmers(gap_kmer, reverse_gap_kmer))
				{
					item = gap_kmer;
				}
				else
				{
					item = reverse_gap_kmer;
				}
				kmerset.insert(item);
			}
		}
	}
}
