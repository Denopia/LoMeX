#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <map>
#include <random>

using namespace std;


/*
	Write locations into a file
*/
std::string write_locations(map<__uint128_t, vector<int> > &kmer2positions, std::string work_dir, int file_number, int fixed_length);


/*
	Write occurrences into a file

	SUPPORT ONLY FOR k-MERS UP TO k = 255*4 = 1020 

*/
std::string write_occurrences_binary(map<__uint128_t, vector<std::string> > &kmer2occurrences, std::string work_dir, int file_number, int fixed_length, int total_length);


/*
	Format the file number
*/
std::string format_number(int number);


/*
	Maps 4 nucleotide characters into a byte
*/
uint8_t map_4nucs2byte(char nucs[4]);


/*
	Gives a random nucleotide character
*/
char random_nucleotide();


/*
	Map a byte into four nucleotides
*/
std::string map_byte2fournucs(uint8_t byte);

#endif