/*

	File writing functionality

*/

#include "file_writer.hpp"
#include "fun_kmers.hpp"
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <map>
#include <random>
#include <math.h>

using namespace std;


std::string write_locations(map<__uint128_t, vector<int> > &kmer2positions, std::string work_dir, int file_number, int fixed_length)
{
	std::string file_path = work_dir + "/spaced_kmer_locations_" + format_number(file_number);

	std::ofstream location_file(file_path);

	__uint128_t p_spaced_kmer;
	std::vector<int> p_kmer_occurrences;
	int kmer_i = 0;
	int kmer_read, kmer_pos, kmer_reversed;
	std::string current_spaced_kmer;

	for (const auto &pair : kmer2positions)
	{
		p_spaced_kmer = pair.first;
		p_kmer_occurrences = pair.second;

		current_spaced_kmer = map_int2str(p_spaced_kmer, fixed_length); // get new spaced k-mer
		
		location_file << current_spaced_kmer;

		// Start going through the current spaced k-mer occurrences in the data and fill the nucleotide count matrix
		while (kmer_i < p_kmer_occurrences.size())
		{
			kmer_read = p_kmer_occurrences[kmer_i]; // next occurrence read
			kmer_i += 1;
			kmer_pos = p_kmer_occurrences[kmer_i]; // next occurrence read position
			kmer_i += 1;
			kmer_reversed = p_kmer_occurrences[kmer_i]; // next occurrence reversed flag
			kmer_i += 1;

			location_file << "\t" << std::to_string(kmer_read) << "\t" << std::to_string(kmer_pos) << "\t" << std::to_string(kmer_reversed);
		}

	location_file << std::endl;
	kmer_i = 0;
	}

	location_file.close(); // close file
	location_file.clear(); // clear flags
	return file_path;	
}


std::string write_occurrences_binary(map<__uint128_t, vector<std::string> > &kmer2occurrences, std::string work_dir, int file_number, int fixed_length, int total_length)
{
	std::string file_path = work_dir + "/spaced_kmer_occurrences_" + format_number(file_number) + ".bin";	// Determine file path
	uint8_t bytebuffer[1];
	std::ofstream location_file(file_path, ios::out | ios::binary);										// Open the file
	float tl = static_cast<float>(total_length);
	uint8_t single_occurrence_bytes = static_cast<uint8_t>(ceil(tl / 4.0));							// How many bytes a single occurrence takes
	
	bytebuffer[0] = single_occurrence_bytes;
	location_file.write((char *) & bytebuffer, 1);														// Write the number of needed bytes for a single occurrence at the beginning of the file

	__uint128_t spaced_kmer;
	//std::vector<std::string> kmer_occurrences;

	char fournucset[4];
	__uint128_t last_kmer = 0;

	for (const auto &pair : kmer2occurrences)		// Go through every spaced k-mer and its occurrences
	{
		spaced_kmer = pair.first;

		//std::cout << "Current spaced k-mer: " << map_int2str(spaced_kmer, fixed_length) << std::endl;

		if (last_kmer > spaced_kmer)
		{
			std::cout << " -------------------------" << std::endl;
			std::cout << "| HORRIBLE ORDERING ERROR |" << std::endl;
			std::cout << " -------------------------" << std::endl;
		}
		last_kmer = spaced_kmer;

		//kmer_occurrences = pair.second;

		//uint8_t number_of_occurrences = kmer_occurrences.size();											// How many occurrences the current k-mer has
		uint8_t number_of_occurrences = kmer2occurrences[spaced_kmer].size();

		bytebuffer[0] = number_of_occurrences;
		location_file.write((char *) & bytebuffer, 1);														// First, write down the number of occurrences 

		for (int i = 0; i < kmer2occurrences[spaced_kmer].size(); i++)	// In this loop we write every occurrence as a set of bytes (number of bytes indicated by the first byte of the file)
		{
			std::string current_occurrence = kmer2occurrences[spaced_kmer][i];

			//std::cout << "Current occurrence: " << current_occurrence << std::endl;
			
			int nucs_in_set = 0;
			for (char & c : current_occurrence)		// In this loop we write a single occurrence byte by byte
			{
				fournucset[nucs_in_set] = c;
				if (map_nuc2int(c) > 3)
				{
					fournucset[nucs_in_set] = random_nucleotide();
					std::cout << "Randomized" << std::endl;
				}

				nucs_in_set += 1;

				if (nucs_in_set == 4)
				{
					uint8_t nucsetbyte = map_4nucs2byte(fournucset);

					//std::cout << "Nucleotides ";
					//for (int hitsi=0; hitsi<4; hitsi++){std::cout << fournucset[hitsi];}
					//std::cout << " map to byte integer " << (int)nucsetbyte << std::endl;

					bytebuffer[0] = nucsetbyte;
					location_file.write((char *) & bytebuffer, 1);
					nucs_in_set = 0;
				}
			}
			while (nucs_in_set != 0)				// If the we started a byte but did not write it down yet, fill it with garbage and write it
			{
				fournucset[nucs_in_set] = 'B';
				nucs_in_set += 1;
				if (nucs_in_set == 4)
				{
					uint8_t nucsetbyte = map_4nucs2byte(fournucset);

					//std::cout << "Nucleotides ";
					//for (int hitsi=0; hitsi<4; hitsi++){std::cout << fournucset[hitsi];}
					//std::cout << " map to byte integer " << (int)nucsetbyte << std::endl;

					bytebuffer[0] = nucsetbyte;
					location_file.write((char *) & bytebuffer, 1);
					nucs_in_set = 0;
				}
			}
		}
	}

	uint8_t worldenddominator = 0;				// Initialize the ender byte
	bytebuffer[0] = worldenddominator;
	location_file.write((char *) & bytebuffer, 1);	// Write down the file ending byte

	location_file.close(); // close file
	location_file.clear(); // clear flags
	return file_path;	
}


uint8_t map_4nucs2byte(char nucs[4])
{
	uint8_t single_nuc;
	uint8_t nucbyte = 0;
	
	for(int i = 0; i < 4; i++)
	{
		nucbyte = nucbyte << 2;
		single_nuc = map_nuc2int(nucs[i]);
		if (single_nuc < 4)
		{
			nucbyte = nucbyte | single_nuc;			
		}
	}
	return nucbyte;
}

std::string map_byte2fournucs(uint8_t byte)
{
	int byteint = (int)byte;
	std::string fournucs = "";
	uint8_t three = 3;
	for(int i = 0; i < 4; i++)
	{
		uint8_t next_char = byte & three;
		fournucs = map_int2nuc(next_char) + fournucs;
		byte = byte >> 2;
	}
	//std::cout << "Integer " << byteint << " as nucleotides is " << fournucs << std::endl;
	//std::cout << "Got these nucs: " << fournucs << std::endl;
	return fournucs;
}




std::string format_number(int number)
{
	std::string formatted_number;
	if (number > 99)
	{
		formatted_number = std::to_string(number);
	}
	else if (number > 9)
	{
		formatted_number = "0" + std::to_string(number);
	}
	else 
	{
		formatted_number = "00" + std::to_string(number);	
	}
	return formatted_number;
}


char random_nucleotide()
{
	uint8_t randnucint = static_cast<uint8_t>(std::rand() % 4);
	return map_int2nuc(randnucint);
}
