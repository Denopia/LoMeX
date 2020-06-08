/*

	File writing functionality

*/

#include "file_writer.hpp"
#include "fun_kmers.hpp"
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
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
	// Define file path
	std::string file_path = work_dir + "/spaced_kmer_occurrences_" + format_number(file_number) + ".bin";
	// Buffer that holds a single byte that will be written to the file
	uint8_t bytebuffer[1];
	// Open the file
	std::ofstream location_file(file_path, ios::out | ios::binary);
	// Get the k-mer total length
	float tl = static_cast<float>(total_length);
	// Determine how many bytes does a single k-mer need to be represented with 2-bit characters
	uint8_t single_occurrence_bytes = static_cast<uint8_t>(ceil(tl / 4.0));							
	// Put the number of needed bytes for a single k-mer in the byte buffer
	bytebuffer[0] = single_occurrence_bytes;
	// Write the content of the byte buffer into the file
	location_file.write((char *) & bytebuffer, 1);
	// Room to store a spaced k-mer
	__uint128_t spaced_kmer;
	// 4 slots to store nuclotides (= 4*2 bits = byte)
	char fournucset[4];
	// Initialize last k-mer variable, for debugging
	//__uint128_t last_kmer = 0;

	for (const auto &pair : kmer2occurrences)		// Go through every spaced k-mer and its occurrences
	{
		spaced_kmer = pair.first;

		//std::cout << "Current spaced k-mer: " << map_int2str(spaced_kmer, fixed_length) << std::endl;
		//if (last_kmer > spaced_kmer)
		//{
		//	std::cout << " -------------------------" << std::endl;
		//	std::cout << "| HORRIBLE ORDERING ERROR |" << std::endl;
		//	std::cout << " -------------------------" << std::endl;
		//}
		//last_kmer = spaced_kmer;

		//kmer_occurrences = pair.second;
		//uint8_t number_of_occurrences = kmer_occurrences.size();
				
		// How many regular k-mers correspond to the current spaced k-mer
		int32_t number_of_occurrences = kmer2occurrences[spaced_kmer].size();

		// Put this here for safety, should not never get in here though
		if(number_of_occurrences == 0){continue;}

		//bytebuffer[0] = number_of_occurrences;
		//location_file.write((char *) & bytebuffer, 1);

		// First, write down the number of regular k-mers 
		location_file.write((char*)(&number_of_occurrences), sizeof(number_of_occurrences));

		// In this loop write every regular k-mer as a set of bytes (number of bytes indicated by the first byte of the file)
		for (int i = 0; i < number_of_occurrences; i+=1)	
		{
			std::string current_occurrence = kmer2occurrences[spaced_kmer][i];

			//std::cout << "Current occurrence: " << current_occurrence << std::endl;
			
			int nucs_in_set = 0;
			// In this loop, write a single occurrence byte by byte
			for (char & c : current_occurrence)		
			{
				fournucset[nucs_in_set] = c;
				if (map_nuc2int(c) > 3)
				{
					fournucset[nucs_in_set] = random_nucleotide();
					//std::cout << "Randomized" << std::endl;
				}

				nucs_in_set += 1;

				if (nucs_in_set == 4)
				{
					uint8_t nucsetbyte = map_4nucs2byte(fournucset);
					bytebuffer[0] = nucsetbyte;
					location_file.write((char *) & bytebuffer, 1);
					nucs_in_set = 0;
				}
			}
			// If the we started a byte but did not write it down yet, fill it with garbage and write it
			while (nucs_in_set != 0)
			{
				fournucset[nucs_in_set] = random_nucleotide();
				nucs_in_set += 1;
				if (nucs_in_set == 4)
				{
					uint8_t nucsetbyte = map_4nucs2byte(fournucset);
					bytebuffer[0] = nucsetbyte;
					location_file.write((char *) & bytebuffer, 1);
					nucs_in_set = 0;
				}
			}
		}
	}

	int32_t worldenddominator = 0;				// Initialize the end marker
	//bytebuffer[0] = worldenddominator;
	//location_file.write((char *) & bytebuffer, 1);	// Write down the file ending byte
	location_file.write((char*)(&worldenddominator), sizeof(worldenddominator));
	location_file.close(); // close file
	location_file.clear(); // clear flags
	return file_path;	
}


std::string write_regular_kmers_ready_binary(map<__uint128_t, vector<uint8_t> > & spaced2regular, std::string work_dir, int file_number, int fixed_length, int total_length)
{
	// Define file path
	std::string file_path = work_dir + "/spaced_kmer_occurrences_" + format_number(file_number) + ".bin";
	// Buffer that holds a single byte that will be written to the file
	uint8_t bytebuffer[1];
	// Open the file
	std::ofstream location_file(file_path, ios::out | ios::binary);
	// Get the k-mer total length
	float tl = static_cast<float>(total_length);
	// Determine how many bytes does a single k-mer need to be represented with 2-bit characters
	uint8_t single_occurrence_bytes = static_cast<uint8_t>(ceil(tl / 4.0));							
	// Put the number of needed bytes for a single k-mer in the byte buffer
	bytebuffer[0] = single_occurrence_bytes;
	// Write the content of the byte buffer into the file
	location_file.write((char *) & bytebuffer, 1);
	// Room to store a spaced k-mer
	__uint128_t spaced_kmer;
	// 4 slots to store nuclotides (= 4*2 bits = byte)
	//char fournucset[4];
	// Initialize last k-mer variable, for debugging
	//__uint128_t last_kmer = 0;

	for (const auto &pair : spaced2regular)		// Go through every spaced k-mer and its occurrences
	{
		spaced_kmer = pair.first;
		float number_of_bytes = spaced2regular[spaced_kmer].size();
		int regular_bytes = static_cast<int>(number_of_bytes);
		float number_of_regulars = number_of_bytes / single_occurrence_bytes;
		int32_t regulars = static_cast<int32_t>(number_of_regulars);

		// First, write down the number of regular k-mers 
		location_file.write((char*)(&regulars), sizeof(regulars));

		// In this loop write every regular k-mer as a set of bytes (number of bytes indicated by the first byte of the file)
		for (int i = 0; i < regular_bytes; i+=1)	
		{
			uint8_t next_byte = spaced2regular[spaced_kmer][i];
			location_file.write((char *) & (next_byte), sizeof(next_byte));
		}
	}
	int32_t worldenddominator = 0;				// Initialize the end marker
	location_file.write((char*)(&worldenddominator), sizeof(worldenddominator));
	location_file.close(); // close file
	location_file.clear(); // clear flags
	return file_path;	
}








void put_kmer_in_buffer_forward(std::vector<char> & read_vector, int ri, int total_length, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer)
{
	uint8_t single_byte = 0;
	char next_nuc;
	uint8_t next_nuc_bits;
	int chars_in_byte = 0;

	for (int i = 0; i < total_length; i++)
	{
		single_byte = single_byte << 2;
		next_nuc = read_vector[ri+i];
		next_nuc_bits = map_nuc2int(next_nuc);
		single_byte = single_byte | next_nuc_bits;
		chars_in_byte += 1;
		if (chars_in_byte == 4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(single_byte);
			single_byte = 0;
			chars_in_byte = 0;
		}
	}
	while (chars_in_byte != 0)
	{
		single_byte = single_byte << 2;
		chars_in_byte += 1;
		if (chars_in_byte == 4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(single_byte);
			single_byte = 0;
			chars_in_byte = 0;
		}
	}
}


void put_kmer_in_buffer_backward(std::vector<char> & read_vector, int ri, int total_length, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer)
{
	uint8_t single_byte = 0;
	char next_nuc;
	uint8_t next_nuc_bits;
	int chars_in_byte = 0;

	for (int i = 0; i < total_length; i++)
	{
		single_byte = single_byte << 2;
		next_nuc = reverse_complement_nucchar(read_vector[ri+total_length-1-i]);
		next_nuc_bits = map_nuc2int(next_nuc);
		single_byte = single_byte | next_nuc_bits;
		chars_in_byte += 1;
		if (chars_in_byte == 4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(single_byte);
			single_byte = 0;
			chars_in_byte = 0;
		}
	}
	while (chars_in_byte != 0)
	{
		single_byte = single_byte << 2;
		chars_in_byte += 1;
		if (chars_in_byte == 4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(single_byte);
			single_byte = 0;
			chars_in_byte = 0;
		}
	}
}

void put_kmer_in_buffer(std::deque<uint8_t> & push_kmer, map<__uint128_t, vector<uint8_t> > & spaced2regular, __uint128_t stored_read_spaced_kmer)
{
	uint8_t push_byte = 0;
	int nucs = 0;
	for (int i = 0; i < push_kmer.size(); i+=1)
	{
		push_byte<<=2;
		push_byte|=push_kmer[i];
		nucs += 1;
		if (nucs==4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(push_byte);
			push_byte = 0;
			nucs = 0;
		}
	}
	while (nucs != 0)
	{
		push_byte<<=2;
		nucs+=1;
		if (nucs == 4)
		{
			spaced2regular[stored_read_spaced_kmer].push_back(push_byte);
			push_byte = 0;
			nucs = 0;
		}
	}
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
	//int byteint = (int)byte;
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
