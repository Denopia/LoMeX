/*

	File reading functionality

*/

#include "file_reader.hpp"
#include "fun_kmers.hpp"
#include "file_writer.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <random>
#include <math.h>

using namespace std;


/*
	Class for reading fastq files
*/

//FastqFileReader::FastqFileReader();

void FastqFileReader::initialize_me(std::string file_path)
{
	reads_left = true;
	read_line = 0;
	fastq_line = "";
	fastq_file.open(file_path, std::ifstream::in);
}

void FastqFileReader::kill_me()
{
	fastq_file.close();
	fastq_file.clear();
}

void FastqFileReader::roll_to_next_read()
{
	// Read fastq file line by line
	bool new_read_gotten = false;
	char current_nucleotide;
	char push_nucleotide;
	current_read.clear();
	while (std::getline(fastq_file, fastq_line, '\n'))
	{
		// Line with the actual read encountered, it is the second line
		if (read_line == 1)
		{
			for (int i = 0; i < fastq_line.length(); i++)
			{
				current_nucleotide = fastq_line.at(i);
				if (current_nucleotide == 'C' || current_nucleotide == 'c'){push_nucleotide = 'C';}
				else if (current_nucleotide == 'A' || current_nucleotide == 'a'){push_nucleotide = 'A';}
				else if (current_nucleotide == 'T' || current_nucleotide == 't'){push_nucleotide = 'T';}
				else if (current_nucleotide == 'G' || current_nucleotide == 'g'){push_nucleotide = 'G';}
				else {push_nucleotide = 'N';}
				current_read.push_back(push_nucleotide);
			}
			new_read_gotten = true;
			reads_left = true;
		}
		// Update line counters
		read_line += 1;
		if (read_line == 4)
		{
			read_line = 0;
		}
		if (new_read_gotten){break;}
	}
	if (!new_read_gotten)
	{
		reads_left = false;
		kill_me();
	}
}



/*
	Class for handling a single temporary k-mer file
*/
TmpFileManager::TmpFileManager(int n, std::string path, int total_length, std::vector<bool> & character_status, int delfile)
{
	is_fixed_character = character_status;
	kmer_remains = true;
	kmer_length = total_length;
	my_file_path = path;
	sudoku = delfile;
	location_file.open(path, ios::in | ios::binary);
	read_file_info();
	read_next_line();	
}

void TmpFileManager::read_file_info()
{	
	location_file.read((char*)(&occurrence_bytes), sizeof(occurrence_bytes));
}


void TmpFileManager::read_next_line()
{
	clear_current_kmer();

	//std::string debug_last_spaced_kmer = "";
	//td::string debug_kmer = "";

	if (!kmer_remains){return;}

	int32_t kmer_occurrences = 0;
	location_file.read((char*) (&kmer_occurrences), sizeof(kmer_occurrences));

	if(kmer_occurrences == 0)
	{
		close_file();
		if (sudoku == 1)
		{
			delete_file();
		}
		return;
	}

	uint8_t occbyte;
	std::string current_kmer_occurrence_string;
	for (int occi = 0; occi < kmer_occurrences; occi+=1)
	{
		current_kmer_occurrence_string = "";
		for (int occb = 0; occb < occurrence_bytes; occb+=1)
		{
			location_file.read((char *) & bytebuffer, 1);
			occbyte = bytebuffer[0];
			current_kmer_occurrence_string = current_kmer_occurrence_string + map_byte2fournucs(occbyte);
		}
		current_kmer_occurrence_string = current_kmer_occurrence_string.substr(0, kmer_length);
		current_kmer_occurrences_string.push_back(current_kmer_occurrence_string);

		// DEBUGGING, can be removed if no problems occur soon
		//debug_kmer = extract_spaced_kmer(current_kmer_occurrence_string, is_fixed_character);
		//if (debug_last_spaced_kmer.length() == 0){debug_last_spaced_kmer = debug_kmer;}
		//else 
		//{
		//	if (debug_last_spaced_kmer.compare(debug_kmer) != 0)
		//	{
		//		std::cout << debug_last_spaced_kmer << std::endl;
		//		std::cout << debug_kmer << std::endl;
		//		std::cout << "ERROR WITH BINARY STORED REGULAR K-MERS !!!!";
		//		exit(9);
		//	}
		//}
	}
	current_spaced_kmer_string = extract_spaced_kmer(current_kmer_occurrences_string[0], is_fixed_character);
	current_spaced_kmer_int = map_str2int(current_spaced_kmer_string);
}


void TmpFileManager::close_file()
{
	kmer_remains = false;
	location_file.close();
	location_file.clear();
}

void TmpFileManager::delete_file()
{
	const char * removed_file = my_file_path.c_str();
	int removed = remove(removed_file);
	if (removed == 0) {std::cout << "Location file " << my_file_path << " permanently deleted" << std::endl;}
	else {std::cout << "Location file deletion was unsuccessful for some reason..." << std::endl;}
}

void TmpFileManager::clear_current_kmer()
{
	current_spaced_kmer_string = "";
	current_spaced_kmer_int = 0;
	current_kmer_occurrences_string.clear();
}


/*
	Class for merging multiple temporary k-mer files
*/

//TmpFileMerger::TmpFileMerger();


void TmpFileMerger::initialize_me(int n_files, vector<std::string> file_paths, uint64_t spaced_length, int total_length, std::vector<bool> & character_status, int delfiles)
{
	kmers = n_files;
	kmer_remains = true;
	fixed_length = spaced_length;
	kmer_length = total_length;
	is_fixed_character = character_status;
	cleanup = delfiles;

	for (int i = 0; i < n_files; i+=1)
	{
		if (boost::algorithm::ends_with(file_paths[i], ".bin"))
		{	
			files.push_back(TmpFileManager(i, file_paths[i], kmer_length, is_fixed_character, cleanup));
		}
	}
	//solve_next_kmer();
}

void TmpFileMerger::delete_location_files()
{
	std::cout << "Deleting location files" << std::endl;
	for(int i = 0; i < kmers; i+=1){files[i].delete_file();}
}


void TmpFileMerger::solve_next_kmer()
{
	// First find the smallest k-mer in the input files
	__uint128_t smallest_kmer = 0;
	bool first_in = false;
	for (int i = 0; i<kmers; i+=1)
	{
		if (files[i].get_kmer_remains())
		{
			__uint128_t file_smallest_kmer = files[i].get_current_kmer_int();
			if (!first_in)
			{
				smallest_kmer = file_smallest_kmer;
				first_in = true;
			}
			else if (file_smallest_kmer < smallest_kmer)
			{
				smallest_kmer = file_smallest_kmer;
			}
		}
	}
	// If no file has more k-mers, take note
	if (!first_in)
	{
		kmer_remains = false;
		current_spaced_kmer = 0;
		current_spaced_kmer_occurrences.clear();
	}
	// Otherwise update the stored smallest k-mer and its locations
	else
	{
		current_spaced_kmer = smallest_kmer;
		current_spaced_kmer_string = map_int2str(current_spaced_kmer, fixed_length);
		current_spaced_kmer_occurrences.clear();

		for (int i = 0; i<kmers; i+=1)
		{
			bool file_kmer_remains = files[i].get_kmer_remains();
			if (file_kmer_remains)
			{
				__uint128_t file_current_kmer = files[i].get_current_kmer_int();
				if (file_current_kmer == current_spaced_kmer)
				{
					std::vector<std::string> append_occurrences = files[i].get_current_kmer_occurrences();
					// Tell the file manager to read next line discarding the currently asked information
					int info_amount = append_occurrences.size();

					for (int info_counter = 0; info_counter < info_amount; info_counter+=1)
					{
						current_spaced_kmer_occurrences.push_back(append_occurrences[info_counter]);
					}
					files[i].read_next_line();
				}
			}
		}
	}
}

tuple<std::string, std::vector<std::string> > TmpFileMerger::get_next_kmer()
{
	return make_tuple(current_spaced_kmer_string, current_spaced_kmer_occurrences);
}

