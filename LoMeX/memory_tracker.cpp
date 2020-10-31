
#include <vector>
#include <map>
#include <iostream>

//#define DEBUG

// Memory usage
//long long memuse;

// Custom memory allocator to be used with stl containers
#include "memory_tracker.hpp"
#include "file_writer.hpp"

// Annoying, try to get rid of this
using namespace std;

/*
KmerToOccurrenceMap::KmerToOccurrenceMap(long long * buffer_tracker)
{
	tracker_address = buffer_tracker;
	map_allocator = count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > >(tracker_address);
	vector_allocator = count_allocator<uint8_t>(tracker_address);
	kmer2vector_map = std::map<__uint128_t, std::vector<uint8_t, count_allocator<uint8_t> >, std::less<__uint128_t>, count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > > >(map_allocator);
}
*/


void KmerToOccurrenceMap::put_kmer_in_buffer(std::deque<uint8_t> & push_kmer, __uint128_t spaced_kmer)
{
	if (kmer2vector_map.count(spaced_kmer) == 0)
	{
		std::vector<uint8_t, count_allocator<uint8_t> > newkmervec(vector_allocator);
        kmer2vector_map.insert( std::pair<__uint128_t, std::vector<uint8_t, count_allocator<uint8_t> > >(spaced_kmer, newkmervec) );	
	}

	uint8_t push_byte = 0;
	int nucs = 0;
	int pushkmersize = (int)push_kmer.size();
	for (int i = 0; i < pushkmersize; i+=1)
	{
		push_byte<<=2;
		push_byte|=push_kmer[i];

		if (push_kmer[i] < 0 || push_kmer[i] > 3){std::cout << "ILLEGAL CHARACTER GOING TO BUFFER!" << std::endl;}
		
		nucs += 1;
		if (nucs==4)
		{
			kmer2vector_map.at(spaced_kmer).push_back(push_byte);
			//spaced2regular[stored_read_spaced_kmer].push_back(push_byte);
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
			kmer2vector_map.at(spaced_kmer).push_back(push_byte);
			//spaced2regular[stored_read_spaced_kmer].push_back(push_byte);
			push_byte = 0;
			nucs = 0;
		}
	}
}

uint8_t KmerToOccurrenceMap::get_kmer_byte_from_position(__uint128_t spaced_kmer, int position)
{
	return kmer2vector_map.at(spaced_kmer).at(position);
}

int KmerToOccurrenceMap::get_number_of_stored_kmer_bytes(__uint128_t spaced_kmer)
{
	return kmer2vector_map.at(spaced_kmer).size();
}

void KmerToOccurrenceMap::clear_everything()
{
	//std::cout << std::format("BUFFER IS FULL WITH: {}", &tracker_address);
	std::cout << "BUFFER IS FULL WITH: " << tracker_address << " ### " <<*tracker_address << std::endl;
	kmer2vector_map.clear();
	std::cout << "AFTER CLEARING IT IS: " << tracker_address << " ### " <<*tracker_address << std::endl;
	//std::cout << std::format("AFTER CLEARING IT IS: {}", &tracker_address);
}



//
//	MAKE THIS DUPESNT !!!!!!!!!!!!!!!!!!!!!!!!
//
/*
uint64_t KmerToOccurrenceMap::write_buffer_to_disk(std::string & work_dir, int file_number, int fixed_length, int total_length, int thread, vector<bool> & character_status)
{
	float number_of_bytes = 0;
	float number_of_regulars = 0;
	int regular_bytes = 0;
	int32_t regulars = 0;
	uint8_t next_byte = 0;
	uint64_t file_size = 0;
	// Define file path
	std::string file_path = work_dir + "/spaced_kmer_occurrences_" + format_number(thread) + "_" + format_number(file_number) + ".bin";
	// Open the file
	std::ofstream location_file(file_path, ios::out | ios::binary);
	// Get the k-mer total length
	float tl = static_cast<float>(total_length);
	// Determine how many bytes does a single k-mer need to be represented with 2-bit characters
	uint8_t single_occurrence_bytes = static_cast<uint8_t>(ceil(tl / 4.0));							
	// Write the content of the byte buffer into the file
	location_file.write((char *) & (single_occurrence_bytes), sizeof(single_occurrence_bytes));
	// Room to store a spaced k-mer
	__uint128_t spaced_kmer;

	for (const auto &pair : kmer2vector_map)		// Go through every spaced k-mer and its occurrences
	{
		spaced_kmer = pair.first;
		number_of_bytes = kmer2vector_map.at(spaced_kmer).size();
		regular_bytes = static_cast<int>(number_of_bytes);
		number_of_regulars = number_of_bytes / single_occurrence_bytes;
		regulars = static_cast<int32_t>(number_of_regulars);

		// First, write down the number of regular k-mers 
		location_file.write((char*)(&regulars), sizeof(regulars));
		// In this loop write every regular k-mer as a set of bytes (number of bytes indicated by the first byte of the file)
		for (int i = 0; i < regular_bytes; i+=1)	
		{
			next_byte = kmer2vector_map.at(spaced_kmer)[i];
			location_file.write((char *) & (next_byte), sizeof(next_byte));
			file_size+=1;
		}
	}
	int32_t worldenddominator = 0;				// Initialize the end marker
	location_file.write((char*)(&worldenddominator), sizeof(worldenddominator));
	location_file.close(); // close file
	location_file.clear(); // clear flags
	//return make_tuple(file_path, file_size);
	return file_size;	
}
*/

uint64_t KmerToOccurrenceMap::write_buffer_to_disk(std::string & work_dir, int file_number, int fixed_length, int total_length, int thread, vector<bool> & character_status)
{
	uint64_t file_size = 0;
	// Define file path
	std::string file_path = work_dir + "/spaced_kmer_occurrences_" + format_number(thread) + "_" + format_number(file_number) + ".bin";
	// Open the file
	std::ofstream location_file(file_path, ios::out | ios::binary);
	// Get the k-mer total length
	float tl = static_cast<float>(total_length);
	// Determine how many bytes does a single k-mer need to be represented with 2-bit characters
	uint8_t single_occurrence_bytes = static_cast<uint8_t>(ceil(tl / 4.0));							
	// Write the content of the byte buffer into the file
	location_file.write((char *) & (single_occurrence_bytes), sizeof(single_occurrence_bytes));

	__uint128_t spaced_kmer;

	for (const auto &pair : kmer2vector_map)		// Go through every spaced k-mer and its occurrences
	{
		spaced_kmer = pair.first;

		float number_of_bytes = kmer2vector_map.at(spaced_kmer).size();
		float number_of_regulars = number_of_bytes / single_occurrence_bytes;
		int32_t regulars = static_cast<int32_t>(number_of_regulars);

		std::vector<int> skip_these_regular_kmers(regulars, 0);
		std::vector<int> regular_kmer_counts(regulars, 0);

		int32_t non_duplicates = 0;

		for (int curkm = 0; curkm < regulars; curkm++)
		{
			if (skip_these_regular_kmers[curkm] == 1){continue;}

			regular_kmer_counts[curkm] += 1;
			
			for (int othkm = curkm+1; othkm < regulars; othkm++)
			{
				if (skip_these_regular_kmers[othkm] == 1){continue;}
				int issame = 1;
				for (uint bytepos = 0; bytepos < (int)single_occurrence_bytes; bytepos+=1)
				{
					if(kmer2vector_map.at(spaced_kmer)[curkm*single_occurrence_bytes+bytepos] != kmer2vector_map.at(spaced_kmer)[othkm*single_occurrence_bytes+bytepos])
					{
						issame = 0;
						break;
					}
				}

				if (issame == 1)
				{
					skip_these_regular_kmers[othkm] = 1;
					regular_kmer_counts[curkm] += 1;
				}
			}

			if (skip_these_regular_kmers[curkm] == 0){non_duplicates += 1;}
		}

		// First, write down the number of regular k-mers 
		location_file.write((char*)(&non_duplicates), sizeof(non_duplicates));

		// Next the spaced kmer
		//location_file.write((char*)(&spaced_kmer), sizeof(spaced_kmer));

		// In this loop write every regular k-mer as a set of bytes (number of bytes indicated by the first byte of the file)
		uint8_t next_byte;

		int written_kmers = 0;

		for (int i = 0; i < regulars; i++)
		{
			if (skip_these_regular_kmers[i] == 1){continue;}
			written_kmers += 1;

			for (int j = i*single_occurrence_bytes; j < (i+1)*single_occurrence_bytes; j++)
			{
				next_byte = kmer2vector_map.at(spaced_kmer)[j];
				location_file.write((char *)&next_byte, sizeof(next_byte));
			}
			
			// Add count
			int32_t regular_counts_to_file = static_cast<int32_t>(regular_kmer_counts[i]);
			
			location_file.write((char *)&regular_counts_to_file, sizeof(regular_counts_to_file));
			
			file_size += (single_occurrence_bytes*regular_kmer_counts[i]);
		}
	}

	int32_t worldenddominator = 0;				// Initialize the end marker
	//file_size += 4;
	location_file.write((char*)&worldenddominator, sizeof(worldenddominator));
	location_file.close(); // close file
	location_file.clear(); // clear flags

	//return make_tuple(file_path, file_size);
	return file_size;
}