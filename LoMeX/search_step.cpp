/*

	Search step

*/

#include "fun_consensus.hpp"
#include "fun_kmers.hpp"
#include "file_reader.hpp"
#include "file_writer.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <tuple>

//namespace po = boost::program_options;

using namespace std;


void run_search_step(std::string work_dir, std::string kmers_path, std::string reads_path, int buffer_size, int total_length, int fixed_length, vector<bool> & character_status, int kmer_min, int curite, int maxite)
{
	/*
	 * 1. Read SpacedSeedSqueakr k-mers into memory
	*/

	// Things needed for reading SpacedSeedSqueakr k-mers into memory
	std::ifstream infile;
	int qualified_kmers, unqualified_kmers, occurrences;
	std::string kmer;
	__uint128_t bin_kmer;
	map<__uint128_t, int> kmer2occ;											

	std::cout << "Reading spaced kmers from SpacedSeedSqueakr output" << std::endl;

	// Start reading spaced k-mers and store them in a the map
	infile.open(kmers_path, std::ifstream::in);

	qualified_kmers = 0; // number of spaced k-mers with enough occurrences
	unqualified_kmers = 0; // number of spaced k-mers with enough occurrences


	while (infile >> kmer >> occurrences)
	{
		if(occurrences >= kmer_min) // if occurrences equal or greater than the minimum number of occurrences
		{
			//bin_kmer = map_str2int(kmer); // k-mer as bit representation
			map_str2intR(kmer, bin_kmer); // k-mer as bit representation

			//std::cout << map_int2str(bin_kmer, fixed_length) << std::endl;

			if (bin_kmer % maxite != curite){continue;}

			kmer2occ[bin_kmer] = occurrences; // store k-mer and its occurrences in the map
			qualified_kmers += 1;
		}
		else
		{
			unqualified_kmers += 1;
		}
	}

	infile.close(); // close file
	infile.clear(); // clear flags	

	std::cout << "Number of qualified spaced k-mers: " << qualified_kmers << std::endl;
	std::cout << "Number of unqualified spaced k-mers: " << unqualified_kmers << std::endl;


	
	/*
	 * 2. Find spaced k-mers in the reads
	*/

	// Things needed to search the k-mers in reads
	FastqFileReader fastq_reader;
	int stored_counter, file_counter, read_counter;
	//std::deque<char> read_deque;
	__uint128_t spaced_kmer_as_bits, spaced_kmer_as_bits_rev;
	bool read_spaced_kmer_ok;
	uint8_t current_read_char;
	uint8_t rando_nuc;
	std::deque<bool> nasty_nucs;
	std::deque<uint8_t> long_kmer_as_bits, long_kmer_as_bits_rev;
	map<__uint128_t, vector<uint8_t> > spaced2regular;
	std::string current_file_path;
	//std::vector<std::string> kmer_occurrences_file_paths;


	std::cout << "Reading fastq file reads one by one and writing k-mer occurrences in separate temporary files" << std::endl;

	// Create a fastq file reader
	fastq_reader.initialize_me(reads_path);

	stored_counter = 0;
	file_counter = 0;
	read_counter = 0;

	// Go through all stored reads
	while(fastq_reader.get_reads_left())
	{
		// Tell file reader to find the next read in the file
		fastq_reader.roll_to_next_read();
		// Fetch the next read
		//read_deque = fastq_reader.get_next_read();
		//fastq_reader.get_next_readR(read_deque);

		if (read_counter % 100000 == 0){std::cout << "Looking at read number: " << read_counter << std::endl;}
		read_counter += 1;

		// If read is too small to contain any kmers we need to skip it completely
		//if (read_deque.size() < total_length){continue;}
		if (fastq_reader.get_read_length() < total_length){continue;}

		long_kmer_as_bits.clear();
		long_kmer_as_bits_rev.clear();

		// Go through all read positions

		//while (!read_deque.empty())
		while (!fastq_reader.is_empty())
		{
			spaced_kmer_as_bits = 0;
			spaced_kmer_as_bits_rev = 0;
			read_spaced_kmer_ok = true;
			//current_read_char = map_nuc2int(read_deque.front());
			current_read_char = map_nuc2int(fastq_reader.get_next_char());
			if (current_read_char > 3)
			{
				rando_nuc = random_nucleotide();
				long_kmer_as_bits.push_back(rando_nuc);
				long_kmer_as_bits_rev.push_front(reverse_complement_nucint(rando_nuc));
				nasty_nucs.push_back(true);
			}
			if (current_read_char <= 3)
			{
				long_kmer_as_bits.push_back(current_read_char);
				long_kmer_as_bits_rev.push_front(reverse_complement_nucint(current_read_char));
				nasty_nucs.push_back(false);
			}

			//read_deque.pop_front();
			fastq_reader.pop_front_char();
			if (long_kmer_as_bits.size() < total_length){continue;}
			if (long_kmer_as_bits.size() > total_length)
			{
				long_kmer_as_bits.pop_front();
				long_kmer_as_bits_rev.pop_back();
				nasty_nucs.pop_front();
			}
			for (int ii = 0; ii < total_length; ii+=1)
			{
				if (!character_status[ii]){continue;}
				if (nasty_nucs[ii]){read_spaced_kmer_ok=false;break;}
				spaced_kmer_as_bits <<= 2;
				spaced_kmer_as_bits |= long_kmer_as_bits[ii];
			}
			if (!read_spaced_kmer_ok){continue;}
			spaced_kmer_as_bits_rev = reverse_complement_seqint(spaced_kmer_as_bits, fixed_length);

			if (kmer2occ.count(spaced_kmer_as_bits) == 0 && kmer2occ.count(spaced_kmer_as_bits_rev) == 0) {continue;}

			// Check which is bigger, forward or reverse spaced k-mer, and add to buffer
			if (compare_seqs(spaced_kmer_as_bits, spaced_kmer_as_bits_rev)){put_kmer_in_buffer(long_kmer_as_bits, spaced2regular, spaced_kmer_as_bits);}
			else{put_kmer_in_buffer(long_kmer_as_bits_rev, spaced2regular, spaced_kmer_as_bits_rev);}

			stored_counter += 1;
			if(stored_counter % 100000 == 0){std::cout << "Regular k-mers in buffer: " << stored_counter << std::endl;}

			if (stored_counter >= buffer_size)
			{
				current_file_path = write_regular_kmers_ready_binary(spaced2regular, work_dir, file_counter, fixed_length, total_length);
				//kmer_occurrences_file_paths.push_back(current_file_path);
				stored_counter = 0;
				file_counter += 1;
				spaced2regular.clear();
				std::cout << "File written: " << file_counter << std::endl;
			}
		}
	}

	current_file_path = write_regular_kmers_ready_binary(spaced2regular, work_dir, file_counter, fixed_length, total_length);
	//kmer_occurrences_file_paths.push_back(current_file_path);
	stored_counter = 0;
	file_counter += 1;
	spaced2regular.clear();
	std::cout << "File written: " << file_counter << std::endl;
	
	fastq_reader.kill_me();

	//return kmer_occurrences_file_paths;
	//return 0;
}