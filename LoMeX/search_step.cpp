/*

	Search step

*/

#include "fun_consensus.hpp"
#include "fun_kmers.hpp"
#include "file_reader.hpp"
#include "file_writer.hpp"
#include "hashing.hpp"
#include "memory_tracker.hpp"
#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <set>

//namespace po = boost::program_options;

using namespace std;


tuple<uint64_t, int> run_search_step(std::string work_dir, std::string kmers_path, std::string reads_path, int buffer_size, int total_length, int fixed_length, vector<bool> & character_status, int kmer_min, int curite, int maxite, int threads, int number_of_reads)
{

	long long USED_MEMORY_K2A = 0;

	count_allocator<std::pair<__uint128_t, int>> allocator_kmer2abundance(&USED_MEMORY_K2A);


	int written_kmers_count = 0;
	/*
	 * 1. Read SpacedSeedSqueakr k-mers into memory
	*/

	// Things needed for reading SpacedSeedSqueakr k-mers into memory
	std::ifstream infile;
	int qualified_kmers, unqualified_kmers, occurrences;
	std::string kmer;
	__uint128_t bin_kmer;
	//map<__uint128_t, int> kmer2occ;

	//int buckets = 1000000000; // --- TODO --- What is the best value for this??
	//uint buckets = 2147483647;
	int buckets = 5;
	
	//std::unordered_map<__uint128_t, int, std::function<__uint128_t(__uint128_t)>, std::function<bool(__uint128_t, __uint128_t)> > kmer2occ(buckets, hash_spaced_kmer, equal_spaced_kmers);
	

	// BELOW WORKS
	//std::unordered_map<__uint128_t, __uint128_t, std::function<__uint128_t(__uint128_t)>, std::function<bool(__uint128_t, __uint128_t)> > kmer2occ(buckets, hash_spaced_kmer, equal_spaced_kmers);
	std::unordered_map<__uint128_t, int, std::function<__uint128_t(__uint128_t)>, std::function<bool(__uint128_t, __uint128_t)>, count_allocator<std::pair<__uint128_t, int> > > kmer2occ(buckets, hash_spaced_kmer, equal_spaced_kmers, allocator_kmer2abundance);
	
	//std::map<int, std::vector<int, count_allocator<int> >, std::less<int>, count_allocator<std::pair<int, std::vector<int> > > > mymap3(ca2);


	//std::set<__uint128_t> kmer2occ2;


	//std::unordered_map< __uint128_t, __uint128_t > kmer2occ(5, hash_spaced_kmer, equal_spaced_kmers);

	std::cout << "Reading spaced kmers from SpacedSeedSqueakr output" << std::endl;

	// Start reading spaced k-mers and store them in a the map
	infile.open(kmers_path, std::ifstream::in);

	qualified_kmers = 0; // number of spaced k-mers with enough occurrences
	unqualified_kmers = 0; // number of spaced k-mers with enough occurrences


	int all_spaced_kmers = 0;
	int current_spaced_kmers = 0;


	

	while (infile >> kmer >> occurrences)
	{
		if(occurrences >= kmer_min) // if occurrences equal or greater than the minimum number of occurrences
		{
			all_spaced_kmers += 1;
			//bin_kmer = map_str2int(kmer); // k-mer as bit representation
			map_str2intR(kmer, bin_kmer); // k-mer as bit representation

			__uint128_t bin_kmer_rev = reverse_complement_seqint(bin_kmer, static_cast<uint64_t>(fixed_length));

			__uint128_t add_kmer;
			if (bin_kmer < bin_kmer_rev)
			{
				add_kmer = bin_kmer;
				//std::cout << "NO KMER SWAP AHAHA" << std::endl;
			}
			else
			{
				add_kmer = bin_kmer_rev;	
			}
			
			//int tester = (int)bin_kmer;
			//if (tester < 0) {tester = 0 - tester;}
			//std::cout << bin_kmer << std::endl;
			//std::cout << tester << std::endl;
			//std::cout << "Tester: " << tester << std::endl;
			//std::cout << "Modulo: " << tester % maxite << std::endl;
			//std::cout << map_int2str(bin_kmer, fixed_length) << std::endl;
			//if (bin_kmer % maxite != tester % maxite){std::cout << "AY DIOS MIO!!" << std::endl;}
			//if (bin_kmer % maxite != curite){continue;}
			if (add_kmer % maxite < 0)
			{
				std::cout << "MODULO FAILURE!" << std::endl;
			}
			if (add_kmer % maxite == curite)
			{
				kmer2occ[add_kmer] = occurrences; // store k-mer and its occurrences in the map
				//kmer2occ2.insert(bin_kmer);
				qualified_kmers += 1;
				current_spaced_kmers += 1;
			}
		}
		else
		{
			unqualified_kmers += 1;
		}
	}

	//std::cout << "THIS ITERATION HAS: " << current_spaced_kmers << " SPACED K-MERS OUT OF ALL QUALIFIED " << all_spaced_kmers << " K-MERS." << std::endl;


	infile.close(); // close file
	infile.clear(); // clear flags	

	std::cout << "This many bytes used to store k-mers and their counts: " << USED_MEMORY_K2A << std::endl;


	std::cout << "Number of qualified spaced k-mers: " << qualified_kmers << std::endl;
	//std::cout << "Number of unqualified spaced k-mers: " << unqualified_kmers << std::endl;

	//std::cout << "The number of reads is: " << number_of_reads << std::endl;
	std::cout << "Splitting the job between " << threads << " threads" << std::endl;

	int break_points[threads+1];
	break_points[0] = 0;
	int read_block = std::floor((double)number_of_reads / (double)threads);
	for (int g = 1; g < threads+1; g+=1){break_points[g] = break_points[g-1] + read_block;}
	break_points[threads] = number_of_reads;

	//for (int & boo : break_points){std::cout << boo << "-";}

	uint64_t file_sizes[threads];


	long long TOTAL_BUFFER_SIZE = buffer_size;
	TOTAL_BUFFER_SIZE -= USED_MEMORY_K2A;

	long long THREAD_BUFFER_SIZE = TOTAL_BUFFER_SIZE / threads;

	std::cout << "BUFFER SIZE FOR ONE THREAD IS " << THREAD_BUFFER_SIZE << std::endl;

	if (THREAD_BUFFER_SIZE < 1000) 
	{
		exit(1);
	}

	#pragma omp parallel for
	for (int thread = 0; thread < threads; thread += 1)
	{

		// ======================================================
		// ========== ALLOCATOR STUFF BEING DONE BELOW ==========

		long long THREAD_BUFFER_USAGE = 0;

		KmerToOccurrenceMap mykmer2occ(&THREAD_BUFFER_USAGE);


		//count_allocator<std::pair<int, std::vector<uint8_t> > > thread_allocator_kmer2occurrences(&THREAD_BUFFER_USAGE);

		//std::map<__uint128_t, std::vector<uint8_t, count_allocator<uint8_t> >, count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > > > spaced2regular(thread_allocator_kmer2occurrences);


		// ========== ALLOCATOR STUFF BEING DONE ABOVE ==========
		// ======================================================

		file_sizes[thread] = 0;
		//std::cout << "Started thread: " << thread << std::endl;
		/*
		 * 2. Find spaced k-mers in the reads
		*/
		// Things needed to search the k-mers in reads
		FastqFileReader fastq_reader;
		int stored_counter, file_counter, read_counter;
		__uint128_t spaced_kmer_as_bits, spaced_kmer_as_bits_rev;
		bool read_spaced_kmer_ok;
		uint8_t current_read_char;
		uint8_t rando_nuc;
		std::deque<bool> nasty_nucs;
		std::deque<uint8_t> long_kmer_as_bits, long_kmer_as_bits_rev;

		// REMADE THIS WITH CUSTOM ALLOCATOR
		//std::map<__uint128_t, vector<uint8_t> > spaced2regular;

		std::string current_file_path;
		int my_buffer;
		int file_size;

		//std::cout << "Reading fastq file reads one by one and writing k-mer occurrences in separate temporary files" << std::endl;

		// Create a fastq file reader
		int start_position = break_points[thread];
		int end_position = break_points[thread+1];
		fastq_reader.initialize_me(reads_path, start_position, end_position);

		my_buffer = (int)std::floor((double)buffer_size / (double)threads);
		stored_counter = 0;
		file_counter = 0;
		read_counter = 0;

		// Go through all stored reads
		while(fastq_reader.get_reads_left())
		{
			// Tell file reader to find the next read in the file
			fastq_reader.roll_to_next_read();
			read_counter += 1;

			if ((int)fastq_reader.get_read_length() < total_length){continue;}

			long_kmer_as_bits.clear();
			long_kmer_as_bits_rev.clear();
			nasty_nucs.clear();

			// Go through all read positions
			while (!fastq_reader.is_empty())
			{
				spaced_kmer_as_bits = 0;
				spaced_kmer_as_bits_rev = 0;
				read_spaced_kmer_ok = true;
				//current_read_char = map_nuc2int(read_deque.front());
				current_read_char = map_nuc2int(fastq_reader.get_next_char());

				//if (current_read_char > 3){std::cout << "CHARACTER FROM READ TOO BIG" << std::endl;}
	
				if (current_read_char > 3)
				{
					//rando_nuc = map_nuc2int(random_nucleotide());
					rando_nuc = random_nucleotide_int();
					
					long_kmer_as_bits.push_back(rando_nuc);
					long_kmer_as_bits_rev.push_front(reverse_complement_nucint(rando_nuc));
					//if (rando_nuc < 0){std::cout << "RANDOM NUC IS NEGATIVE" << std::endl;}
					//if (rando_nuc > 3){std::cout << "RANDOM NUC IS TOO BIG" << std::endl;}

					nasty_nucs.push_back(true);
				}
				else if (current_read_char <= 3)
				{
					//if (current_read_char < 0){std::cout << "CHARACTER FROM READ IS NEGATIVE" << std::endl;}
					long_kmer_as_bits.push_back(current_read_char);
					long_kmer_as_bits_rev.push_front(reverse_complement_nucint(current_read_char));
					nasty_nucs.push_back(false);
				}

				fastq_reader.pop_front_char();

				if ((int)long_kmer_as_bits.size() < total_length){continue;}
				
				if ((int)long_kmer_as_bits.size() > total_length)
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

				if (spaced_kmer_as_bits <= spaced_kmer_as_bits_rev)
				{
					//put_kmer_in_buffer(long_kmer_as_bits, spaced2regular, spaced_kmer_as_bits);
					mykmer2occ.put_kmer_in_buffer(long_kmer_as_bits, spaced_kmer_as_bits);
				}
				else
				{
					//put_kmer_in_buffer(long_kmer_as_bits_rev, spaced2regular, spaced_kmer_as_bits_rev);
					mykmer2occ.put_kmer_in_buffer(long_kmer_as_bits_rev, spaced_kmer_as_bits_rev);
				}

				stored_counter += 1;
				
				if (stored_counter >= my_buffer)
				//if (THREAD_BUFFER_USAGE >= THREAD_BUFFER_SIZE)
				{
					//tie(current_file_path, file_size) = write_regular_kmers_ready_binary(spaced2regular, work_dir, file_counter, fixed_length, total_length, thread);

					//tie(current_file_path, file_size) = write_regular_kmers_ready_binary_dupesnt(spaced2regular, work_dir, file_counter, fixed_length, total_length, thread, character_status);
					//spaced2regular.clear();
					
					file_size = mykmer2occ.write_buffer_to_disk(work_dir, file_counter, fixed_length, total_length, thread, character_status);
					mykmer2occ.clear_everything();

					file_sizes[thread] += file_size;
					written_kmers_count += stored_counter;
					stored_counter = 0;
					file_counter += 1;	
				}
			}
		}
		//tie(current_file_path, file_size) = write_regular_kmers_ready_binary(spaced2regular, work_dir, file_counter, fixed_length, total_length, thread);
		if (stored_counter > 0)
		{
			//tie(current_file_path, file_size) = write_regular_kmers_ready_binary_dupesnt(spaced2regular, work_dir, file_counter, fixed_length, total_length, thread, character_status);
			//spaced2regular.clear();

			file_size = mykmer2occ.write_buffer_to_disk(work_dir, file_counter, fixed_length, total_length, thread, character_status);
			mykmer2occ.clear_everything();

			file_sizes[thread] += file_size;
			written_kmers_count += stored_counter;
			stored_counter = 0;
			file_counter += 1;
		}

		fastq_reader.kill_me();
	}

	// FILE SIZE DEBUGGING STUFF
	/* 
	int total_file_size = 0;
	for (int i = 0; i < threads; i++)
	{
		total_file_size += file_sizes[i];
	}
	*/

	//return make_tuple(total_file_size, written_kmers_count);	
	//return kmer_occurrences_file_paths;
	//return 0;
	return make_tuple(0, written_kmers_count);
}