#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <regex>
#include <string>
#include <tuple>
#include <vector>
#include <deque>
#include <random>
#include <math.h>
#include "fun_kmers.hpp"
#include "fun_consensus.hpp"
#include "file_writer.hpp"
#include "file_reader.hpp"
#include "search_step.hpp"
#include "consensus_step.hpp"


namespace po = boost::program_options;
namespace bfs = boost::filesystem;
//using namespace std;
//using namespace boost::filesystem;


int main(int argc, char *argv[])
{

	/*
		Initialize random
	*/

	int random_seed = 999;
	std::srand(random_seed);


	/*
	 * Step 1: Parse input arguments
	 *
	*/

	// Initialize arguments
	string kmers_path, reads_path, output_path, work_dir, spaced_seed_pattern;
	int kmer_min, buffer_size, delete_files, skip_occurrences, skip_consensus, abs_min_nuc, iterations;
	float relative_nuc_threshold;

	// Parse arguments
	po::options_description desc("k-mer index options");

	desc.add_options()
		("help,h", "Give help")
		("k-mers,k", po::value<std::string>(& kmers_path)->default_value("na"), "Path to the k-mer file")
		("reads,r", po::value<std::string>(& reads_path)->default_value("na"), "Path to the read file")
		("output,o", po::value<std::string>(& output_path)->default_value("out_default.txt"), "Path to the output file")
		("work-dir,w", po::value<std::string>(& work_dir)->default_value("work_tmp"), "Path to a working directory")
		("spaced-seed,s", po::value<std::string>(& spaced_seed_pattern)->default_value("na"), "Spaced seed pattern")
		("min-occ,m", po::value<int>(& kmer_min)->default_value(2), "Minimum number of k-mer occurrences required")
		("min-nuc,n", po::value<int>(& abs_min_nuc)->default_value(2), "Minimum nucleotide threshold")
		("rel-nuc,t", po::value<float>(& relative_nuc_threshold)->default_value(0.1), "Relative nucleotide threshold")
		("buffer-size,b", po::value<int>(& buffer_size)->default_value(5000000), "Buffer size for k-mer occurrences")
		("delete-files,d", po::value<int>(& delete_files)->default_value(0), "Delete temporary location files")
		("skip-occurrences,f", po::value<int>(& skip_occurrences)->default_value(0), "Skip occurrence finding step")
		("skip-consensus,c", po::value<int>(& skip_consensus)->default_value(0), "Skip consensus finding step")
		("iterations,i", po::value<int>(& iterations)->default_value(1), "Iterations");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
 		std::cout << "Print help message here" << std::endl;
    	return 1;
	}

	std::cout << "Initialization started" << std::endl;


	/*
	 * Step 2: Check spaced seed pattern characteristics
	 *
	*/

	// Initialize k-mer stuff
	vector<bool> character_status; // bit vector for spaced seed pattern
	int total_length, fixed_length; // pattern total length, pattern fixed characters
	vector<int> fixed_block_positions, fixed_block_lengths, fixed_nuc_positions; // fixed character block start positions and lengths

	// Get gapped k-mer characteristics, std::tie
	tie(character_status, total_length, fixed_length, fixed_block_positions, fixed_block_lengths, fixed_nuc_positions) = interpret_spaced_seed_pattern(spaced_seed_pattern);


	// Pattern length cannot be even, throw an exception and kill the program
	if (total_length % 2 == 0)
	{
		std::cout << "The spaced seed pattern length cannot be even. Choose a pattern with an odd length." << std::endl;
		return 0;
	}

	std::cout << "Initialization done" << std::endl;

	std::cout << "Long k-mer extraction" << std::endl;
	
	for (int i = 0; i < iterations; i+=1)
	{

		std::cout << "Iteration: " << i+1 << std::endl;

		/*
		 * Step 3: Run search step if it is not skipped
		 *
		*/
		if(skip_occurrences != 1)
		{
			run_search_step(work_dir, kmers_path, reads_path, buffer_size, total_length, fixed_length, character_status, kmer_min, i, iterations);
		}

		/*
		 * Step 4: Run consensus step if it is not skipped
		 *
		*/
		if (skip_consensus != 1)
		{
			run_consensus_step(work_dir, output_path, total_length, fixed_length, relative_nuc_threshold, abs_min_nuc, delete_files, character_status);
		}
	}
	
	/*
	 *
	 * Step 5: End of the program
	 *
	*/	
	std::cout << "Program run finished successfully. Thank you for using LoMeX." << std::endl;
	return 0;

}
