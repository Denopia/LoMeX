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
#include <random>
#include <math.h>
#include "fun_kmers.hpp"
#include "fun_consensus.hpp"
#include "file_writer.hpp"


namespace po = boost::program_options;
namespace bfs = boost::filesystem;
//using namespace std;
//using namespace boost::filesystem;



/*

	Class for reading fastq files

*/
class FastqFileReader
{
	public:

	std::ifstream fastq_file;
	std::string fastq_line;
	std::vector<char> current_read;
	bool reads_left;
	int read_line;

	/*
	FastqFileReader(std::string file_path)
	{
		reads_left = true;
		read_line = 0;
		fastq_line = "";
		fastq_file.open(file_path, std::ifstream::in);
		//roll_to_next_read();
	}
	*/

	void initialize_me(std::string file_path)
	{
		reads_left = true;
		read_line = 0;
		fastq_line = "";
		fastq_file.open(file_path, std::ifstream::in);
		//roll_to_next_read();
	}

	void kill_me()
	{
		fastq_file.close();
		fastq_file.clear();
	}

	void roll_to_next_read()
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


	std::vector<char> get_next_read()
	{
		return current_read;
	}

	bool get_reads_left()
	{
		return reads_left;
	}

	void flip_reads_left()
	{
		reads_left = !reads_left;
	}
};



/*
	Class for handling a single k-mer location file
*/
class LocationFileManager
{
	public:

	int id;
	int kmer_length;
	uint8_t occurrence_bytes;
	bool kmer_remains;
	uint8_t bytebuffer[1];
	int sudoku;

	std::ifstream location_file;

	std::vector<bool> is_fixed_character;

	std::string current_spaced_kmer_string;
	__uint128_t current_spaced_kmer_int;
	vector<std::string> current_kmer_occurrences_string;

	std::string my_file_path;


	LocationFileManager(int n, std::string path, int total_length, std::vector<bool> & character_status, int delfile)
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


	void read_file_info()
	{	
		//char ok[3];
		//location_file.read(ok, 1);
		//std::cout << "This many bytes in one occurrence: " << bytebuffer[0] << std::endl;
		//location_file.read((char *) & bytebuffer, 1);
		location_file.read((char*)(&occurrence_bytes), sizeof(occurrence_bytes));
		//std::cout << "Bytes per k-mer: " << occurrence_bytes << std::endl;
		//occurrence_bytes = bytebuffer[0];
		//occurrence_bytes = (int)bytebuffer[0];
		//std::cout << "This many bytes in one occurrence: " << occurrence_bytes << std::endl;
	}


	void read_next_line()
	{
		clear_current_kmer();

		std::string debug_last_spaced_kmer = "";
		std::string debug_kmer = "";

		if (!kmer_remains){return;}

		int32_t kmer_occurrences = 0;
		location_file.read((char*) (&kmer_occurrences), sizeof(kmer_occurrences));

		//kmer_occurrences = (int)bytebuffer[0];

		//std::cout << "The current k-mer has this many occurrences: " << kmer_occurrences << std::endl;

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
			debug_kmer = extract_spaced_kmer(current_kmer_occurrence_string, is_fixed_character);
			if (debug_last_spaced_kmer.length() == 0){debug_last_spaced_kmer = debug_kmer;}
			else 
			{
				if (debug_last_spaced_kmer.compare(debug_kmer) != 0)
				{
					std::cout << debug_last_spaced_kmer << std::endl;
					std::cout << debug_kmer << std::endl;
					std::cout << "ERROR WITH BINARY STORED REGULAR K-MERS !!!!";
					exit(9);
				}
			}

		}

		current_spaced_kmer_string = extract_spaced_kmer(current_kmer_occurrences_string[0], is_fixed_character);
		current_spaced_kmer_int = map_str2int(current_spaced_kmer_string);

	}


	void close_file()
	{
		kmer_remains = false;
		location_file.close();
		location_file.clear();
	}

	void delete_file()
	{
		const char * removed_file = my_file_path.c_str();
		int removed = remove(removed_file);
		if (removed == 0)
		{
			std::cout << "Location file " << my_file_path << " permanently deleted" << std::endl;
		}
		else
		{
			std::cout << "Location file deletion was unsuccessful for some reason..." << std::endl;
		}
	}


	void clear_current_kmer()
	{
		current_spaced_kmer_string = "";
		current_spaced_kmer_int = 0;
		current_kmer_occurrences_string.clear();
	}


	bool get_kmer_remains()
	{
		return kmer_remains;
	}


	__uint128_t get_current_kmer_int()
	{
		return current_spaced_kmer_int;
	}

	std::string get_current_kmer_str()
	{
		return current_spaced_kmer_string;
	}


	vector<std::string> get_current_kmer_occurrences(){
		return current_kmer_occurrences_string;
	}

};


/*
	Class for merging multiple k-mer location files
*/
class LocationMerger
{

	public:

	bool kmer_remains;
	int kmers;
	uint64_t fixed_length;
	int kmer_length;
	std::vector<bool> is_fixed_character;
	int cleanup;

	__uint128_t current_spaced_kmer;
	std::string current_spaced_kmer_string;

	//std::vector<int> current_spaced_kmer_locations;
	std::vector<std::string> current_spaced_kmer_occurrences;

	std::vector<LocationFileManager> files;
	/*
	LocationMerger(int n_files, vector<std::string> file_paths, uint64_t spaced_length, int total_length, std::vector<bool> & character_status)
	{
		kmers = n_files;
		kmer_remains = true;
		fixed_length = spaced_length;
		kmer_length = total_length;
		is_fixed_character = character_status;

		//std::cout << "Start building file merger" << std::endl;

		for (int i = 0; i < n_files; i+=1)
		{
			//std::cout << "First file manager for: " << file_paths[i] << std::endl;
			files.push_back(LocationFileManager(i, file_paths[i], kmer_length, is_fixed_character));
			//std::cout << "File manager created successfully " << std::endl;
			
		}
		//std::cout << "Start solving next k-mer" << std::endl;
		solve_next_kmer();
	}
	*/

	void initialize_me(int n_files, vector<std::string> file_paths, uint64_t spaced_length, int total_length, std::vector<bool> & character_status, int delfiles)
	{
		kmers = n_files;
		kmer_remains = true;
		fixed_length = spaced_length;
		kmer_length = total_length;
		is_fixed_character = character_status;
		cleanup = delfiles;

		//std::cout << "Start building file merger" << std::endl;

		for (int i = 0; i < n_files; i+=1)
		{
			if (boost::algorithm::ends_with(file_paths[i], ".bin"))
			{	
				//std::cout << "First file manager for: " << file_paths[i] << std::endl;
				files.push_back(LocationFileManager(i, file_paths[i], kmer_length, is_fixed_character, cleanup));
				//std::cout << "File manager created successfully " << std::endl;
			}
		}
		//std::cout << "Start solving next k-mer" << std::endl;
		solve_next_kmer();
	}


	void delete_location_files()
	{
		std::cout << "Deleting location files" << std::endl;

		for(int i = 0; i < kmers; i+=1)
		{
			files[i].delete_file();
		}
	}


	void solve_next_kmer()
	{
		// First find the smallest k-mer in the input files

		//std::cout << "k-mer solving started" << std::endl;
		__uint128_t smallest_kmer = 0;
		bool first_in = false;
		//std::cout << "First find the smallest k-mer" << std::endl;
		for (int i = 0; i<kmers; i+=1)
		{
			if (files[i].get_kmer_remains())
			{
				__uint128_t file_smallest_kmer = files[i].get_current_kmer_int();
				if (!first_in)
				{
					//std::cout << "At least one valid file found: " << map_int2str(file_smallest_kmer,fixed_length) << std::endl;
					smallest_kmer = file_smallest_kmer;
					first_in = true;
				}
				else if (file_smallest_kmer < smallest_kmer)
				{
					//std::cout << "Even smaller one found: " << map_int2str(file_smallest_kmer,fixed_length) << std::endl;
					smallest_kmer = file_smallest_kmer;
				}
			}
		}
		//std::cout << "Smallest k-mer found: " << std::endl;
		// If no file has more k-mers, take note
		if (!first_in)
		{
			//std::cout << "No file has more k-mers" << std::endl;
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
			//std::cout << "Update current k-mer" << std::endl;

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


	tuple<std::string, std::vector<std::string> > get_next_kmer()
	{
		return make_tuple(current_spaced_kmer_string, current_spaced_kmer_occurrences);
	}


	bool get_kmer_remains()
	{
		return kmer_remains;
	}

};




int main(int argc, char *argv[])
{

	/*
 	*
 	*   # # # # # # # # # # # # # # # 
 	*   #  Phase 1: Initialization  # 
 	*   # # # # # # # # # # # # # # #
 	*
 	*/

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
	int kmer_min, buffer_size, delete_files, skip_occurrences, skip_consensus, abs_min_nuc;
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
		("skip-consensus,c", po::value<int>(& skip_consensus)->default_value(0), "Skip consensus finding step");

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
	vector<int> fixed_block_positions, fixed_block_lengths; // fixed character block start positions and lengths

	// Get gapped k-mer characteristics, std::tie
	tie(character_status, total_length, fixed_length, fixed_block_positions, fixed_block_lengths) = interpret_spaced_seed_pattern(spaced_seed_pattern);


	// Pattern length cannot be even, throw an exception and kill the program
	if (total_length % 2 == 0)
	{
		std::cout << "The spaced seed pattern length cannot be even. Choose a pattern with an odd length." << std::endl;
		return 0;
	}

	std::cout << "Initialization done" << std::endl;

	/*
	 *
	 * INITIALIZE A BUNCH OF STUFF HERE
	 *
	 *
	*/

	std::ifstream infile;
	int qualified_kmers;
	int unqualified_kmers;
	__uint128_t bin_kmer;
	std::string kmer; //k-mer here
	int occurrences; // its number of occurrences here
	map<__uint128_t, int> kmer2occ; // k-mer to its occurrences map here
	std::vector<std::string> kmer_occurrences_file_paths;						// Place to store occurrence file names
	std::string current_file_path;												// The name of the current written file
	int stored_counter;															// Counter for how many k-mer locations have been stored into memory
	int file_counter;															// Counter for how many location files have been created
	int stored_buffer;
	map<__uint128_t, vector<std::string> > kmer2nonfixed;						// k-mer and its occurrences in a map
	__uint128_t current_read_spaced_kmer;										// Current spaced k-mer as integer
	__uint128_t reversed_current_read_spaced_kmer;								// Current spaced k-mer reverse complement as integer
	__uint128_t stored_read_spaced_kmer;										// Current greater spaced k-mer (forward or reverse) as integer
	bool read_spaced_kmer_ok;													// Current spaced k-mer contains no Ns
	uint8_t current_read_char;													// Current read character is ok (is A, C, G or T)
	std::vector<char> read_vector;												// Store current read vector here
	std::string kmer_loc_path;
	int read_counter;
	FastqFileReader fastq_reader; 

	

	if(skip_occurrences == 1)
	{
		goto occurrences_skipped;
	}



	/*
	 *
	 *  # # # # # # # # # # # # # # # # # # # # # # # # # # #
	 *  #  Phase 2: Read spaced k-mers from Squeakr output  #
	 *  # # # # # # # # # # # # # # # # # # # # # # # # # # #
	 *
	*/

	std::cout << "Reading spaced kmers" << std::endl;

	// Start reading spaced k-mers and store them in a the map
	//std::ifstream infile(kmers_path);
	infile.open(kmers_path, std::ifstream::in);

	qualified_kmers = 0; // number of spaced k-mers with enough occurrences
	unqualified_kmers = 0; // number of spaced k-mers with enough occurrences


	while (infile >> kmer >> occurrences)
	{
		if(occurrences >= kmer_min) // if occurrences equal or greater than the minimum number of occurrences
		{
			bin_kmer = map_str2int(kmer); // k-mer as bit representation
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


	/*
	 *
	 *  # # # # # # # # # # # # # # # # # # # # # # # # # # #
	 *  # Phase 4: Find spaced k-mer locations in the reads #
	 *  # # # # # # # # # # # # # # # # # # # # # # # # # # #
	 *
	*/
	
	std::cout << "Reading fastq file reads one by one and writing k-mer occurrences in separate files" << std::endl;

	// Create a fastq file reader
	
	//fastq_reader = FastqFileReader(reads_path);
	fastq_reader.initialize_me(reads_path);

	stored_buffer = buffer_size;
	stored_counter = 0;
	file_counter = 0;
	read_counter = 0;

	// Go through all stored reads
	while(fastq_reader.get_reads_left())
	//for (int i = 0; i < read_vector.size(); i++)
	{
		// Tell file reader to find the next read in the file
		fastq_reader.roll_to_next_read();
		// Fetch the next read
		read_vector = fastq_reader.get_next_read();

		if (read_counter % 100000 == 0){std::cout << "Looking at read number: " << read_counter << std::endl;}
		read_counter += 1;
		//for (char c : read_vector){std::cout << c;}
		//std::cout << std::endl;

		// If read is too small to contain any kmers we need to skip it completely
		if (read_vector.size() < total_length){continue;}

		// Go through all read positions
		for (int ri = 0; ri + total_length - 1 < read_vector.size(); ri++)
		{
			// Solve the current read spaced k-mer
			current_read_spaced_kmer = 0;
			read_spaced_kmer_ok = true;
			for (int j = 0; j < total_length; j++)
			{
				if(character_status[j])
				{
					current_read_char = map_nuc2int(read_vector[ri+j]);
					if (current_read_char > 3)
					{
						read_spaced_kmer_ok = false;
						break;
					}
					current_read_spaced_kmer = current_read_spaced_kmer << 2;
					current_read_spaced_kmer = current_read_spaced_kmer | current_read_char;
				}
			}

			// Current read spaced k-mer solved, check if it is valid and add to map if so

			// Check that there are no N characters
			if (!read_spaced_kmer_ok){continue;}

			// Check that the spaced k-mer was reported by Squeakr
			reversed_current_read_spaced_kmer = reverse_complement_seqint(current_read_spaced_kmer, fixed_length);
			if (kmer2occ.count(current_read_spaced_kmer) == 0 && kmer2occ.count(reversed_current_read_spaced_kmer) == 0) {continue;}

			// Check which is bigger, forward or reverse spaced k-mer
			if (compare_seqs(current_read_spaced_kmer, reversed_current_read_spaced_kmer))
			{
				stored_read_spaced_kmer = current_read_spaced_kmer; // First is bigger
				kmer2nonfixed[stored_read_spaced_kmer].push_back(std::string(read_vector.begin()+ri, read_vector.begin()+ri+total_length));
			}
			else
			{	
				stored_read_spaced_kmer = reversed_current_read_spaced_kmer; // Second is bigger
				std::string fortemp = std::string(read_vector.begin()+ri, read_vector.begin()+ri+total_length);
				std::string fortemprev = reverse_complement_seqstr(fortemp);
				kmer2nonfixed[stored_read_spaced_kmer].push_back(fortemprev);
			}

			// Finally add k-mer to buffer 

			stored_counter += 1;

			if(stored_counter % 100000 == 0){std::cout << "This many occurrences in buffer: " << stored_counter << std::endl;}

			if (stored_counter >= stored_buffer)
			{
				current_file_path = write_occurrences_binary(kmer2nonfixed, work_dir, file_counter, fixed_length, total_length);
				kmer_occurrences_file_paths.push_back(current_file_path);
				stored_counter = 0;
				file_counter += 1;
				kmer2nonfixed.clear();
				std::cout << "File written: " << file_counter << std::endl;
			}
		}
	}

	current_file_path = write_occurrences_binary(kmer2nonfixed, work_dir, file_counter, fixed_length, total_length);
	kmer_occurrences_file_paths.push_back(current_file_path);
	stored_counter = 0;
	file_counter += 1;
	kmer2nonfixed.clear();
	std::cout << "File written: " << file_counter << std::endl;
	
	fastq_reader.kill_me();


occurrences_skipped:


	/*
	 *
	 *  INITIALIZE MORE STUFF IN HERE
	 *
	 *
	*/

	// Initialize long k-mer specific nucleotide occurrences 
	int kmer_nucleotide_occurrences[4][total_length]; // nucleotide position specific counts

	vector<int> consensus_nuc_init(total_length, -1);
	vector<vector<int> > consensus_nucleotides(4, consensus_nuc_init);
	string current_long_kmer, current_spaced_kmer, reverse_long_kmer;
	char curr_kmer_char;
	// Counts for different types of consensus k-mers
	int unambiguous_consensus = 0;
	int simple_ambiguous_consensus = 0;
	int complex_ambiguous_consensus = 0;
	int undecided_consensus = 0;
	bool undecided_kmer;
	std::string consensus_spaced_kmer;
	std::vector<std::string> consensus_spaced_kmer_occurrences;
	int kmer_pattern_occurrences; // counter for how many times a specific spaced k-mer appears in the data
	//cout << "Debug 3" << endl;
	// Go through each gapped k-mer and build the corresponding consensus long k-mer(s)
	kmer_occurrences_file_paths.clear();
	LocationMerger loc_merger;
	int solved_kmer_counter;
	std::ofstream output_file;
	bfs::path wd{work_dir};



	if (skip_consensus == 1)
	{
		goto consensus_skipped;
	}


	/*
	 *
	 *  # # # # # # # # # # # # # # # # # # # # # #
	 *  # Phase 5: Finding long consensus k-mers  #
	 *  # # # # # # # # # # # # # # # # # # # # # #
	 *
	*/

	std::cout << "Finding consensus sequences" << std::endl;
	std::cout << "START" << std::endl;

	// Open a file to write the k-mers

	output_file.open(output_path, ios::out);
	//std::ofstream output_file(output_path);

	//int consensus_nucleotides[4][total_length]; // 1s and 0s only, access to know which characters are trusted by consensus

	//cout << "Debug 1" << endl;

	// Reset k-mer specific nuclotide occurrences
	for (int row = 0; row < 4; row++)
	{			
		for (int column = 0; column < total_length; column++)
		{
			kmer_nucleotide_occurrences[row][column] = 0;
			consensus_nucleotides[row][column] = 0;
		}
	}
	
	// Get all file names

	//path p("data/run/work");

	std::cout << "Finding k-mer location files" << std::endl;

    if(is_directory(wd)) {
        std::cout << wd << " is a directory containing:\n";
        for(auto& entry : boost::make_iterator_range(bfs::directory_iterator(wd), {}))
        {
            std::cout << entry << std::endl;
            kmer_occurrences_file_paths.push_back(entry.path().string());
        }
    }
    else
    {
    	std::cout << "No k-mer location files found" << std::endl;
    	return 0;
    }

    // Create file merger
	
	//LocationMerger loc_merger(kmer_occurrences_file_paths.size(), kmer_occurrences_file_paths, fixed_length, total_length, character_status);

    loc_merger.initialize_me(kmer_occurrences_file_paths.size(), kmer_occurrences_file_paths, fixed_length, total_length, character_status, delete_files);


	std::cout << "File merger successfully created" << std::endl;

	//return 0;

	solved_kmer_counter = 0;
	while(loc_merger.get_kmer_remains())
	{
		// Debug printing
		if (solved_kmer_counter % 10000 == 0)
		{
			std::cout << "Solved k-mers: " << solved_kmer_counter << std::endl;
		}
		
		solved_kmer_counter += 1;

		tie(consensus_spaced_kmer, consensus_spaced_kmer_occurrences) = loc_merger.get_next_kmer(); // FIX THIS SO THAT IT RETURN THE DATA IN CORRECT FORM PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
		loc_merger.solve_next_kmer();

		/*
		std::cout << "Occurrences start" << std::endl;
		for (const std::string &asdf : consensus_spaced_kmer_occurrences)
		{
			std::cout << asdf << std::endl;
		}
		std::cout << "Occurrences end" << std::endl;
		*/
		
		//std::cout << map_int2str(p_spaced_kmer, fixed_length) << std::endl;
		//std::cout << consensus_spaced_kmer << std::endl;

		undecided_kmer = false;

		//current_spaced_kmer = map_int2str(p_spaced_kmer, fixed_length); // get new spaced k-mer
		current_spaced_kmer = consensus_spaced_kmer;

		kmer_pattern_occurrences = consensus_spaced_kmer_occurrences.size(); // reset the new spaced k-mer counter to zero

		//cout << "Debug 4" << endl;

		// Start going through the current spaced k-mer occurrences in the data and fill the nucleotide count matrix

		for (int kmerlen = 0; kmerlen < total_length; kmerlen++)
		{
			for (int occs = 0; occs < consensus_spaced_kmer_occurrences.size(); occs++)
			{
				curr_kmer_char = consensus_spaced_kmer_occurrences[occs].at(kmerlen);
				if (curr_kmer_char == 'C'){kmer_nucleotide_occurrences[0][kmerlen] += 1;}
				else if (curr_kmer_char == 'A'){kmer_nucleotide_occurrences[1][kmerlen] += 1;}
				else if (curr_kmer_char == 'T'){kmer_nucleotide_occurrences[2][kmerlen] += 1;}
				else if (curr_kmer_char == 'G'){kmer_nucleotide_occurrences[3][kmerlen] += 1;}
			}
		}

		/*
		for (int yty = 0; yty < 4; yty++)
		{
			for (int wew = 0; wew < total_length; wew++)
			{
				std::cout << kmer_nucleotide_occurrences[yty][wew] << " ";
			}
			std::cout << "\n";
		}
		*/




		// Now find which characters are trusted in each gap position

		int ambiguous_count = 0; // Count the number of ambiguous consensus positions
		int ambiguous_positions[total_length]; // Array to specify which positions are ambiguous
		int relative_threshold; // Required character count for it to be trusted
		double rcm = relative_nuc_threshold; // Required occurrence proportion 
		int absolute_min_char_count = abs_min_nuc; // Absolute minimum required character count 	

		relative_threshold = static_cast<int>(ceil(kmer_pattern_occurrences * rcm)); 

		//cout << "Debug 6" << endl;

		// Set the minimum count threshold
		if (relative_threshold < absolute_min_char_count)
		{
			relative_threshold = absolute_min_char_count;
		}

		// Reset ambiguous positions to zero
		for (int i = 0; i < total_length; i++)
		{
			ambiguous_positions[i] = 0;
		}

		int trusted_nucs; // How many trusted characters a position has

		// Find out the trusted nuclotides for each gap position

		//cout << "Debug 7" << endl;

		// Loop through all spaced k-mer positions
		for (int consensus_i = 0; consensus_i < total_length; consensus_i++)
		{
			trusted_nucs = 0; // Reset trusted nucleotides
			// Loop through all 4 nucleotides
			for (int nuc_i = 0; nuc_i < 4; nuc_i++)
			{
				if (kmer_nucleotide_occurrences[nuc_i][consensus_i] >= relative_threshold)
				{
					consensus_nucleotides[nuc_i][consensus_i] = 1;
					trusted_nucs += 1;
				}
				if (trusted_nucs == 2)
				{
					//ambiguous_count += 1;
					ambiguous_positions[consensus_i] = 1;
				}
			
			// If none of the nucleotides appears enough (more than once)	
			}
			if (trusted_nucs == 0)
			{
				undecided_kmer = true;
			
			}
		}

		for (const int &ambistatus : ambiguous_positions){ambiguous_count+=ambistatus;}

		//std::cout << "FIND CONSENSUS CASE" << endl;

		// Finally find the consensus sequenses


		if (undecided_kmer)
		{

			//std::cout << "UNDECIDED" << endl;

			undecided_consensus += 1;

			//std::cout << "UNDECIDED FINISHED" << endl;
			goto consensus_cleanup;
		}


		/*
			-- First case: unambiguous consensus --

			Print the consensus k-mer if it is unambiguous
		*/

		else if (ambiguous_count == 0)
		{
			//std::cout << "UNAMBIGUOUS" << endl;
			unambiguous_consensus = determine_unambiguous_consensus(output_file, total_length, consensus_nucleotides, unambiguous_consensus);
			//std::cout << "UNAMBIGUOUS FINISHED" << endl;
		}


		/*
			-- Second case: simple ambiguous consensus --

			i.e. Print the consensus k-mers if they differ in only one position
		*/

		else if (ambiguous_count == 1)
		{
			//std::cout << "SIMPLE AMBIGUOUS" << endl;
			bool unfinished_consensus = true;
			
			while(unfinished_consensus)
			{
				tie(simple_ambiguous_consensus, unfinished_consensus) = determine_simple_ambiguous_consensus(output_file, total_length, consensus_nucleotides, simple_ambiguous_consensus);
			}
			//std::cout << "SIMPLE AMBIGUOUS FINISHED" << endl;
		}


		/*
			-- Third case: complex ambiguous consensus --

			Solve the more ambiguous consensus k-mers with different means
		*/

		else if (ambiguous_count > 1)
		{
			//std::cout << "COMPLEX AMBIGUOUS" << std::endl;
			//std::cout << "AMBIGUOUS COUNT: " << ambiguous_count << std::endl;
			//std::cout << "TOTAL LENGTH: " << total_length<< std::endl;
			vector<vector<char> > ambiguous_patterns;
			ambiguous_patterns = extract_ambiguous_position_matrix(consensus_spaced_kmer_occurrences.size(), ambiguous_count, consensus_spaced_kmer_occurrences, ambiguous_positions, total_length);

			complex_ambiguous_consensus = determine_complex_ambiguous_consensus(output_file, total_length, consensus_nucleotides, kmer_pattern_occurrences, ambiguous_positions,
										  										ambiguous_patterns, ambiguous_count, complex_ambiguous_consensus);
			//std::cout << "COMPLEX AMBIGUOUS FINISHED" << std::endl;
		}

		//std::cout << "SONNA??" << std::endl;

consensus_cleanup:

		// Reset matrices for next iteration
		//std::cout << "RESET MATRIX" << endl;
		for (int row = 0; row < 4; row++)
		{			
			for (int column = 0; column < total_length; column++)
			{
				kmer_nucleotide_occurrences[row][column] = 0;
				consensus_nucleotides[row][column] = 0;
			}
		}
		//std::cout << "RESET MATRIX FINISHED" << endl;
		// Reset k-mer specific occurrence counter
		//kmer_i = 0;
	} 

	output_file.close(); // close file
	output_file.clear(); // clear flags	

	std::cout << "END" << std::endl;
	std::cout << "Print consensus stats" << std::endl;
	std::cout << "Total number of consensus k-mers: " << unambiguous_consensus + simple_ambiguous_consensus + complex_ambiguous_consensus << std::endl;
	std::cout << "Unambiguous consensus k-mers: " << unambiguous_consensus << std::endl;
	std::cout << "Simple ambiguous consensus k-mers: " << simple_ambiguous_consensus << std::endl;
	std::cout << "Complex ambiguous consensus k-mers: " << complex_ambiguous_consensus << std::endl;


consensus_skipped:

	std::cout << "Program run finished successfully. Thank you for using LoMeX." << std::endl;

	return 0;

}
