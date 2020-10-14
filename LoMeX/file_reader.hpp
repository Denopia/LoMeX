#ifndef FILE_READER_H
#define FILE_READER_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <random>

using namespace std;

class FastqFileReader
{

private:
	std::ifstream fastq_file;
	std::string fastq_line;
	bool reads_left;
	int read_line;
	int first;
	int last;
	int current_read_number;
	std::deque<char> current_read;

public:
    //FastqFileReader(){}
 
    void initialize_me(std::string file_path, int start_position, int end_position);

	void kill_me();

	void roll_to_next_read();

	bool get_reads_left(){return reads_left;}

	std::deque<char> & get_next_read(){return current_read;}

	int get_read_length(){return current_read.size();}

	bool is_empty(){return current_read.empty();}

	char get_next_char(){return current_read.front();}

	void pop_front_char(){current_read.pop_front();}

};


class TmpFileManager
{

private:
	int id;
	int kmer_length;
	uint8_t occurrence_bytes;
	bool kmer_remains;
	uint8_t bytebuffer[1];
	int sudoku;
	std::ifstream location_file;
	std::vector<bool> is_fixed_character;
	std::string current_spaced_kmer_string;
	std::string current_spaced_kmer_string_from_file;
	__uint128_t current_spaced_kmer_int;
	vector<std::string> current_kmer_occurrences_string;
	std::string my_file_path;
	uint64_t fixed_length;
	
public:

	TmpFileManager(int n, std::string path, int total_length, std::vector<bool> & character_status, int delfile, uint64_t fixed_length_to_save);

	void read_file_info();

	void read_next_line();

	void read_next_line_dupesnt();

	void read_next_line_dupesnt_DEBUG();

	void close_file();

	void delete_file();

	void clear_current_kmer();

	bool get_kmer_remains(){return kmer_remains;}

	__uint128_t get_current_kmer_int(){return current_spaced_kmer_int;}

	std::string get_current_kmer_str(){return current_spaced_kmer_string;}

	vector<std::string> & get_current_kmer_occurrences(){return current_kmer_occurrences_string;}

};



class TmpFileMerger
{

private:
	bool kmer_remains;
	int tmp_kmer_files;
	int my_thread;
	int all_threads;
	int my_iteration;
	int all_iterations;
	uint64_t fixed_length;
	int kmer_length;
	std::vector<bool> is_fixed_character;
	int cleanup;
	__uint128_t current_spaced_kmer;
	__uint128_t previous_spaced_kmer;
	std::string current_spaced_kmer_string;
	std::vector<std::string> current_spaced_kmer_occurrences;
	std::vector<TmpFileManager> files;
	int my_spaced_kmers;

public:

	void initialize_me(int n_files, vector<std::string> file_paths, uint64_t spaced_length, int total_length, std::vector<bool> & character_status, int delfiles, int iterations, int iteration, int threads, int thread);

	void delete_location_files();

	bool solve_next_kmer();

	tuple<std::string, std::vector<std::string> > get_next_kmer();

	std::string get_next_spaced_kmer(){return current_spaced_kmer_string;}

	std::vector<std::string> & get_next_spaced_kmer_regular_kmers(){return current_spaced_kmer_occurrences;}

	bool get_kmer_remains(){return kmer_remains;}
};


int count_fastq_lines(std::string reads_path);


#endif