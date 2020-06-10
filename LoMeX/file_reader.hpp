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
	std::deque<char> current_read;

public:
    //FastqFileReader(){}
 
    void initialize_me(std::string file_path);

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
	__uint128_t current_spaced_kmer_int;
	vector<std::string> current_kmer_occurrences_string;
	std::string my_file_path;
	
public:

	TmpFileManager(int n, std::string path, int total_length, std::vector<bool> & character_status, int delfile);

	void read_file_info();

	void read_next_line();

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
	int kmers;
	uint64_t fixed_length;
	int kmer_length;
	std::vector<bool> is_fixed_character;
	int cleanup;
	__uint128_t current_spaced_kmer;
	std::string current_spaced_kmer_string;
	std::vector<std::string> current_spaced_kmer_occurrences;
	std::vector<TmpFileManager> files;

public:

	void initialize_me(int n_files, vector<std::string> file_paths, uint64_t spaced_length, int total_length, std::vector<bool> & character_status, int delfiles);

	void delete_location_files();

	void solve_next_kmer();

	tuple<std::string, std::vector<std::string> > get_next_kmer();

	std::string get_next_spaced_kmer(){return current_spaced_kmer_string;}

	std::vector<std::string> & get_next_spaced_kmer_regular_kmers(){return current_spaced_kmer_occurrences;}

	bool get_kmer_remains(){return kmer_remains;}
};


#endif