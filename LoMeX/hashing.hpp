#ifndef HASHING_H
#define HASHING_H

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;



/*

	Custom hash function for spaced k-mers

*/
__uint128_t hash_spaced_kmer(__uint128_t kmer);




/*

	Custom hash comparison function for spaced k-mers

*/
bool equal_spaced_kmers(__uint128_t f, __uint128_t s);

#endif
