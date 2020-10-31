#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

#include <iostream>
#include <deque>

using namespace std; 


template <typename T>
class count_allocator: public std::allocator<T> {

private:
  long long *mem_tracker;

public:
  
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  template<typename _Tp1>
  struct rebind {
    typedef count_allocator<_Tp1> other;
  };

  pointer allocate(size_type n, const void *hint=0) {
#ifdef DEBUG
    fprintf(stderr, "Alloc %lu bytes.\n", n*sizeof(T));
#endif
    //memuse += n*sizeof(T);
    *mem_tracker += n*sizeof(T);
#ifdef DEBUG
    std::cout << "Now I have this much memory used " <<  *mem_tracker << " ::: " <<mem_tracker << std::endl;
#endif
    return std::allocator<T>::allocate(n, hint);
  }

  void deallocate(pointer p, size_type n) {
#ifdef DEBUG
    fprintf(stderr, "Dealloc %lu bytes (%p).\n", n*sizeof(T), p);
#endif
    //memuse -= n*sizeof(T);
    *mem_tracker -= n*sizeof(T);
#ifdef DEBUG
    std::cout << "Now I have this much memory used " <<  *mem_tracker << " ::: " <<mem_tracker << std::endl;
#endif
    return std::allocator<T>::deallocate(p, n);
  }

  long long * get_mem_tracker() const
  {
#ifdef DEBUG
    std::cout << "Hello I am getting my memory usage through allocator access:  " << *mem_tracker << std::endl;
#endif
    return mem_tracker;
  }

  count_allocator(long long * mem_pointer) throw(): std::allocator<T>() {
    mem_tracker = mem_pointer;
#ifdef DEBUG
    fprintf(stderr, "Hello allocator!\n"); 
    //std::cout << "I am initialized with this much memory usage " << *mem_tracker << " ::: " << mem_tracker << "\n";
#endif
  }

  count_allocator(const count_allocator &a) throw(): std::allocator<T>(a) {
    //mem_tracker = a.mem_tracker;
    mem_tracker = a.get_mem_tracker();
  }
  
  template <typename U>                    
  count_allocator(const count_allocator<U> &a) throw(): std::allocator<T>(a) {
    //mem_tracker = a.mem_tracker;
    mem_tracker = a.get_mem_tracker();
  }
  ~count_allocator() throw() { }
};



class KmerToOccurrenceMap {

private:

	long long *tracker_address;
	count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > > map_allocator;;
	count_allocator<uint8_t> vector_allocator;
	std::map<__uint128_t, std::vector<uint8_t, count_allocator<uint8_t> >, std::less<__uint128_t>, count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > > > kmer2vector_map;


public:

	KmerToOccurrenceMap(long long * buffer_tracker) : tracker_address(buffer_tracker), map_allocator(count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > >(buffer_tracker)), vector_allocator(count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > >(buffer_tracker)),kmer2vector_map(std::map<__uint128_t, std::vector<uint8_t, count_allocator<uint8_t> >, std::less<__uint128_t>, count_allocator<std::pair<__uint128_t, std::vector<uint8_t> > > >(buffer_tracker))
	{
	
	}

	void put_kmer_in_buffer(std::deque<uint8_t> & push_kmer, __uint128_t spaced_kmer);

	uint8_t get_kmer_byte_from_position(__uint128_t spaced_kmer, int position);

	int get_number_of_stored_kmer_bytes(__uint128_t spaced_kmer);

	uint64_t write_buffer_to_disk(std::string & work_dir, int file_counter, int fixed_length, int total_length, int thread, vector<bool> & character_status);

	void clear_everything();

};



#endif

