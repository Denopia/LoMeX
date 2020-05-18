/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include "SqueakrFS.h"
#include <sys/stat.h>
#include <iostream>
#include <algorithm>
#include <dirent.h>

namespace squeakr {
  namespace fs {

		bool has_suffix(const std::string& s, const std::string& suffix)
		{
			return (s.size() >= suffix.size()) && equal(suffix.rbegin(),
																									suffix.rend(), s.rbegin());
		}

		std::string GetDir(std::string str) {
			uint64_t found = str.find_last_of('/');
			return str.substr(0, found);
		}

		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool FileExists(const char* path) {
			struct stat fileStat;
			if (stat(path, &fileStat)) {
				return false;
			}
			if (!S_ISREG(fileStat.st_mode)) {
				return false;
			}
			return true;
		}

		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool DirExists(const char* path) {
			struct stat fileStat;
			if (stat(path, &fileStat)) {
				return false;
			}
			if (!S_ISDIR(fileStat.st_mode)) {
				return false;
			}
			return true;
		}

		void MakeDir(const char* path) { mkdir(path, ACCESSPERMS); }

		std::vector<std::string> GetFilesExt(const char *dir, const char *ext) {
			DIR *folder = opendir(dir);

			if (!folder) {
				std::cerr << "Directory doesn't exist " << dir << std::endl;
				exit(1);
			}

			std::vector<std::string> ret;
			dirent *entry;
			while((entry = readdir(folder)) != NULL)
			{
				if(has_suffix(entry->d_name, ext))
				{
					std::string filename(entry->d_name);
					std::string dirname(dir);
					ret.push_back(std::string(dirname + filename));
				}
			}

			return ret;
		}
  }
}
