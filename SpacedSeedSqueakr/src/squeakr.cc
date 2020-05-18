/*
 * This is a modified version of the original file appearing in Squeakr. 
 * Modification by Miika Leinonen
*/

/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <exception>

#include "ProgOpts.h"
#include "clipp.h"
#include "SqueakrFS.h"

template <typename T>
void explore_options_verbose(T& res) {
  if(res.any_error()) { std::cerr << "error\n"; }

  //aggregated errors
  if(res.unmapped_args_count()) { std::cerr << "error unmapped args count\n"; /* ... */ }
  if(res.any_bad_repeat()) { std::cerr << "error bad repeat \n"; /* ... */ }
  if(res.any_blocked())    { std::cerr << "error blocked \n"; /* ... */ }
  if(res.any_conflict())   { std::cerr << "error conflict\n"; /* ... */ }

  for(const auto& m : res.missing()) { 
    std::cerr << "missing " << m.param() << " after index " << m.after_index() << '\n';
  }

  //per-argument mapping
  for(const auto& m : res) {
    std::cerr << m.index() << ": " << m.arg() << " -> " << m.param();
    std::cerr << " repeat #" << m.repeat();
    if(m.blocked()) std::cerr << " blocked";
    if(m.conflict()) std::cerr << " conflict";
    std::cerr<< '\n';
  }
}

int query_main (QueryOpts& opt);
int count_main (CountOpts& opt);
int inner_prod_main (InnerProdOpts& opt);
int list_main (ListOpts& opt);
int info_main (InfoOpts& opt);

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int main ( int argc, char *argv[] ) {
  using namespace clipp;
  enum class mode {count, query, inner_prod, list, info, help};
  mode selected = mode::help;

  auto console = spdlog::stdout_color_mt("squeakr_console");

  CountOpts countopt;
  QueryOpts queryopt;
  InnerProdOpts innerprodopt;
	ListOpts listopt;
	InfoOpts infoopt;
  countopt.console = console;
  queryopt.console = console;
  innerprodopt.console = console;
	listopt.console = console;
	infoopt.console = console;

  auto ensure_file_exists = [](const std::string& s) -> bool {
    bool exists = squeakr::fs::FileExists(s.c_str());
    if (!exists) {
      std::string e = "The required input file " + s + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

  auto ensure_parent_dir_exists = [](const std::string& s) -> bool {
		std::string parent_dir = squeakr::fs::GetDir(s);
    bool exists = squeakr::fs::DirExists(parent_dir.c_str());
    if (!exists) {
      std::string e = "The required input directory " + parent_dir + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

	auto enusure_size_is_specified = [](const CountOpts countopt) -> bool {
		if (!countopt.setqbits &&  countopt.numthreads == 0) {
			std::string e = "Size option is required if the thread count is greater than 1.";
			throw std::runtime_error{e};
		}
		return true;
	};

	auto count_mode = (
									command("count").set(selected, mode::count),

									required("-p", "--pattern") & value("pattern", countopt.pattern) %
									"gapped k-mer pattern",

									option("-e", "--exact").set(countopt.exact, 1) %
									"squeakr-exact (default is Squeakr approximate)",

									required("-k","--kmer") & value("k-size", countopt.ksize) %
									"length of k-mers to count",

									option("-c","--cutoff") & value("cutoff", countopt.cutoff) %
									"only output k-mers with count greater than or equal to cutoff (default = 1)",

									option("-n","--no-counts").set(countopt.contains_counts, 0) %
									"only output k-mers and no counts (default = false)",

									option("-s","--log-slots").set(countopt.setqbits, true) &
									value("log-slots", countopt.qbits) % "log of number of slots in the CQF. (Size argument is only optional when numthreads is exactly 1.)",
									option("-t","--threads") & value("num-threads",
																									 countopt.numthreads) %
									"number of threads to use to count (default = number of hardware threads)",
									required("-o","--output-file") &
									value(ensure_parent_dir_exists, "out-file",
												countopt.output_file) %
									"file in which output should be written",
									values(ensure_file_exists, "files", countopt.filenames) % "list of files to be counted (supported files: fastq and compressed gzip or bzip2 fastq files)"
									//option("-h", "--help")      % "show help"
						 );

	auto query_mode = (
							command("query").set(selected, mode::query),
							required("-f", "--squeakr-file") & value(ensure_file_exists, "squeakr-file",
																									 queryopt.squeakr_file) % "input squeakr file",
							required("-q","--query-file") & value(ensure_file_exists, "query-file",
																								queryopt.queryfile) % "input query file",
							required("-o", "--output-file") &
							value(ensure_parent_dir_exists, "output-file",
										queryopt.output_file) % "output file"
							//option("-h", "--help")  % "show help"
						 );

	auto inner_prod_mode = (
							command("inner_prod").set(selected, mode::inner_prod),
							value(ensure_file_exists, "first-input",
																												 innerprodopt.squeakr_filea)
							% "first input squeakr file",
							value(ensure_file_exists, "second-input",
																													innerprodopt.squeakr_fileb)
							% "second input squeakr file"
							//option("-h", "--help")  % "show help"
						 );

	auto list_mode = (
							command("list").set(selected, mode::list),
							required("-f", "--squeakr-file-file") & value(ensure_file_exists, "squeakr-file",
																									 listopt.squeakr_file) % "input squeakr file",
							required("-o", "--output-file") &
							value(ensure_parent_dir_exists, "output-file",
										listopt.output_file) % "output file"
							//option("-h", "--help")  % "show help"
							);

 auto info_mode = (
							command("info").set(selected, mode::info),
							required("-f", "--squeakr-file-file") & value(ensure_file_exists, "squeakr-file",
																									 infoopt.squeakr_file) % "input squeakr file"
							//option("-h", "--help")  % "show help"
							);

  auto cli = (
							(count_mode | query_mode | inner_prod_mode | list_mode |
							 info_mode |
							 command("help").set(selected,mode::help) ),
							option("-v", "--version").call([]{std::cout << "version 0.7\n\n";}).doc("show version")
							);

  assert(count_mode.flags_are_prefix_free());
  assert(query_mode.flags_are_prefix_free());
  assert(inner_prod_mode.flags_are_prefix_free());
  assert(list_mode.flags_are_prefix_free());
  assert(info_mode.flags_are_prefix_free());

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
		std::cout << "\n\nParsing command line failed with exception: " <<
			e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "squeakr");
    return 1;
  }

  if(res) {
		switch(selected) {
			case mode::count:
				try {
					enusure_size_is_specified(countopt);
				} catch (std::exception& e) {
					std::cout << "\n\nParsing command line failed with exception: " <<
						e.what() << "\n";
					std::cout << "\n\n";
					std::cout << make_man_page(cli, "squeakr");
					return 1;
				}
				count_main(countopt);
				break;
			case mode::query: query_main(queryopt);  break;
			case mode::inner_prod: inner_prod_main(innerprodopt);  break;
			case mode::list: list_main(listopt);  break;
			case mode::info: info_main(infoopt);  break;
			case mode::help:  break;
		}
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "count") {
        std::cout << make_man_page(count_mode, "squeakr");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "squeakr");
      } else if (b->arg() == "inner_prod") {
        std::cout << make_man_page(inner_prod_mode, "squeakr");
      } else if (b->arg() == "list") {
        std::cout << make_man_page(list_mode, "squeakr");
      } else if (b->arg() == "info") {
        std::cout << make_man_page(info_mode, "squeakr");
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "squeakr") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "squeakr") << '\n';
    }
  }

  return 0;
}
