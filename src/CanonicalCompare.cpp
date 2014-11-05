#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <functional>
#include <cmath>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

#include <tbb/parallel_sort.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "MinceUtils.hpp"

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

std::vector<std::string> parseReadFile(const std::string& input) {
	const int nb_threads      = 10;
	const int concurrent_file = 1;   // Number of files to read simultaneously
	const int max_read_group  = 100; // Number of reads in each "job" group

	using std::string;
	std::vector<string> reads;
	char* inputs[] = { const_cast<char*>(input.c_str()) };
	stream_manager  streams(inputs, inputs + 1, concurrent_file);
	sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);

	bool testFastq{true};
	bool isFastq{false};

	size_t totReads{0};
	while(true) {
		sequence_parser::job j(parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
		if(j.is_empty()) break;          // If we have nothing, quit the loop

		for(size_t i = 0; i < j->nb_filled; ++i) { // For all the reads we have
			if (totReads++ % 1000000 == 0) {
				std::cerr << "\r\rprocessed " << totReads << " reads";
			}
			if (testFastq) {
				if (j->data[i].seq.length() == j->data[i].qual.length()) {
					isFastq = true;
				}
				testFastq = false;
			}
			string s(j->data[i].seq);
			string rcs(s);
			mince::utils::reverseComplement(rcs);

			j->data[i].seq = (rcs < s) ? rcs : s;
			//std::replace(j->data[i].seq.begin(), j->data[i].seq.end(), 'N', 'A');
			reads.push_back(j->data[i].seq);
		}
	}

	std::cerr << "sorting the reads\n";
	tbb::parallel_sort(reads.begin(), reads.end(),
			[](const string& r1, const string& r2) -> bool {
			return r1 < r2;
			});
	return reads;
}


int main(int argc, char *argv[]) {
  const int nb_threads      = 10;
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  using std::string;

  uint32_t bucketStringLength;

  string input1;
  string input2;

  po::options_description generic("canonicalize options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      ("input1,1", po::value<string>(&input1), "input file 1")
      ("input2,2", po::value<string>(&input2), "input file 2");

  po::variables_map vm;
  try {
      auto orderedOptions = po::command_line_parser(argc, argv).options(generic).run();
      po::store(orderedOptions, vm);

      if (vm.count("help")) {
          auto hstring = R"(
Compare Canonical
=================
)";
          std::cerr << hstring << "\n";
          std::cerr << generic << "\n";
          std::exit(1);
      }

      po::notify(vm);

      std::vector<string> reads1 = parseReadFile(input1);
      std::vector<string> reads2 = parseReadFile(input2);
      if (reads1.size() != reads2.size()) {
	std::cerr << "files have a different number of reads\n";
	std::cerr << "file 1 has " << reads1.size() << " reads\n";
	std::cerr << "file 2 has " << reads2.size() << " reads\n";
	std::exit(1);
      }

      for (size_t i = 0; i < reads1.size(); ++i) {
	      auto& r1 = reads1[i];
	      auto& r2 = reads2[i];
	      if (r1 != r2) {
		      std::cerr << "files first differed at read " << i << "\n";
		      std::cerr << "file 1 had: " << r1 << "\n";
		      std::cerr << "file 2 had: " << r2 << "\n";
		      std::exit(1);
	      }
      }
   } // end try
  catch (po::error &e) {
      std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
      std::exit(1);
  } catch (std::exception& e) {
      std::cerr << "Exception : [" << e.what() << "]\n";
      std::cerr << argv[0] << " was invoked improperly.\n";
      std::cerr << "For usage information, try " << argv[0] << " --help\nExiting.\n";
  }

  return 0;
}


