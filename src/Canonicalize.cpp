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

int main(int argc, char *argv[]) {
  const int nb_threads      = 10;
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  using std::string;

  uint32_t bucketStringLength;

  string input;
  string output;

  po::options_description generic("canonicalize options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      ("input,i", po::value<string>(&input), "input file")
      ("output,o", po::value<string>(&output), "output file");

  po::variables_map vm;
  try {
      auto orderedOptions = po::command_line_parser(argc, argv).options(generic).run();
      po::store(orderedOptions, vm);

      if (vm.count("help")) {
          auto hstring = R"(
Canonicalize
============
)";
          std::cerr << hstring << "\n";
          std::cerr << generic << "\n";
          std::exit(1);
      }

      po::notify(vm);


      char* input[] = { const_cast<char*>(vm["input"].as<string>().c_str()) };

      stream_manager  streams(input, input + 1, concurrent_file);
      sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);

      std::ostream* of;
      std::ofstream ofile;
      if (output == "-") {
        of = &std::cout;
      } else {
        ofile.open(output);
        of = &ofile;
      }


      bool testFastq{true};
      bool isFastq{false};

      std::vector<jellyfish::header_sequence_qual> reads;
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
              std::replace(j->data[i].seq.begin(), j->data[i].seq.end(), 'N', 'A');
              reads.push_back(j->data[i]);
          }
      }

      std::cerr << "sorting the reads\n";
      tbb::parallel_sort(reads.begin(), reads.end(),
                [](const jellyfish::header_sequence_qual& r1,const jellyfish::header_sequence_qual& r2) -> bool {
                    return r1.seq < r2.seq;
                });


      if (isFastq) {
          for (auto& r : reads) {
              *of << "@" << r.header << "\n";
              *of << r.seq << "\n";
              *of << "+" << "\n";
              *of << r.qual<< "\n";
          }
      } else {
          for (auto& r : reads) {
              *of << ">" << r.header << "\n";
              *of << r.seq << "\n";
          }
      }

      if (output != "-") {
          ofile.close();
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


