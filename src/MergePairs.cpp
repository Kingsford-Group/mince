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
#include "PairSequenceParser.hpp"

using pair_parser = pair_sequence_parser<char**>;

int main(int argc, char *argv[]) {
    const int nb_threads      = 10;
    const int concurrent_file = 1;   // Number of files to read simultaneously
    const int max_read_group  = 100; // Number of reads in each "job" group

    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;
    using std::string;

    uint32_t bucketStringLength;

    string output;

    po::options_description generic("merge-pairs options");
    generic.add_options()
        ("version,v", "print version information")
        ("help,h", "print help message")
        ("libtype,l", po::value<string>(), "library format string")
        ("mates1,1", po::value<std::vector<string>>(), "mate file 1")
        ("mates2,2", po::value<std::vector<string>>(), "mate file 2")
        ("output,o", po::value<string>(&output), "output file");

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc, argv).options(generic).run();
        po::store(orderedOptions, vm);

        if (vm.count("help")) {
            auto hstring = R"(
                Merge Pairs
                ============
                )";
            std::cerr << hstring << "\n";
            std::cerr << generic << "\n";
            std::exit(1);
        }

        po::notify(vm);

        auto readLibraries = mince::utils::extractReadLibraries(orderedOptions);

        for (auto& rl : readLibraries) {
            if (rl.format().type == ReadType::SINGLE_END) {
                std::cerr << "Cannot merge single end library; skipping\n";
                continue;
            }
            //std::cerr << "Merging read library " << "\n";
            std::cerr << "Merging ";
            std::cerr << const_cast<LibraryFormat&>(rl.format());
            std::cerr << " format library\n";
            std::vector<const char*> inputs;
            auto& m1 = rl.mates1();
            auto& m2 = rl.mates2();
            for (size_t i = 0; i < m1.size(); ++i) {
                inputs.push_back(m1[i].c_str());
                std::cerr << "merging " << inputs.back() << " with ";
                inputs.push_back(m2[i].c_str());
                std::cerr << inputs.back() << "\n";
            }

            std::function<void(std::pair<header_sequence_qual, header_sequence_qual>&)> mergeFn;
            if (rl.format().orientation == ReadOrientation::TOWARD) {
                mergeFn = [](std::pair<header_sequence_qual, header_sequence_qual>& rp) -> void {
                    mince::utils::reverseComplement(rp.second.seq);
                    std::reverse(rp.second.qual.begin(), rp.second.qual.end());
                };
            } else if (rl.format().orientation == ReadOrientation::AWAY) {
                mergeFn = [](std::pair<header_sequence_qual, header_sequence_qual>& rp) -> void {
                    mince::utils::reverseComplement(rp.second.seq);
                    std::reverse(rp.second.qual.begin(), rp.second.qual.end());
                };
            } else if (rl.format().orientation == ReadOrientation::SAME) {
                // If the reads are oriented in the same direction, they're from the same strand --- do nothing
                mergeFn = [](std::pair<header_sequence_qual, header_sequence_qual>& rp) -> void {
                    return;
                };
            }


            char** start = const_cast<char**>(&(*inputs.begin()));
            char** stop = const_cast<char**>(&(*inputs.end()));
            pair_parser parser(4 * nb_threads, max_read_group, concurrent_file,
                               start, stop);
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

            size_t totReads{0};
            while(true) {
                pair_parser::job j(parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
                if(j.is_empty()) break;          // If we have nothing, quit the loop

                for(size_t i = 0; i < j->nb_filled; ++i) { // For all the reads we have
                    if (totReads++ % 1000000 == 0) {
                        std::cerr << "\r\rprocessed " << totReads << " reads";
                    }
                    if (testFastq) {
                        if (j->data[i].first.seq.length() == j->data[i].first.qual.length()) {
                            isFastq = true;
                        }
                        testFastq = false;
                    }

                    mergeFn(j->data[i]);

                    if (isFastq) {
                        *of << "@" << j->data[i].first.header << "\n";
                        *of << j->data[i].first.seq << j->data[i].second.seq << "\n";
                        *of << "+\n";
                        *of << j->data[i].first.qual << j->data[i].second.qual << "\n";
                    } else {
                        *of << ">" << j->data[i].first.header << "\n";
                        *of << j->data[i].first.seq << j->data[i].second.seq << "\n";
                    }

                } // end of job
            } // end of while


            if (output != "-") {
                ofile.close();
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


