#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <functional>
#include <cmath>

#include <boost/pending/property.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

#include <tbb/concurrent_unordered_map.h>

#include "MinceUtils.hpp"
#include "FindPartition.hpp"

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

enum Direction { FORWARD, REVERSE };

struct CountList {
    uint32_t count;
    std::vector<uint32_t> rlist;
};

struct BucketedString {
    std::string str;
    uint8_t offset;
};

using CountMap = tbb::concurrent_unordered_map<uint32_t, uint32_t>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
using MapT = std::map<my_mer, std::vector<std::tuple<std::string, uint8_t>>>;

std::string unpermute(std::string& permS, std::string& key, size_t offset) {
    size_t l = permS.length();
    return permS.substr(l - offset) + key + permS.substr(0, l - offset);
}


std::string permute(std::string& s, size_t offset, size_t kl) {
    return s.substr(offset + kl) + s.substr(0, offset);
}

void decode(std::string& ifname, std::string& ofname) {
    std::ifstream ifile;
    uint8_t rl{0};
    uint8_t kl{0};
    size_t numOnsies{0};

    ifile.open(ifname, std::ios::in | std::ios::binary);
    ifile.read(reinterpret_cast<char*>(&kl), sizeof(kl));
    ifile.read(reinterpret_cast<char*>(&rl), sizeof(rl));
    ifile.read(reinterpret_cast<char*>(&numOnsies), sizeof(numOnsies));

    std::ofstream ofile;
    ofile.open(ofname);

    size_t readLength = rl;
    uint64_t bucketMer{0};
    std::string bucketString('X', kl);
    my_mer bucketKmer;
    my_mer::k(kl);

    size_t effectiveReadLength = std::ceil((readLength) / 4.0);

    size_t i{0};
    char* onsieRead = new char[effectiveReadLength];

    std::cerr << "kl = " << +kl << "\n";
    std::cerr << "readLength = " << readLength << "\n";
    std::cerr << "num onsies = " << numOnsies << "\n";
    // decode all of the onsies
    for (size_t j = 0; j < numOnsies; ++j) {
        ifile.read(reinterpret_cast<char*>(onsieRead), effectiveReadLength);
        std::string recon = mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(onsieRead), readLength);
        ofile << ">" << i << "\n" << recon << std::endl;
    }

    std::cerr << "done with onsies\n";
    delete onsieRead;

    uint32_t maxBucketSize{256};

    effectiveReadLength = std::ceil((readLength - kl) / 4.0);
    std::vector<std::string> reads(maxBucketSize, std::string(effectiveReadLength, 'X'));
    std::vector<uint32_t> offsets(maxBucketSize, 0);
    uint8_t offset{0};
    size_t subBucketSize{0};

    while (ifile.good()) {
        uint8_t bsize{0};
        uint8_t bucketStringLength{0};
        uint32_t numBytes;

        // Read the length of the bucket string (in characters)
        ifile.read(reinterpret_cast<char*>(&bucketStringLength), sizeof(bucketStringLength));
        numBytes = std::ceil(bucketStringLength / 4.0);

        // Read the actual bucket string
        std::vector<uint8_t> coreStrBytes(numBytes);
        ifile.read(reinterpret_cast<char*>(&coreStrBytes[0]), numBytes);
        bucketString = mince::utils::twoBitDecode(&coreStrBytes.front(), bucketStringLength);

        // Read the bucket size
        ifile.read(reinterpret_cast<char*>(&bsize), sizeof(bsize));

        subBucketSize = bsize + 1;

        //std::cerr << "bucket size = " << subBucketSize << "\n";

        effectiveReadLength = std::ceil((readLength - bucketStringLength) / 4.0);

        size_t subBucketCount{0};
        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            //std::cerr << "next read\n";
            //std::cerr << "effectiveReadLength = " << effectiveReadLength << "\n";
            reads[subBucketCount].resize(effectiveReadLength, 'X');
            ifile.read(reinterpret_cast<char*>(&reads[subBucketCount][0]), effectiveReadLength);
        }

        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            ifile.read(reinterpret_cast<char*>(&offsets[subBucketCount]), sizeof(uint8_t));
        }
        // Transform offset deltas into absolute offsets
        uint8_t base = offsets.front();
        for (size_t i = 1; i < subBucketSize; ++i) {
            offsets[i] = base + offsets[i];
            base = offsets[i];
        }

        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            const char* read = reads[subBucketCount].c_str();
            uint32_t offset = offsets[subBucketCount];

                /*
                std::cerr << "HERE6\n";
                std::cerr << "bucket size = " << +subBucketSize << "\n";
                std::cerr << "bucket-string = " << bucketString << "\n";
                std::cerr << "bucket-string length = " << +bucketStringLength << "\n";
                */
            std::string s = mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(read), readLength - bucketStringLength);

                /*
                std::cerr << "[" << s << "]( " << offset << ")\n";
                std::cerr << "i = " << i << "\n";
                std::cerr << "HERE7\n";
                */

            if (offset <= s.length()) {
                //std::cerr << "HERE8\n";
              std::string recon = unpermute(s, bucketString, offset);

              if (recon.length() < readLength) {
                  std::cerr << "recon = " << recon << "\n";
                  std::cerr << "bucket string = " << bucketString << "\n";
                  std::cerr << "offset = " << +offset << "\n";
                  std::cerr << "s = " << s << "\n";
                  std::exit(1);
              }

              ofile << ">" << i << "\n" << recon << std::endl;
                //std::cerr << "recon = " << recon << "\n";
            } else {
                ofile << ">" << i << "\nERROR RECONSTRUCTING THIS READ" << std::endl;
                std::cerr << "s.length() = " << s.length() << "\n";
                std::cerr << "offset = " << +offset << "\n";
                std::cerr << "HERE!!!!!!!!!\n";
                std::exit(1);
            }

            ++i;
            if (i % 100000 == 0) {
                std::cerr << "\r\rwrote read " << i;
            }
        }
    }
    if (ifile.eof()) {
        std::cerr << "reached EOF\n";
    }


    ifile.close();
    ofile.close();
}

template <typename MerT>
void bucketReads(sequence_parser* parser, std::atomic<uint64_t>* total,
        std::atomic<uint64_t>& totReads,
        MapT& buckets,
        CountMap& countMap,
        uint32_t bucketStringLength,
        bool noRC) {

    uint64_t count = 0;
    unsigned int kl{bucketStringLength};
    unsigned int cmlen{0};

    while(true) {
        sequence_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;          // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
            if (totReads++ % 1000000 == 0) {
                std::cerr << "\r\rprocessed " << totReads << " reads";
            }
            std::string s(j->data[i].seq);

            count += s.size();        // Add up the size of the sequence

            cmlen = 0;
            MerT mer(kl);
            MerT rmer(kl);
            MerT minmer(kl);
            mer.polyT();
            rmer.polyT();
            minmer.polyT();

            size_t readLength{s.length()};
            std::map<uint32_t, std::tuple<uint8_t, Direction>> merOffsetMap;

            size_t offset{0};
            size_t roffset{readLength};
            size_t firstOffset{offset};
            while (offset < s.size()) {
                // Get the code for the next base
                int c = jellyfish::mer_dna::code(s[offset]);

                // If it's not a valid DNA code
                if (jellyfish::mer_dna::not_dna(c)) {
                    // Switch it to an 'A' in the mer
                    c = jellyfish::mer_dna::code('A');
                    // Mark it as 'N' in the string (will be changed to 'A' later)
                    s[offset] = 'N';
                    mer.shift_left(c);
                    rmer.shift_right(c);
                } else { // Otherwise, base is legit
                    mer.shift_left(c);
                    rmer.shift_right(jellyfish::mer_dna::complement(c));
                }

                ++offset;
                --roffset;
                // If we've read a full k-mer
                if (++cmlen >= kl) {
                    cmlen = kl;

                    auto key = mer.get_bits(0, 2*kl);
                    auto rkey = rmer.get_bits(0, 2*kl);

                    // If mer hasn't been seen yet, then record it
                    auto it = merOffsetMap.find(key);
                    if (it == merOffsetMap.end()) {
                        merOffsetMap[key] = std::make_tuple(offset - kl, Direction::FORWARD);
                    }

                    if (!noRC) {
                        // If rmer hasn't been seen yet, then record it
                        it = merOffsetMap.find(rkey);
                        if (it == merOffsetMap.end()) {
                            merOffsetMap[rkey] = std::make_tuple(roffset, Direction::REVERSE);
                        }
                    }

                }
            }

            // Decide the heaviest bucket --- this is where the read will be placed
            uint32_t maxBucketKey{merOffsetMap.begin()->first};
            uint32_t maxBucketValue{0};
            Direction maxBucketDirection{std::get<1>(merOffsetMap.begin()->second)};
            std::atomic<uint32_t> count;
            // Iterate over all k-mers and reverse complement k-mers
            for (auto& kv : merOffsetMap) {
                auto cit = countMap.find(kv.first);
                // If the bucket for this key is the heaviest so far
                // then it becomes the new bucket key.
                if (cit != countMap.end() and cit->second > maxBucketValue) {
                    maxBucketValue = cit->second;
                    maxBucketDirection = std::get<1>(kv.second);
                    maxBucketKey = kv.first;
                }
            }

            // If the heaviest bucket contains the RC mer,
            // then RC the read
            if (maxBucketDirection == Direction::REVERSE) {
                mince::utils::reverseComplement(s);
            }
            // Increment the count of the assigned bucket
            countMap[maxBucketKey]++;

            // Encode the read accounting for the shared string
            minmer.set_bits(0, 2*kl, maxBucketKey);
            // The first appearance of the key in the read
            firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);

            std::replace(s.begin(), s.end(), 'N', 'A');
            // orig: x . key . y
            // new : y . x
            std::string reord = permute(s, firstOffset, kl);

            std::string kstr = minmer.to_str();

#ifdef DEBUG
            if ( s != unpermute(reord, kstr, firstOffset)) {
                std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
            }
            if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "A string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
            }
#endif
            // Put the reordered string, offset tuple in the bucket
            buckets[minmer].push_back(make_tuple(reord, firstOffset));
        } // for all of the reads in this job
    } // for all jobs of this thread
}

void reassignOnsies(
        std::vector<std::tuple<my_mer, std::string>>& onsies,
        CountMap& countMap,
        MapT& buckets,
        uint32_t bucketStringLength,
        bool noRC
        ) {

    for (auto& kr : onsies) {
        std::string s = std::get<1>(kr);
        size_t readLength{s.length()};
        size_t offset{0};
        size_t roffset{readLength};
        size_t cmlen{0};
        size_t firstOffset{offset};
        size_t kl{bucketStringLength};
        my_mer mer;
        my_mer rmer;
        my_mer minmer;
        my_mer::k(kl);
        mer.polyT();
        rmer.polyT();

        std::map<uint32_t, std::tuple<uint8_t, Direction>> merOffsetMap;
        while (offset < s.size()) {

            int ca = jellyfish::mer_dna::code(s[offset]);
            /*
            if (jellyfish::mer_dna::not_dna(ca)) {
                ca = static_cast<int>('A');
                s[offset] = 'A';
            }
            */

            mer.shift_left(ca);
            rmer.shift_right(jellyfish::mer_dna::complement(ca));

            ++offset;
            --roffset;
            if (++cmlen >= kl) {
                cmlen = kl;
                auto key = mer.get_bits(0, 2*kl);
                auto rkey = rmer.get_bits(0, 2*kl);

                auto it = merOffsetMap.find(key);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[key] = std::make_tuple(offset - kl, Direction::FORWARD);
                    if (s.substr(offset - kl, kl) != mer.to_str()) {
                        std::cerr << "putting " << mer.to_str() << " into the map, but I really have " << s.substr(offset - kl, kl) << "\n";
                    }
                }

                if (!noRC) {

                    // If rmer hasn't been seen yet, then record it
                    it = merOffsetMap.find(rkey);
                    if (it == merOffsetMap.end()) {
                        merOffsetMap[rkey] = std::make_tuple(roffset, Direction::REVERSE);
                    }
                }
            }
        }

        uint32_t maxBucketKey{merOffsetMap.begin()->first};
        uint32_t maxBucketValue{0};
        Direction maxBucketDirection{std::get<1>(merOffsetMap.begin()->second)};
        std::atomic<uint32_t> count;
        for (auto& kv : merOffsetMap) {
            auto cit = countMap.find(kv.first);
            if (cit != countMap.end() and cit->second > maxBucketValue) {
                maxBucketValue = cit->second;
                maxBucketDirection = std::get<1>(kv.second);
                maxBucketKey = kv.first;
            }
        }

        if (maxBucketDirection == Direction::REVERSE) {
            mince::utils::reverseComplement(s);
        }

        countMap[maxBucketKey]++;
        minmer.set_bits(0, 2*kl, maxBucketKey);
        firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);

        std::string reord = permute(s, firstOffset, kl);

        std::string kstr = minmer.to_str();

#ifdef DEBUG
        if ( s != unpermute(reord, kstr, firstOffset)) {
            std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
        }
#endif
        buckets[minmer].push_back(std::make_tuple(reord, firstOffset));

    }

}

int main(int argc, char *argv[]) {
  const int nb_threads      = 20;
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  using std::string;

  bool doDecode{false};
  bool doEncode{false};
  bool noRC{false};

  uint32_t bucketStringLength;

  po::options_description generic("mince options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      ("input,i", po::value<string>(), "input file; FASTA/Q if encoding, MINCE if decoding")
      ("output,o", po::value<string>(), "output base file")
      ("blength,b", po::value<uint32_t>(&bucketStringLength)->default_value(15) , "length of the bucket string [1,16]")
      ("norc,n", po::bool_switch(&noRC)->default_value(false), "don't allow reverse complementing when encoding / bucketing")
      ("encode,e", po::bool_switch(&doEncode)->default_value(false), "encode the input file")
      ("decode,d", po::bool_switch(&doDecode)->default_value(false), "decode the input file");

  po::variables_map vm;
  try {
      auto orderedOptions = po::command_line_parser(argc, argv).options(generic).run();
      po::store(orderedOptions, vm);

      if (vm.count("help")) {
          auto hstring = R"(
Mince
=====
)";
          std::cerr << hstring << "\n";
          std::cerr << generic << "\n";
          std::exit(1);
      }

      po::notify(vm);

      for (auto& opt : orderedOptions.options) {
          std::cerr << "[ " << opt.string_key << "] => {";
          for (auto& val : opt.value) {
              std::cerr << " " << val;
          }
          std::cerr << " }\n";
      }

      if (bucketStringLength < 1 or bucketStringLength > 16) {
          std::stringstream errstr;
          errstr << "The length of the bucket string must be between 1 and 16 (inclusive)";
          throw std::invalid_argument(errstr.str());
      }

      if (doEncode and doDecode) {
          std::stringstream errstr;
          errstr << "You cannot pass both the encode and decode flag in a single run";
          throw std::invalid_argument(errstr.str());
      }
      if (!doEncode and !doDecode) {
          std::stringstream errstr;
          errstr << "You must select either encoding or decoding by passing the -e or -d flag";
          throw std::invalid_argument(errstr.str());
      }

      if (doDecode) {
          std::string input = vm["input"].as<string>();
          std::string output = vm["output"].as<string>();
          decode(input, output);
          std::exit(0);
      }

      // Currently only single-end or pre-merged paired-end reads
      // No reason we can't do the merging on the fly and accept paired-end reads directly
      char* input[] = { const_cast<char*>(vm["input"].as<string>().c_str()) };

      stream_manager  streams(input, input + 1, concurrent_file);
      sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);

      std::string outfname = vm["output"].as<string>();

      std::vector<std::thread> threads;
      std::atomic<uint64_t> total(0);
      std::atomic<uint64_t> totReads{0};
      std::atomic<uint64_t> distinct{0};
      using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

      std::ofstream outfileSeqs;
      outfileSeqs.open(outfname+".seqs", std::ios::out | std::ios::binary);
      std::ofstream outfileOffsets;
      outfileOffsets.open(outfname+".offs", std::ios::out | std::ios::binary);
      unsigned int kl{bucketStringLength};
      my_mer::k(kl);

      CountMap countMap;
      std::vector<MapT> maps(nb_threads);
      // Spawn of the read parsing threads
      for(int i = 0; i < nb_threads; ++i) {
          threads.push_back(std::thread(bucketReads<my_mer>, &parser, &total,
                                        std::ref(totReads), std::ref(maps[i]),
                                        std::ref(countMap), bucketStringLength, noRC));
      }

      // Join all of the read-parsing threads
      for(int i = 0; i < nb_threads; ++i) {
          threads[i].join();
      }

      // The union of keys from all maps
      std::set<my_mer> keys;
      for (auto& m : maps) {
          for (auto& kv : m ) {
              keys.insert(kv.first);
          }
      }

      // The length of the reads in this file
      uint8_t readLength = std::get<0>(maps[0].begin()->second.front()).length() + kl;
      std::cerr << "\n\nread length = " << +readLength << "\n";

      size_t totOutput{0};
      size_t totBuckets{0};
      size_t bnum{0};
      size_t totNumReads{0};

      // The total # of buckets
      std::cerr << "Num buckets = " << keys.size() << "\n";

      std::vector<std::tuple<my_mer, std::string>> onsies;
      // For all buckets
      for (auto& k : keys) {
          std::vector<std::tuple<my_mer, std::string, size_t>> combined;
          bool isOnsie{true};

          // For all maps
          for (auto& m : maps) {
              // If this map contains this bucket, and this bucket is of size one
              // in this map, then the contained read is a potential onsie
              auto kit = m.find(k);
              if (kit != m.end()) {
                  if (kit->second.size() == 1) {
                      combined.push_back(std::forward_as_tuple(k,
                                  std::get<0>(kit->second.front()),
                                  std::get<1>(kit->second.front())
                                  ));
                  } else {
                      isOnsie = false;
                  }
                  totNumReads += kit->second.size();
              }
          }

          // If, over all maps, this read appears in a singleton bucket
          if (isOnsie and combined.size() == 1) {
              // The string on which this read was originally bucketed
              my_mer bucketMer = std::get<0>(combined[0]);
              auto bucketString = bucketMer.to_str();

              // The split / swapped string
              std::string s(std::get<1>(combined[0]));
              // Recover the original string
              size_t offset = std::get<2>(combined[0]);
              std::string recon = unpermute(s, bucketString, offset);

              // Push back the original key and the reconstructed string
              onsies.push_back(std::forward_as_tuple(k, recon));
          }
      }

      // For all onsies
      for (auto& kr : onsies) {
          // erase this guy from the map and the set of keys
          for (auto& m : maps) {
              auto key = std::get<0>(kr);
              auto kit = m.find(key);
              if (kit != m.end()) {
                  if (kit->second.size() > 1) {
                      std::cerr << "erasing non-singleton bucket (" << kit->second.size() << ")!\n";
                  }
                  m.erase(kit);
                  keys.erase(key);
                  countMap[key.get_bits(0, 2*kl)]--;
              }
          }
      }

      std::cerr << "|keys| = " << keys.size() << "\n";
      std::cerr << "there were " << totNumReads << " original reads\n";
      std::cerr << "there were " << onsies.size() << " original onsies\n";

      // Re-assign the onsies to see if they can be placed in any non-empty
      // buckets that didn't exist when they were first considered.
      MapT collatedMap;
      reassignOnsies(onsies, countMap, collatedMap, bucketStringLength, noRC);

      size_t newOnsies{0};
      size_t newTotal{0};
      std::set<my_mer> onsieKeys;
      std::vector<std::tuple<my_mer, std::string>> onsieVec;
      for (auto& kv : collatedMap) {
          keys.insert(kv.first);
          if (kv.second.size() == 1) {
              ++newOnsies;
              onsieKeys.insert(kv.first);
              std::string bucketString = kv.first.to_str();
              auto& so = kv.second.front();
              std::string recon = unpermute(std::get<0>(so), bucketString, std::get<1>(so));
              onsieVec.emplace_back(std::forward_as_tuple(kv.first, recon));
          }
          newTotal += kv.second.size();
      }

      std::cerr << "|keys| = " << keys.size() << "\n";
      std::cerr << "now there are " << newOnsies << " onsies\n";
      std::cerr << "rebucketed " << newTotal << " reads\n";


      // The length of the string on which we'll bucket
      uint8_t bslen = bucketStringLength;
      outfileSeqs.write(reinterpret_cast<char*>(&bslen), sizeof(bslen));
      // The length of the reads
      outfileSeqs.write(reinterpret_cast<char*>(&readLength), sizeof(readLength));

      std::sort(onsieVec.begin(), onsieVec.end(),
              [](const std::tuple<my_mer, std::string>& o1, const std::tuple<my_mer, std::string>& o2) -> bool {
                return std::get<1>(o1) < std::get<1>(o2);
              });

      // First, write out all remaining onsies
      outfileSeqs.write(reinterpret_cast<char*>(&newOnsies), sizeof(newOnsies));
      for (auto& ks : onsieVec) {
          auto& k = std::get<0>(ks);
          if (collatedMap[k].size() != 1) { std::cerr << "onsie was not a onsie!!!\n"; std::exit(1); }
          auto& recon = std::get<1>(ks);
          auto bytes = mince::utils::twoBitEncode(recon);
          outfileSeqs.write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
          totOutput++;
          collatedMap.erase(k);
      }

      // Push back the remaining onsies map
      maps.push_back(collatedMap);

      std::cerr << "wrote out onsies\n";

      // Go through all of the buckets
      for (auto& k : keys ) {
          // Collect all of the reads for this bucket
          std::vector<std::tuple<std::string, size_t>> reads;
          for (auto& m : maps) {
              auto kit = m.find(k);
              if (kit != m.end()) {
                  reads.insert(
                          reads.end(),
                          std::make_move_iterator(kit->second.begin()),
                          std::make_move_iterator(kit->second.end())
                          );

              }
          }

          // Sort the reads; first by the offset of the bucket string then lexicographically
          std::sort(reads.begin(), reads.end(),
                    [&](const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
                        int cmp = std::get<1>(a) - std::get<1>(b);
                        if (cmp == 0) {
                            return std::get<0>(a) < std::get<0>(b);
                        } else {
                            return cmp < 0;
                        }
                  });

          // Report bucket progress
          ++bnum;
          if (bnum % 100 == 0) {
              std::cerr << "\r\rprocessed " << bnum << "buckets";
          }

          // The number of sub-buckets into which this bucket will be split
          size_t subBuckets = std::ceil(reads.size() / 256.0);

          // Real optimal partition
          //std::vector<mince::BucketInfo> optimalPartition = mince::findPartition(reads, bucketStringLength, readLength, 255);

          // Dummy partition
          std::vector<mince::BucketInfo> optimalPartition;
          size_t N = reads.size();
          for (size_t i = 0; i < N; i += 256) {
              uint32_t m = std::min(N-1, i+255);
              optimalPartition.push_back(mince::BucketInfo(1, i, m, -1, 0));
          }

          // Write out all of the sub-buckets
          size_t nextReadIdx = 0;
          for (size_t i = 0; i < optimalPartition.size(); ++i) {
              auto& subBucket = optimalPartition[i];
              // current bucket is [nextRead, optimalPartition[i]]
              size_t nextReadIdx = subBucket.startIdx;
              size_t lastReadIdx = subBucket.stopIdx;
              size_t sbsize = subBucket.size();
              uint32_t lcpLen = subBucket.prefLen;

              if (nextReadIdx > lastReadIdx) {
                  std::cerr << "next read index = " << nextReadIdx << " > "
                            << "last read index = " << lastReadIdx << "\n";
                  std::exit(1);
                  continue;
              }

              std::string& firstRead = std::get<0>(reads[nextReadIdx]);
              std::string& lastRead = std::get<0>(reads[lastReadIdx]);

              size_t offset1 = std::get<1>(reads[nextReadIdx]);
              size_t offset2 = std::get<1>(reads[lastReadIdx]);

              if (nextReadIdx == lastReadIdx) { lcpLen = 0; }
              std::string keyStr = k.to_str() + std::get<0>(reads[nextReadIdx]).substr(0, lcpLen);
              auto keyStrLen = keyStr.length();

              if (offset1 > firstRead.length()) {
                  std::cerr << "offset1 = " << offset1 << " > "
                            << "read length = " << firstRead.length() << "\n";
                  std::exit(1);
              }
              if (offset2 > lastRead.length()) {
                  std::cerr << "offset2 = " << offset2 << " > "
                            << "read length = " << lastRead.length() << "\n";
                  std::exit(1);
              }

              std::vector<uint8_t> keyBytes = mince::utils::twoBitEncode(keyStr);
              uint8_t kstrLen = keyStr.length();
              // If this is the first sub-bucket to use this string, then write out the core string
              if (i == 0) {
                  // The size of the bucket string
                  outfileSeqs.write(reinterpret_cast<char*>(&kstrLen), sizeof(kstrLen));
                  // The bucket string
                  outfileSeqs.write(reinterpret_cast<char*>(&keyBytes[0]), sizeof(uint8_t) * keyBytes.size());
              } else { // otherwise, write out the sentinel value, 0, and assume this bucket uses the same core
                  uint8_t zero{0};
                  outfileSeqs.write(reinterpret_cast<char*>(&zero), sizeof(zero));
              }
              // The size of this sub-bucket
              uint8_t subBucketSize = static_cast<uint8_t>(sbsize - 1);
              outfileSeqs.write(reinterpret_cast<char*>(&subBucketSize), sizeof(subBucketSize));

              size_t encReadLength = readLength - keyStr.length();
              size_t extraBasesPerRead = encReadLength % 4;
              std::stringstream sstream;
              extraBasesPerRead = 0;
              // Write out the reads of this sub-bucket
              for (size_t j = nextReadIdx; j <= lastReadIdx; ++j) {
                  auto& r = reads[j];
                  std::string rstr = std::get<0>(r);
                  uint8_t offset = std::get<1>(r);

                  try {
                      std::string rsub = rstr.substr(lcpLen, encReadLength - extraBasesPerRead);
                      if (offset > rsub.size()) {
                          std::cerr << "bucket[" << nextReadIdx << ", " << lastReadIdx << "]\n";
                          std::cerr << "lcpLen = " << lcpLen << "\n";
                          std::cerr << "encountered an offset " << +offset << " that is greater than the string length " << rsub.length() << "\n";
                      }

                      std::vector<uint8_t> bytes = mince::utils::twoBitEncode(rsub);
                      outfileSeqs.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
                      if (extraBasesPerRead > 0) {
                          sstream << rstr.substr(encReadLength - extraBasesPerRead);
                      }

                  } catch (std::exception& e ) {
                      std::cerr << "EXCEPTION: " << e.what() << rstr << "\n";
                  }
              }

              // Write out the offsets (delta-encoded) of the bucket strings within the reads
              std::vector<uint8_t> offsets;
              offsets.reserve(sbsize);
              for (size_t j = nextReadIdx; j <= lastReadIdx; ++j) {
                  auto& r = reads[j];
                  uint8_t l = static_cast<uint8_t>(std::get<1>(r));

                  if (l > std::get<0>(r).length()) {
                    std::cerr << "offset is " << l << ", which is greater than "
                              << "text length " << std::get<0>(r).length() << "\n";
                  }
                  offsets.push_back(l);
                  totOutput++;
              }

              std::vector<uint8_t> deltas;
              deltas.reserve(sbsize);
              deltas.push_back(offsets[0]);
              for (size_t j = 1; j < sbsize; ++j) {
                  int diff = (offsets[j] - offsets[j-1]);
                  if (diff < 0) {
                      std::cerr << "offsets[" << j << "] = " << +offsets[j]
                                << "offsets[" << j-1 << "] = " << +offsets[j-1] << "\n";
                      std::exit(1);
                  }
                  deltas.push_back(static_cast<uint8_t>(diff));
              }
              outfileOffsets.write(reinterpret_cast<char*>(&deltas[0]), deltas.size());
          }

          // Increment the total bucket count
          totBuckets++;
      }

      outfileSeqs.close();
      outfileOffsets.close();
      std::cerr << "Wrote: " << totOutput << " reads, in " << totBuckets << " buckets\n";
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


