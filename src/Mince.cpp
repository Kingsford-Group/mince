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

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

enum Direction { FORWARD, REVERSE };

struct CountList {
    uint32_t count;
    std::vector<uint32_t> rlist;
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
    ifile.open(ifname, std::ios::in | std::ios::binary);
    ifile.read(reinterpret_cast<char*>(&kl), sizeof(kl));
    ifile.read(reinterpret_cast<char*>(&rl), sizeof(rl));
    std::ofstream ofile;
    ofile.open(ofname);

    size_t readLength = rl;

    uint64_t bucketMer{0};
    std::string bucketString('X', kl);
    my_mer bucketKmer;
    my_mer::k(kl);

    uint32_t maxBucketSize{256};

    size_t effectiveReadLength = std::ceil((readLength - kl) / 4.0);
    std::vector<std::string> reads(maxBucketSize, std::string(effectiveReadLength, 'X'));
    std::vector<uint32_t> offsets(maxBucketSize, 0);
    uint8_t offset{0};
    size_t i{0};
    size_t subBucketSize{0};
    while (ifile.good()) {
        ifile.read(reinterpret_cast<char*>(&bucketMer), sizeof(bucketMer));
        uint8_t bsize{0};
        ifile.read(reinterpret_cast<char*>(&bsize), sizeof(bsize));
        subBucketSize = bsize + 1;
        bucketKmer.set_bits(0, 2*kl, bucketMer);
        bucketString = bucketKmer.to_str();

        size_t subBucketCount{0};
        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            //std::cerr << "next read\n";
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

            std::string s(mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(read), readLength - kl));
            std::string recon = unpermute(s, bucketString, offset);
            ofile << ">" << i << "\n" << recon << std::endl;
            /*
            if (offset == 0) {
                std::string recon = bucketString + s;
                ofile << ">" << i << "\n" << recon << std::endl;
            }
            else if (offset < s.size()) {
                std::string recon = s.substr(newOffset) + bucketString + s.substr(0, newOffset);
                //std::cerr << "A recon.size() " << recon.size() << "\n";
                ofile << ">" << i << "\n" << recon << std::endl;
            } else {
                std::string recon = s + bucketString;
                //std::cerr << "B recon.size() " << recon.size() << "\n";
                ofile << ">" << i << "\n" << recon << std::endl;
            }
            */
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
        uint32_t bucketStringLength) {

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

                    // If rmer hasn't been seen yet, then record it
                    it = merOffsetMap.find(rkey);
                    if (it == merOffsetMap.end()) {
                        merOffsetMap[rkey] = std::make_tuple(roffset, Direction::REVERSE);
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
            if ( s != unpermute(reord, kstr, firstOffset)) {
                std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
            }
            if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "A string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
            }
            // Put the reordered string, offset tuple in the bucket
            buckets[minmer].push_back(make_tuple(reord, firstOffset));
            /*
            // If the first appearance isn't at the end if the read
            if (firstOffset + kl < s.size()) {
                // orig: x . key . y
                // new : y . x
                std::string reord = permute(s, firstOffset, kl);

                std::string kstr = minmer.to_str();
                if ( s != unpermute(reord, kstr, firstOffset)) {
                    std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
                }
                if (s.substr(firstOffset, kl) != minmer.to_str()) {
                    std::cerr << "A string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                    std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
                }
                // Put the reordered string, offset tuple in the bucket
                buckets[minmer].push_back(make_tuple(reord, firstOffset));
            } else { // if the first appearance is at the end of the string, then the 'y' substring is empty
                std::string reord = permute(s, firstOffset, kl);

                std::string kstr = minmer.to_str();
                if ( s != unpermute(reord, kstr, firstOffset)) {
                    std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
                }

                buckets[minmer].push_back(make_tuple(reord, firstOffset));
                if (s.substr(firstOffset, kl) != minmer.to_str()) {
                    std::cerr << "B string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                    std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
                }
            }
            */
        } // for all of the reads in this job
    } // for all jobs of this thread
}

void reassignOnsies(
        std::vector<std::tuple<my_mer, std::string>>& onsies,
        CountMap& countMap,
        MapT& buckets,
        uint32_t bucketStringLength
        ) {

    for (auto& kr : onsies) {
        std::string s = std::get<1>(kr);
        size_t readLength{s.length()};
        size_t offset{0};
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
                /*
                it = merOffsetMap.find(rkey);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[rkey] = std::make_tuple(readLength - offset - 1, Direction::REVERSE);
                }
                */
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
        if ( s != unpermute(reord, kstr, firstOffset)) {
            std::cerr << "ERROR!!! unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
        }

        buckets[minmer].push_back(std::make_tuple(reord, firstOffset));

        /*
        if (firstOffset + kl < s.size()) {
            std::string reord(s.substr(firstOffset + kl) + s.substr(0, firstOffset));
            if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "A string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
            }

            buckets[minmer].push_back(std::make_tuple(reord, firstOffset));
        } else {
            buckets[minmer].push_back(std::make_tuple(s.substr(0, firstOffset), firstOffset));
            if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "B string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
            }
        }
        */


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

  uint32_t bucketStringLength;

  po::options_description generic("mince options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      ("input,i", po::value<string>(), "input file; FASTA/Q if encoding, MINCE if decoding")
      ("output,o", po::value<string>(), "output file")
      ("blength,b", po::value<uint32_t>(&bucketStringLength)->default_value(15) , "length of the bucket string [1,16]")
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

      char* input[] = { const_cast<char*>(vm["input"].as<string>().c_str()) };

      stream_manager  streams(input, input + 1, concurrent_file);
      sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);

      std::string outfname = vm["output"].as<string>();

      std::vector<std::thread> threads;
      std::atomic<uint64_t> total(0);
      std::atomic<uint64_t> totReads{0};
      std::atomic<uint64_t> distinct{0};
      using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

      std::ofstream outfile;
      outfile.open(outfname, std::ios::out | std::ios::binary);
      unsigned int kl{bucketStringLength};
      my_mer::k(kl);

      CountMap countMap;
      std::vector<MapT> maps(nb_threads);
      // Spawn of the read parsing threads
      for(int i = 0; i < nb_threads; ++i) {
          threads.push_back(std::thread(bucketReads<my_mer>, &parser, &total, std::ref(totReads), std::ref(maps[i]), std::ref(countMap), bucketStringLength));
      }

      // Join all of the read-parsing threads
      for(int i = 0; i < nb_threads; ++i) {
          threads[i].join();
      }

      std::set<my_mer> keys;
      for (auto& m : maps) {
          for (auto& kv : m ) {
              keys.insert(kv.first);
          }
      }

      uint8_t readLength = std::get<0>(maps[0].begin()->second.front()).length() + kl;
      std::cerr << "\n\nread length = " << +readLength << "\n";
      //outfile << +readLength << "\n";

      uint8_t bslen = bucketStringLength;
      outfile.write(reinterpret_cast<char*>(&bslen), sizeof(bslen));
      outfile.write(reinterpret_cast<char*>(&readLength), sizeof(readLength));

      size_t totOutput{0};
      size_t totBuckets{0};
      size_t bnum{0};
      size_t totNumReads{0};
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
              // Recover the original string
              size_t offset = std::get<2>(combined[0]);
              // The string on which this read was originally bucketed
              my_mer bucketMer = std::get<0>(combined[0]);
              auto bucketString = bucketMer.to_str();

              // The split / swapped string
              std::string s(std::get<1>(combined[0]));
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

      MapT collatedMap;
      reassignOnsies(onsies, countMap, collatedMap, bucketStringLength);
      maps.push_back(collatedMap);

      size_t newOnsies{0};
      size_t newTotal{0};
      for (auto& kv : collatedMap) {
          keys.insert(kv.first);
          if (kv.second.size() == 1) { ++newOnsies; }
          newTotal += kv.second.size();
      }

      std::cerr << "|keys| = " << keys.size() << "\n";
      std::cerr << "now there are " << newOnsies << " onsies\n";
      std::cerr << "rebucketed " << newTotal << " reads\n";

      for (auto& k : keys ) {
          //for (auto& kv : buckets) {
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

          std::sort(reads.begin(), reads.end(),
                  [&](const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
                  int cmp = std::get<1>(a) - std::get<1>(b);
                  if (cmp == 0) {
                  return std::get<0>(a) < std::get<0>(b);
                  } else {
                  return cmp < 0;
                  }

                  });

          ++bnum;
          if (bnum % 100 == 0) {
              std::cerr << "\r\rprocessed " << bnum << "buckets";
          }

          size_t subBuckets = std::ceil(reads.size() / 256.0);

          for (size_t i = 0; i < subBuckets; ++i) {
              auto curr = reads.begin() + 256 * i;
              size_t remaining = std::distance(curr, reads.end());
              size_t sbsize = (remaining > 256) ? 255 : (remaining - 1);
              auto key = k.get_bits(0, 2*kl);
              outfile.write(reinterpret_cast<char*>(&key), sizeof(key) );
              uint8_t subBucketSize = static_cast<uint8_t>(sbsize);
              outfile.write(reinterpret_cast<char*>(&subBucketSize), sizeof(subBucketSize));

              size_t encReadLength = readLength - kl;
              size_t extraBasesPerRead = encReadLength % 4;
              std::stringstream sstream;
              extraBasesPerRead = 0;
              for (size_t j = 0; j <= sbsize; ++j) {
                  auto& r = reads[i*256 + j];
                  std::string rstr = std::get<0>(r);

                  try {
                      //std::vector<uint8_t> bytes = twoBitEncode(sstream.str());
                      std::vector<uint8_t> bytes = mince::utils::twoBitEncode(rstr.substr(0, encReadLength - extraBasesPerRead));
                      outfile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
                      if (extraBasesPerRead > 0) {
                          sstream << rstr.substr(encReadLength - extraBasesPerRead);
                      }
                  } catch (std::exception& e ) {
                      std::cerr << "EXCEPTION: " << e.what() << rstr << "\n";
                  }
              }

              if (extraBasesPerRead > 0) {
                  std::vector<uint8_t> bonusBytes = mince::utils::twoBitEncode(sstream.str());
                  outfile.write(reinterpret_cast<char*>(&bonusBytes[0]), sizeof(uint8_t) * bonusBytes.size());
              }

              std::vector<uint8_t> offsets;
              offsets.reserve(sbsize);
              for (size_t j = 0; j <= sbsize; ++j) {
                  auto& r = reads[i*256 + j];
                  uint8_t l = static_cast<uint8_t>(std::get<1>(r));
                  offsets.push_back(l);
                  totOutput++;
              }

              std::vector<uint8_t> deltas;
              deltas.reserve(sbsize);
              deltas.push_back(offsets[0]);
              for (size_t j = 1; j <= sbsize; ++j) {
                  int diff = (offsets[j] - offsets[j-1]);
                  if (diff < 0) {
                      std::cerr << "offsets[" << j << "] = " << offsets[j]
                          << "offsets[" << j-1 << "] = " << offsets[j-1] << "\n";
                      std::exit(1);
                  }

                  deltas.push_back(static_cast<uint8_t>(diff));
              }
              outfile.write(reinterpret_cast<char*>(&deltas[0]), deltas.size());
          }
          totBuckets++;
      }
      outfile.close();

      std::cerr << "Total bases: " << total << "\n";
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


