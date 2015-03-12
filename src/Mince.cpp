#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <functional>
#include <cmath>
#include <mutex>
#include <unordered_set>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>

#include "pstream.h"

#include "bitfile.h"
#include "MinceConfig.hpp"
#include "MinceUtils.hpp"
#include "FindPartition.hpp"
#include "Decoder.hpp"
#include "BucketModel.hpp"
#include "PairSequenceParser.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"
#include "MinceOpts.hpp"

using pair_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<char**>;
using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;

enum Direction { FORWARD, REVERSE };

struct CountList {
    uint32_t count;
    std::vector<uint32_t> rlist;
};

class BucketedString {
    public:
        BucketedString(std::string& strIn, uint8_t offsetIn, bool rcIn, std::vector<uint8_t>&& nlocsIn = std::vector<uint8_t>()) :
            str(strIn), offset(offsetIn), rc(rcIn ? 1 : 0), nlocs(nlocsIn) {}

        BucketedString(const BucketedString& other) {
            str = other.str;
            offset = other.offset;
            rc = other.rc;
	    nlocs = other.nlocs;
        }


        BucketedString& operator=(BucketedString& other) {
            str = other.str;
            offset = other.offset;
            rc = other.rc;
	    nlocs = other.nlocs;
            return *this;
        }

        BucketedString(BucketedString&& other) {
            std::swap(str, other.str);
            offset = other.offset;
            rc = other.rc;
	    std::swap(nlocs, other.nlocs);
        }

        BucketedString& operator=(BucketedString&& other) {
            std::swap(str, other.str);
            offset = other.offset;
            rc = other.rc;
	    std::swap(nlocs, other.nlocs);
            return *this;
        }


    std::string str;
    uint8_t offset;
    uint8_t rc : 1;
    std::vector<uint8_t> nlocs;
};

using CountMap = tbb::concurrent_unordered_map<uint32_t, BucketModel>;
//using CountMap = tbb::concurrent_unordered_map<uint32_t, uint32_t>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
using MapT = std::map<my_mer, std::vector<BucketedString>>;

std::string getSequence(std::pair<header_sequence_qual, header_sequence_qual>& rpair, const ReadLibrary& rl) {

	switch(rl.format().orientation) {
		case ReadOrientation::TOWARD:
			mince::utils::reverseComplement(rpair.second.seq);
			//std::reverse(rp.second.qual.begin(), rp.second.qual.end());
			return rpair.first.seq + rpair.second.seq;
		case ReadOrientation::AWAY:
			mince::utils::reverseComplement(rpair.second.seq);
			//std::reverse(rp.second.qual.begin(), rp.second.qual.end());
			return rpair.first.seq + rpair.second.seq;
		case ReadOrientation::SAME:
			return rpair.first.seq + rpair.second.seq;
		default:
			return rpair.first.seq + rpair.second.seq;
	}

}

std::string getSequence(jellyfish::header_sequence_qual& hsq, const ReadLibrary& rl) {
	return hsq.seq;
}

template <typename MerT,typename ParserT>
void bucketReads(ParserT* parser, ReadLibrary& rl, std::atomic<uint64_t>* total,
        std::atomic<uint64_t>& totReads,
        MapT& buckets,
        CountMap& countMap,
        uint32_t bucketStringLength,
        bool noRC,
        std::mutex& iomutex) {

    uint64_t count = 0;
    unsigned int kl{bucketStringLength};
    unsigned int cmlen{0};

    while(true) {
        typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;          // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
            if (totReads++ % 1000000 == 0) {
                iomutex.lock();
                std::cerr << "\r\rprocessed " << totReads << " reads";
                iomutex.unlock();
            }
            //std::string s(j->data[i].seq);
	        std::string s(getSequence(j->data[i], rl));

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
            //std::unordered_set<uint16_t> sforward = BucketModel::readHash(s, 8, false);
            //std::unordered_set<uint16_t> sreverse = BucketModel::readHash(s, 8, true);

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
            double maxBucketValue{0};
            Direction maxBucketDirection{std::get<1>(merOffsetMap.begin()->second)};
            //auto featVec = mince::utils::trimerVector(s, false);
            //auto featVecRC = mince::utils::trimerVector(s, true);

            // Iterate over all k-mers and reverse complement k-mers
            for (auto& kv : merOffsetMap) {
                auto cit = countMap.find(kv.first);
                // If the bucket for this key is the heaviest so far
                // then it becomes the new bucket key.
                if (cit != countMap.end()) {
                    auto dir = std::get<1>(kv.second);
                    bool rc = dir == Direction::REVERSE;
                    //auto& fvec = (dir == Direction::REVERSE) ? featVecRC : featVec;
                    /*double score;
                    if (rc) {
                        score = cit->second.scoreOfRead(sreverse, 8);
                    } else {
                        score = cit->second.scoreOfRead(sforward, 8);
                    }*/
                    double score = cit->second.scoreOfRead(s, 8, rc);//fvec);
                    if (score > maxBucketValue) {
                        maxBucketValue = score;
                        maxBucketDirection = std::get<1>(kv.second);
                        maxBucketKey = kv.first;
                    }
                }
            }

            // If the heaviest bucket contains the RC mer,
            // then RC the read
            if (maxBucketDirection == Direction::REVERSE) {
                mince::utils::reverseComplement(s);
            }
            // Increment the count of the assigned bucket
            //countMap[maxBucketKey]++;
            countMap[maxBucketKey].addRead(s, 8);//(maxBucketDirection == Direction::REVERSE) ? featVecRC : featVec);

            // Encode the read accounting for the shared string
            minmer.set_bits(0, 2*kl, maxBucketKey);
            // The first appearance of the key in the read
            firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);

            std::vector<uint8_t> nlocs;
            //std::replace(s.begin(), s.end(), 'N', 'A');
            for (uint32_t i = 0; i < s.size(); ++i) {
                if (s[i] == 'N') {
                    nlocs.emplace_back(static_cast<uint8_t>(i));
                    s[i] = 'A';
                }
            }

            // orig: x . key . y
            // new : y . x
            std::string reord = mince::utils::permute(s, firstOffset, kl);

#ifdef DEBUG
            std::string kstr = minmer.to_str();
            if ( s != mince::utils::unpermute(reord, kstr, firstOffset)) {
                std::cerr << "ERROR!!! mince::utils::unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
            }
            if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "A string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
            }
#endif
            // Put the reordered string, offset tuple in the bucket
            // buckets[minmer].push_back(make_tuple(reord, firstOffset));
            buckets[minmer].push_back({reord, static_cast<uint8_t>(firstOffset), (maxBucketDirection == Direction::REVERSE), std::move(nlocs)});
        } // for all of the reads in this job
    } // for all jobs of this thread
}

void reassignOnsies(
        std::vector<std::tuple<my_mer, BucketedString>>& onsies,
        CountMap& countMap,
        MapT& buckets,
        uint32_t bucketStringLength,
        bool noRC
        ) {

    for (auto& kr : onsies) {
        auto& bs = std::get<1>(kr);
        std::string s = bs.str;
        // pre-compute the hash of the minimers of s and rc(s)
        //std::unordered_set<uint16_t> sforward = BucketModel::readHash(s, 8, false);
        //std::unordered_set<uint16_t> sreverse = BucketModel::readHash(s, 8, true);

        std::vector<uint8_t>& nlocs = bs.nlocs;
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
#ifdef DEBUG
                    if (s.substr(offset - kl, kl) != mer.to_str()) {
                        std::cerr << "putting " << mer.to_str() << " into the map, but I really have " << s.substr(offset - kl, kl) << "\n";
                    }
#endif // DEBUG
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
        double maxBucketValue{0};
        Direction maxBucketDirection{std::get<1>(merOffsetMap.begin()->second)};
        //auto featVec = mince::utils::trimerVector(s, false);
        //auto featVecRC = mince::utils::trimerVector(s, true);

        std::atomic<uint32_t> count;
        for (auto& kv : merOffsetMap) {
            auto cit = countMap.find(kv.first);
            if (cit != countMap.end()) {

                auto dir = std::get<1>(kv.second);
                bool rc = dir == Direction::REVERSE;
                //auto& fvec = (dir == Direction::REVERSE) ? featVecRC : featVec;
                /*
                double score;
                if (rc) {
                    score = cit->second.scoreOfRead(sreverse, 8);
                } else {
                    score = cit->second.scoreOfRead(sforward, 8);
                }*/
                double score = cit->second.scoreOfRead(s, 8, rc);//fvec);
                if (score > maxBucketValue) {
                    maxBucketValue = score;
                    maxBucketDirection = std::get<1>(kv.second);
                    maxBucketKey = kv.first;
                }
            }
        }

        if (maxBucketDirection == Direction::REVERSE) {
            mince::utils::reverseComplement(s);
	    size_t N = s.size() - 1;
	    for (size_t i = 0; i < nlocs.size(); ++i) {
		nlocs[i] = N - nlocs[i];
	    }
        }

        //countMap[maxBucketKey]++;
        countMap[maxBucketKey].addRead(s, 8);//(maxBucketDirection == Direction::REVERSE) ? featVecRC : featVec);

        minmer.set_bits(0, 2*kl, maxBucketKey);
        firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);

        std::string reord = mince::utils::permute(s, firstOffset, kl);

#ifdef DEBUG
        std::string kstr = minmer.to_str();
        if ( s != mince::utils::unpermute(reord, kstr, firstOffset)) {
            std::cerr << "ERROR!!! mince::utils::unpermute( " << reord << ", " << kstr << ", " << offset << ") != " << s << "\n";
        }
#endif
        //buckets[minmer].push_back(std::make_tuple(reord, firstOffset));
	bool canonicalFlip = (bs.rc xor (maxBucketDirection == Direction::REVERSE));
        buckets[minmer].push_back({reord, static_cast<uint8_t>(firstOffset), canonicalFlip, std::move(nlocs)});
    }

}

int main(int argc, char *argv[]) {
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  using std::string;

  MinceOpts minceOpts;
  bool doDecode{false};
  bool doEncode{false};
  bool noRC{false};

  uint32_t bucketStringLength;

  po::options_description generic("mince options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      //("input,i", po::value<string>(), "input file; FASTA/Q if encoding, MINCE if decoding")
      ("threads,p", po::value<uint32_t>(&minceOpts.numThreads)->default_value(20),
                 "number of concurrent threads to use for compression / decompression; "
                 "the minimum value is 4.")
      ("libtype,l", po::value<string>(), "library format string [only for encoding]")
      ("mates1,1", po::value<std::vector<string>>(), "mate file 1 [only for encoding]")
      ("mates2,2", po::value<std::vector<string>>(), "mate file 2 [only for encoding]")
      ("unmated_reads,r", po::value<std::vector<string>>(), "unmated reads [only for encoding]")
      ("input,i", po::value<string>(), "input base file (same as -o option passed in when encoding) [only for decoding]")
      ("output,o", po::value<string>(), "output base file [for both encoding / decoding]")
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
= Mince =
---------
)";
          std::cerr << hstring << "\n";
          std::cerr << generic << "\n";
          std::cerr << "Mince v" << mince::version << "\n";
          std::exit(0);
      }

      if (vm.count("version")) {
          std::cout << "Mince version : " << mince::version << "\n";
          std::exit(0);
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

      std::string outfname = vm["output"].as<string>();

      bfs::path logFilePath = outfname + ".log";
      size_t maxLogQueueSize = 2097152;
      spdlog::set_async_mode(maxLogQueueSize);

      auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(
              logFilePath.string(), true);
      auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
      auto consoleLog = spdlog::create("consoleLog", {consoleSink});
      auto fileLog = spdlog::create("fileLog", {fileSink});
      auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

      minceOpts.fileLog = fileLog;
      minceOpts.jointLog = jointLog;

      if (minceOpts.numThreads < 4) {
          jointLog->info("You requested {} threads, but the minimum is 4. "
                         "Setting # threads to 4.", minceOpts.numThreads);
          minceOpts.numThreads = 4;
      }

      if (doDecode) {
          std::string input = vm["input"].as<string>();
          std::string output = vm["output"].as<string>();
          Decoder d;
          d.decode(minceOpts, input, output);
          std::exit(0);
      }


      /*
      // Currently only single-end or pre-merged paired-end reads
      // No reason we can't do the merging on the fly and accept paired-end reads directly
      char* input[] = { const_cast<char*>(vm["input"].as<string>().c_str()) };

      stream_manager  streams(input, input + 1, concurrent_file);
      sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);
      */


      std::vector<std::thread> threads;
      std::atomic<uint64_t> total(0);
      std::atomic<uint64_t> totReads{0};
      std::atomic<uint64_t> distinct{0};
      using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

      unsigned int kl{bucketStringLength};
      my_mer::k(kl);

      // It's a small likelihood, but prevent us from messing up i/o
      // by writing from multiple threads
      std::mutex iomutex;
      CountMap countMap;

      std::vector<MapT> maps(minceOpts.numThreads);
      jointLog->info("Starting {} processing threads", minceOpts.numThreads);

      auto readLibraries = mince::utils::extractReadLibraries(orderedOptions);
      assert(readLibraries.size() == 1);
      std::vector<const char*> inputs;
      std::unique_ptr<pair_parser> pparser(nullptr);
      std::unique_ptr<stream_manager> streams(nullptr);
	  std::unique_ptr<sequence_parser> sparser(nullptr);

      for (auto& rl : readLibraries) {
	      if (rl.format().type == ReadType::PAIRED_END) {
		      auto& m1 = rl.mates1();
		      auto& m2 = rl.mates2();
		      for (size_t i = 0; i < m1.size(); ++i) {
			      inputs.push_back(m1[i].c_str());
			      std::cerr << "merging " << inputs.back() << " with ";
			      inputs.push_back(m2[i].c_str());
			      std::cerr << inputs.back() << "\n";
		      }

		      char** start = const_cast<char**>(&(*inputs.begin()));
		      char** stop = const_cast<char**>(&(*inputs.end()));
		      pparser.reset(new pair_parser(4 * minceOpts.numThreads, max_read_group, concurrent_file,
				      start, stop));
		      // Spawn off the read parsing threads
		      for(int i = 0; i < minceOpts.numThreads; ++i) {
                  threads.emplace_back(bucketReads<my_mer, pair_parser>,
                          pparser.get(), std::ref(rl), &total,
                          std::ref(totReads), std::ref(maps[i]),
                          std::ref(countMap), bucketStringLength,
                          noRC, std::ref(iomutex));
		      }

	      } else if (rl.format().type == ReadType::SINGLE_END) {
		      auto& m1 = rl.unmated();
		      for (size_t i = 0; i < m1.size(); ++i) {
			      inputs.push_back(m1[i].c_str());
		      }

		      char** start = const_cast<char**>(&(*inputs.begin()));
		      char** stop = const_cast<char**>(&(*inputs.end()));
              streams.reset(new stream_manager(start, stop, concurrent_file));
		      sparser.reset(new sequence_parser(4 * minceOpts.numThreads, max_read_group, concurrent_file, *(streams.get())));

		      // Spawn off the read parsing threads
		      for(int i = 0; i < minceOpts.numThreads; ++i) {
                  threads.emplace_back(bucketReads<my_mer, sequence_parser>,
                          sparser.get(), std::ref(rl), &total,
                          std::ref(totReads), std::ref(maps[i]),
                          std::ref(countMap), bucketStringLength,
                          noRC, std::ref(iomutex));
		      }

	      }


      }

      // Join all of the read-parsing threads
      for(int i = 0; i < minceOpts.numThreads; ++i) {
          threads[i].join();
      }
      jointLog->info("finished parsing files");

      // The union of keys from all maps
      std::set<my_mer> keys;
      for (auto& m : maps) {
          for (auto& kv : m ) {
              keys.insert(kv.first);
          }
      }

      jointLog->info("Marking and collecting onsies");

      // The length of the reads in this file
      uint8_t readLength = maps[0].begin()->second.front().str.length() + kl;
      jointLog->info("Read length (assumed = for all reads) is {}",
            +readLength);

      size_t totOutput{0};
      size_t totBuckets{0};
      size_t bnum{0};
      size_t totNumReads{0};

      // The total # of buckets
      std::cerr << "Num buckets = " << keys.size() << "\n";

      std::vector<std::tuple<my_mer, BucketedString>> onsies;
      // For all buckets
      for (auto& k : keys) {
          std::vector<std::tuple<my_mer, BucketedString>> combined;
          bool isOnsie{true};

          // For all maps
          for (auto& m : maps) {
              // If this map contains this bucket, and this bucket is of size one
              // in this map, then the contained read is a potential onsie
              auto kit = m.find(k);
              if (kit != m.end()) {
                  if (kit->second.size() == 1) {
                      combined.push_back(std::forward_as_tuple(k, kit->second.front()));
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
              std::string s(std::get<1>(combined[0]).str);
              // Recover the original string
              size_t offset = std::get<1>(combined[0]).offset;
	          bool rc = std::get<1>(combined[0]).rc;
    	      std::vector<uint8_t>& nlocs = std::get<1>(combined[0]).nlocs;
              std::string recon = mince::utils::unpermute(s, bucketString, offset);

              // Push back the original key and the reconstructed string
              onsies.push_back(std::forward_as_tuple(k, BucketedString(recon, offset, rc, std::move(nlocs))));
          }
      }

      jointLog->info("erasing onsies from individual maps");

      // For all onsies
      for (auto& kr : onsies) {
          // erase this guy from the map and the set of keys
          for (auto& m : maps) {
              auto key = std::get<0>(kr);
              auto kit = m.find(key);
              if (kit != m.end()) {
                  if (kit->second.size() > 1) {
                      jointLog->critical("Erasing non-singleton bucket ({})!",
                              kit->second.size());
                      throw std::logic_error("Erasing non-singleton bucket!");
                  }
                  m.erase(kit);
                  keys.erase(key);
                  countMap[key.get_bits(0, 2*kl)].subCount();
              }
          }
      }

      jointLog->info("|keys| = {}", keys.size());
      jointLog->info("total number of reads = {}", totNumReads);
      jointLog->info("total number of onsies = {}", onsies.size());

      jointLog->info("re-assigning onsies");

      // Re-assign the onsies to see if they can be placed in any non-empty
      // buckets that didn't exist when they were first considered.
      MapT collatedMap;
      reassignOnsies(onsies, countMap, collatedMap, bucketStringLength, noRC);

      size_t newOnsies{0};
      size_t newTotal{0};
      std::set<my_mer> onsieKeys;
      std::vector<std::tuple<my_mer, BucketedString>> onsieVec;
      for (auto& kv : collatedMap) {
          keys.insert(kv.first);
          if (kv.second.size() == 1) {
              ++newOnsies;
              onsieKeys.insert(kv.first);
              std::string bucketString = kv.first.to_str();
              BucketedString& so = kv.second.front();
              std::string recon = mince::utils::unpermute(so.str, bucketString, so.offset);
              onsieVec.push_back(std::forward_as_tuple(kv.first, std::move(BucketedString(recon, so.offset, so.rc, std::move(so.nlocs)))));
          }
          newTotal += kv.second.size();
      }

      jointLog->info("done re-assigning onsies");

      jointLog->info("|keys| = {}", keys.size());
      jointLog->info("now there are {} onsies", newOnsies);
      jointLog->info("re-bucked {} reads", newTotal);

      jointLog->info("sorting onsies");
      tbb::parallel_sort(onsieVec.begin(), onsieVec.end(),
              [](const std::tuple<my_mer,
                  BucketedString>& o1,
                  const std::tuple<my_mer,
                  BucketedString>& o2) -> bool {
              auto& s1 = const_cast<std::string&>(std::get<1>(o1).str);
              auto& s2 = const_cast<std::string&>(std::get<1>(o2).str);
              for (auto i1 = s1.rbegin(), i2 = s2.rbegin(); i1 != s1.rend() and i2 != s2.rend(); ++i1, ++i2) {
              if (*i1 != *i2) {
              return *i1 < *i2;
              }
              }
              return s1.size() < s2.size();
              });

      jointLog->info("done sorting onsies");

      // The sequence information is output to four streams
      // The sequence & control stream --- ending in .seqs
      /*
      std::ofstream outfileSeqs;
      outfileSeqs.open(outfname+".seqs", std::ios::out | std::ios::binary);
      */
      fmt::MemoryWriter w;
      w.write("plzip -o {} -f -n {} -", outfname+".seqs", minceOpts.numThreads - 3);
      jointLog->info("writing seq buffer with command {}", w.str());
      redi::opstream outfileSeqs(w.str());
      w.clear();

      // the offset stream --- ending in .offs
      /*
      std::ofstream outfileOffsets;
      outfileOffsets.open(outfname+".offs", std::ios::out | std::ios::binary);
      */
      w.write("plzip -o {} -f -n 1 -", outfname+".offs");
      jointLog->info("writing offset buffer with command {}", w.str());
      redi::opstream outfileOffsets(w.str());
      w.clear();

      // the flip stream --- just a bunch of bits denoting if the read was RC'd or not
      const char* flipFileName = (outfname+".flips").c_str();
      bit_file_c outfileFlips(flipFileName, BF_WRITE);

      // the Ns stream --- encodes the locations where an N has been flipped in the
      // encoded sequence.
      /*
      std::ofstream outfileNLocs;
      outfileNLocs.open(outfname+".nlocs", std::ios::out | std::ios::binary);
      */
      w.write("plzip -o {} -f -n 1 -", outfname+".nlocs");
      jointLog->info("writing nloc buffer with command {}", w.str());
      redi::opstream outfileNLocs(w.str());
      w.clear();



      // The read library type
      uint8_t rltype = readLibraries[0].format().formatID();
      outfileSeqs.write(reinterpret_cast<char*>(&rltype), sizeof(rltype));
      // The length of the string on which we'll bucket
      uint8_t bslen = bucketStringLength;
      outfileSeqs.write(reinterpret_cast<char*>(&bslen), sizeof(bslen));
      // The length of the reads
      outfileSeqs.write(reinterpret_cast<char*>(&readLength), sizeof(readLength));

      uint32_t readCtr{0};
      // First, write out all remaining onsies
      outfileSeqs.write(reinterpret_cast<char*>(&newOnsies), sizeof(newOnsies));
      for (auto& bs : onsieVec) {
          auto& k = std::get<0>(bs);
	  auto& br = std::get<1>(bs);
          auto& recon = br.str;
          auto bytes = mince::utils::twoBitEncode(recon);
	  if (br.nlocs.size() > 0) {
	  	outfileNLocs.write(reinterpret_cast<char*>(&readCtr), sizeof(readCtr));
		uint8_t nn = br.nlocs.size();
		outfileNLocs.write(reinterpret_cast<char*>(&nn), sizeof(nn));
		outfileNLocs.write(reinterpret_cast<char*>(&br.nlocs[0]), nn * sizeof(br.nlocs[0]));
	  }
	  ++readCtr;
          // BINARY
          outfileSeqs.write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
          // ASCII
          //outfileSeqs.write(&recon[0], recon.length());

          outfileFlips.PutBit(br.rc);
          totOutput++;
          collatedMap.erase(k);
      }

      // Push back the remaining onsies map
      maps.push_back(collatedMap);

      jointLog->info("wrote out onsies");

      // Go through all of the buckets
      for (auto& k : keys ) {
          // Collect all of the reads for this bucket
          std::vector<BucketedString> reads;
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
                    [&](const BucketedString& a, const BucketedString& b) -> bool {
                        if (a.offset == b.offset) {
                            return a.str < b.str;
                        } else {
                            return a.offset < b.offset;
                        }
                  });

          // Report bucket progress
          ++bnum;
          if (bnum % 1000 == 0) {
              std::cerr << "\r\rprocessed " << bnum << " buckets";
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
              fileLog->info("writing bucket of size {}", sbsize);

              if (nextReadIdx > lastReadIdx) {
                  fmt::MemoryWriter w;
                  w.write("sub-bucket {} had non-positive range!", i);
                  jointLog->critical(w.str());
                  throw std::logic_error(w.str());
              }

              std::string& firstRead = reads[nextReadIdx].str;
              std::string& lastRead = reads[lastReadIdx].str;

              size_t offset1 = reads[nextReadIdx].offset;
              size_t offset2 = reads[lastReadIdx].offset;

              if (nextReadIdx == lastReadIdx) { lcpLen = 0; }
              std::string keyStr = k.to_str() + reads[nextReadIdx].str.substr(0, lcpLen);
              auto keyStrLen = keyStr.length();

              if (offset1 > firstRead.length()) {
                  fmt::MemoryWriter w;
                  w.write("offset1 = {} > read length = {}",
                          offset1, firstRead.length());
                  jointLog->critical(w.str());
                  throw std::logic_error(w.str());
              }

              if (offset2 > lastRead.length()) {
                  fmt::MemoryWriter w;
                  w.write("offset2 = {} > read length = {}",
                          offset2, lastRead.length());
                  jointLog->critical(w.str());
                  throw std::logic_error(w.str());
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
                  std::string rstr = r.str;
                  uint8_t offset = r.offset;

                  try {
                      std::string rsub = rstr;//rstr.substr(lcpLen, encReadLength - extraBasesPerRead);
                      if (offset > rsub.size()) {
                          jointLog->info("bucket [{}, {}]", nextReadIdx, lastReadIdx);
                          jointLog->info("lcpLen = {}", lcpLen);
                          jointLog->info("encountered an offset {} that is greater than the string length",
                                  +offset, rsub.length());
                      }

                      // BINARY
                      std::vector<uint8_t> bytes = mince::utils::twoBitEncode(rsub);
                      outfileSeqs.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
                      // ASCII
                      // outfileSeqs.write(&rsub[0], rsub.length());

        		      outfileFlips.PutBit(r.rc);
                      if (r.nlocs.size() > 0) {
                          outfileNLocs.write(reinterpret_cast<char*>(&readCtr), sizeof(readCtr));
                          uint8_t nn = r.nlocs.size();
                          outfileNLocs.write(reinterpret_cast<char*>(&nn), sizeof(nn));
                          outfileNLocs.write(reinterpret_cast<char*>(&r.nlocs[0]), nn * sizeof(r.nlocs[0]));
                      }
                      ++readCtr;
                      if (extraBasesPerRead > 0) {
                          sstream << rstr.substr(encReadLength - extraBasesPerRead);
                      }

                  } catch (std::exception& e ) {
                      jointLog->critical("EXCEPTION: {}", e.what());
                  }
              }

              // Write out the offsets (delta-encoded) of the bucket strings within the reads
              std::vector<uint8_t> offsets;
              offsets.reserve(sbsize);
              for (size_t j = nextReadIdx; j <= lastReadIdx; ++j) {
                  auto& r = reads[j];
                  uint8_t l = r.offset;

                  if (l > r.str.length()) {
                      jointLog->info("offset is {}, which is greater than text length {}",
                              l, r.str.length());
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
                      fmt::MemoryWriter w;
                      w.write("offsets[{}] = {}"
                              "offsets[{}] = {}",
                              j, +offsets[j],
                              j-1, +offsets[j-1]);
                      jointLog->critical(w.str());
                      throw std::logic_error(w.str());
                  }
                  deltas.push_back(static_cast<uint8_t>(diff));
              }
              outfileOffsets.write(reinterpret_cast<char*>(&deltas[0]), deltas.size());
          }

          // Increment the total bucket count
          totBuckets++;
      }

      outfileSeqs << redi::peof;
      outfileOffsets << redi::peof;
      outfileNLocs << redi::peof;

      jointLog->info("Done with calls to write; waiting for files to close!");

      outfileSeqs.close();
      outfileOffsets.close();
      outfileNLocs.close();
      outfileFlips.Close();

      // The "flips" file is small enough that it's first written
      // uncompressed.  We call plzip on it here.
      {
          fmt::MemoryWriter zipCmd;
          zipCmd.write("plzip -f -n 1 {}", outfname+".flips");
          redi::opstream zipFlipFile(zipCmd.str());
          if (!zipFlipFile.good()) {
              jointLog->critical("Error compressing flip file!");
          }
      }
      // done compressing flip file

      jointLog->info("Wrote {} reads in {} buckets", totOutput, totBuckets);
      jointLog->info("Done writing all buckets");
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


