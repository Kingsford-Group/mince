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

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

#include <tbb/concurrent_unordered_map.h>

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

struct CountList {
    uint32_t count;
    std::vector<uint32_t> rlist;
};
using CountMap = tbb::concurrent_unordered_map<uint32_t, uint32_t>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
using MapT = std::map<my_mer, std::vector<std::tuple<std::string, uint8_t>>>;

std::string twoBitDecode(uint8_t* read, size_t readLength) {

   // std::cerr << "here :: readLength = " << readLength << "\n";
    boost::dynamic_bitset<uint8_t> bitmap(2*readLength, 0);
    size_t effectiveReadLength = bitmap.num_blocks();
    //std::cerr << "effectiveReadLength = " << effectiveReadLength << ", bitmap.size() " << bitmap.size()
    //          << ", bitmap.num_blocks() = " << bitmap.num_blocks() << "\n";
    //std::cerr << "read = " << read << "\n";
    boost::from_block_range(read, read + effectiveReadLength, bitmap);
    std::string readOut(readLength, 'X');

    for (size_t i = 0; i < 2*readLength; i += 2) {
        // A
        if (!bitmap[i] and !bitmap[i+1]) {
            readOut[i/2] = 'A';
        }
        // C
        else if (!bitmap[i] and bitmap[i+1]) {
            readOut[i/2] = 'C';
        }
        // G
        else if (bitmap[i] and !bitmap[i+1]) {
            readOut[i/2] = 'G';
        }
        // T
        else if (bitmap[i] and bitmap[i+1]) {
            readOut[i/2] = 'T';
        }
    }
    //std::cerr << "readOut: " << readOut << "\n";
    return readOut;
}

std::vector<uint8_t> twoBitEncode(const std::string& str) {
    boost::dynamic_bitset<uint8_t> bitmap;
    for (auto c : str) {
        switch(c) {
            case 'A':
            case 'N':
                bitmap.push_back(0);
                bitmap.push_back(0);
                break;
            case 'C':
                bitmap.push_back(0);
                bitmap.push_back(1);
                break;
            case 'G':
                bitmap.push_back(1);
                bitmap.push_back(0);
                break;
            case 'T':
                bitmap.push_back(1);
                bitmap.push_back(1);
                break;
        }
    }
    std::vector<uint8_t> bytes;
    boost::to_block_range(bitmap, std::back_inserter(bytes));
    return bytes;
}

void reverseComplement(std::string& s) {
    std::string revs(s.length(), 'X');
    size_t N = s.length();
    size_t o = N;
    for (size_t i = 0; i < N; ++i) {
        --o;
        switch (s[i]) {
            case 'A':
                revs[o] = 'T';
                break;
            case 'C':
                revs[o] = 'G';
                break;
            case 'G':
                revs[o] = 'C';
                break;
            case 'T':
                revs[o] = 'A';
                break;
            default:
                revs[o] = 'N';
                break;
        }
    }
    std::swap(revs, s);
}

void decode(std::string& ifname, std::string& ofname) {
    size_t kl{15};
    std::ifstream ifile;
    uint8_t rl{0};
    ifile.open(ifname, std::ios::in | std::ios::binary);
    ifile.read(reinterpret_cast<char*>(&rl), sizeof(rl));
    //std::cerr << "READ LENGTH = " << +rl << "\n";
    std::ofstream ofile;
    ofile.open(ofname);

    size_t readLength = rl;

    uint64_t bucketMer{0};
    std::string bucketString('X', kl);
    my_mer bucketKmer;
    my_mer::k(kl);

    size_t effectiveReadLength = std::ceil((readLength - kl) / 4.0);
    uint8_t* read = new uint8_t[effectiveReadLength];
    uint8_t offset{0};
    size_t i{0};
    size_t subBucketSize{0};
    while (ifile.good()) {
        //std::cerr << "cha!\n";
        ifile.read(reinterpret_cast<char*>(&bucketMer), sizeof(bucketMer));
        uint8_t bsize{0};
        ifile.read(reinterpret_cast<char*>(&bsize), sizeof(bsize));
        subBucketSize = bsize + 1;
        bucketKmer.set_bits(0, 2*kl, bucketMer);
        bucketString = bucketKmer.to_str();
        //std::cerr << "new bucket --- kmer = " << bucketString << " count = " << +subBucketSize << "\n";

        size_t subBucketCount{0};
        while (subBucketCount++ < subBucketSize) {
            //std::cerr << "next read\n";
            ifile.read(reinterpret_cast<char*>(read), effectiveReadLength);
            ifile.read(reinterpret_cast<char*>(&offset), sizeof(uint8_t));
            size_t newOffset = readLength - kl - offset;
            /*std::cerr << "is: ";
            for (size_t j = 0; j < effectiveReadLength; ++j) {
                std::cerr << +read[j];
            }
            std::cerr << "\n";
            */
            std::string s(twoBitDecode(read, readLength - kl));
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

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename MerT>
void add_sizes(sequence_parser* parser, std::atomic<uint64_t>* total,
               std::atomic<uint64_t>& totReads,
               MapT& buckets,
               CountMap& countMap) {
  uint64_t count = 0;
  unsigned int kl{15};
  unsigned int cmlen{0};

  //unsigned int kl{15};
  //using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
  //my_mer::k(kl);

  //std::map<MerT, std::vector<std::tuple<std::string, size_t>>> buckets;
  while(true) {
      sequence_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
      if(j.is_empty()) break;          // If got nothing, quit

      for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
          if (totReads % 1000000 == 0) {
              std::cerr << "\r\rprocessed " << totReads << " reads";
          }
          ++totReads;
          std::string s(j->data[i].seq);

          count += s.size();        // Add up the size of the sequence

          cmlen = 0;
          MerT mer(kl);
          MerT minmer(kl);
          MerT rmer(kl);
          minmer.polyT();
          mer.polyT();
          rmer.polyT();

          enum Direction { FORWARD, REVERSE };
          size_t readLength{s.length()};
          std::map<uint32_t, std::tuple<uint8_t, Direction>> merOffsetMap;
          size_t offset{0};
          size_t firstOffset{offset};
          while (offset < s.size()) {
            int c = jellyfish::mer_dna::code(s[offset]);
            /*
            if (jellyfish::mer_dna::not_dna(c)) {
                c = static_cast<int>('A');
                s[offset] = 'A';
            }
            */
            mer.shift_left(c);
            rmer.shift_right(jellyfish::mer_dna::complement(c));
            if (jellyfish::mer_dna::not_dna(c)) {
                cmlen = 0;
                ++offset;
                continue;
            }
            if (++cmlen >= kl) {
                cmlen = kl;

                auto key = mer.get_bits(0, 2*kl);
                auto rkey = rmer.get_bits(0, 2*kl);

                auto it = merOffsetMap.find(key);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[key] = std::make_tuple(offset - (kl - 1), Direction::FORWARD);
                }

                it = merOffsetMap.find(rkey);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[rkey] = std::make_tuple(readLength - offset - 1, Direction::REVERSE);
                }

                /*
                if (mer < minmer) {
                    minmer = mer;
                    if (offset < kl - 1) {
                        std::cerr << "WHAT: offset = " << offset << ", but kl = " << kl << "\n";
                    }
                    firstOffset = offset-(kl - 1);
                }
                */
            }
            ++offset;
          }

          /* decide the bucket */
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
              reverseComplement(s);
          }
          /*
          if (cit == countMap.end()) {
              countMap[maxBucketKey] = {1};
          } else {
              countMap[maxBucketKey]++;
          }*/
          countMap[maxBucketKey]++;
          minmer.set_bits(0, 2*kl, maxBucketKey);
          firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);
          if (firstOffset + kl < s.size()) {
              std::string reord(s.substr(firstOffset + kl) + s.substr(0, firstOffset));
              if (s.substr(firstOffset, kl) != minmer.to_str()) {
                std::cerr << "string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
              }

              buckets[minmer].push_back(make_tuple(reord, firstOffset));
              //buckets[minmer].push_back(make_tuple(s, firstOffset));
              /*
              auto pos = reord.find('*');
              std::string recon = reord.substr(pos + 1) + minmer.to_str() + reord.substr(0, pos);
              if (recon != s) {
                  std::cerr << "recon: " << recon << "\n";
                  std::cerr << "orig:  " << s << "\n";
              }
              */
          } else {

              //buckets[minmer].push_back(make_tuple(s, firstOffset));
              buckets[minmer].push_back(make_tuple(s.substr(0, firstOffset), firstOffset));
              if (s.substr(firstOffset, kl) != minmer.to_str()) {
                  std::cerr << "string contains " << s.substr(firstOffset, kl) << ", but mer is " << minmer.to_str() << ", read = " << s << ", offset = " << firstOffset << "\n";
                  std::cerr << "read is in the " << ((maxBucketDirection == Direction::FORWARD) ? "forward" : "reverse") << " direction\n";
              }


          }
      }
  }

  *total += count;
}

double angularDistance(std::string& a, std::string& b) {

    std::map<uint32_t, uint32_t> keyToCountA;
    std::map<uint32_t, uint32_t> keyToCountB;

    size_t readLength{a.length()};
    size_t offset{0};
    size_t cmlen{0};
    size_t firstOffset{offset};
    size_t kl{5};
    my_mer amer, bmer;
    my_mer::k(kl);
    amer.polyT();
    bmer.polyT();

    while (offset < a.size()) {

        int ca = jellyfish::mer_dna::code(a[offset]);
        int cb = jellyfish::mer_dna::code(b[offset]);

        if (jellyfish::mer_dna::not_dna(ca)) {
            ca = static_cast<int>('A');
            a[offset] = 'A';
        }
        if (jellyfish::mer_dna::not_dna(cb)) {
            cb = static_cast<int>('A');
            b[offset] = 'A';
        }
        amer.shift_left(ca);
        bmer.shift_left(cb);
        if (++cmlen >= kl) {
            cmlen = kl;

            auto keya = amer.get_bits(0, 2*kl);
            auto keyb = bmer.get_bits(0, 2*kl);
            keyToCountA[keya] += 1;
            keyToCountB[keyb] += 1;
        }
        ++offset;
    }

    double normA{0.0};
    double normB{0.0};
    double aDotB{0.0};
    for (auto& kv : keyToCountA) {
        normA += kv.second * kv.second;
        auto it = keyToCountB.find(kv.first);
        if (it != keyToCountB.end()) {
            aDotB += kv.second * it->second;
        }
    }
    for (auto& kv : keyToCountB) {
        normB += kv.second * kv.second;
    }

    normA = std::sqrt(normA);
    normB = std::sqrt(normB);
    double sim = std::min(aDotB / (normA * normB), 1.0);
    const double pi = boost::math::constants::pi<double>();
    std::cerr << "a . b = " << aDotB << ", normA = " << normA << ", normB = " << normB << ", sim = " << sim  << "\n";
    return 1.0 - (1.0 - (std::acos(sim)) / pi);
}

double sqEuclideanDistance(std::string& a, std::string& b) {

    std::map<uint32_t, uint32_t> keyToCountA;
    std::map<uint32_t, uint32_t> keyToCountB;

    size_t readLength{a.length()};
    size_t offset{0};
    size_t cmlen{0};
    size_t firstOffset{offset};
    size_t kl{15};
    my_mer amer, bmer;
    my_mer::k(kl);
    amer.polyT();
    bmer.polyT();

    while (offset < a.size()) {

        int ca = jellyfish::mer_dna::code(a[offset]);
        int cb = jellyfish::mer_dna::code(b[offset]);

        if (jellyfish::mer_dna::not_dna(ca)) {
            ca = static_cast<int>('A');
            a[offset] = 'A';
        }
        if (jellyfish::mer_dna::not_dna(cb)) {
            cb = static_cast<int>('A');
            b[offset] = 'A';
        }
        amer.shift_left(ca);
        bmer.shift_left(cb);
        if (++cmlen >= kl) {
            cmlen = kl;

            auto keya = amer.get_bits(0, 2*kl);
            auto keyb = bmer.get_bits(0, 2*kl);
            keyToCountA[keya] += 1;
            keyToCountB[keyb] += 1;
        }
        ++offset;
    }

    double sqdiff{0.0};
    for (auto& kv : keyToCountA) {
        auto it = keyToCountB.find(kv.first);
        if (it != keyToCountB.end()) {
            auto diff = kv.second - it->second;
            sqdiff += diff * diff;
            keyToCountB.erase(it);
        } else {
            sqdiff += kv.second * kv.second;
        }
    }
    for (auto& kv : keyToCountB) {
        auto it = keyToCountA.find(kv.first);
        if (it != keyToCountA.end()) {
            auto diff = kv.second - it->second;
            sqdiff += diff * diff;
            keyToCountA.erase(it);
        } else {
            sqdiff += kv.second * kv.second;
        }

    }

    return sqdiff;
}



/** from: http://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices **/
template< typename order_iterator, typename value_iterator >
void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        value_t temp = v[s];
        for ( index_t d2; d != s; d = d2 ) {
            std::swap( temp, v[d] );
            std::swap( order_begin[d], d2 = (diff_t) -1 );
            -- remaining;
        }
        v[s] = temp;
    }
}


template <class T>
void reorder(std::vector<T>& vA, std::vector<size_t>& vI)
{
    size_t i, j, k;
    T t;
    for(i = 0; i < vA.size(); i++){
        if(i != vI[i]){
            t = vA[i];
            k = i;
            while(i != (j = vI[k])){
                // every move places a value in it's final location
                vA[k] = vA[j];
                vI[k] = k;
                k = j;
            }
            vA[k] = t;
            vI[k] = k;
        }
    }
}

void approxTSPReorder(std::vector<std::tuple<std::string, size_t>>& bucketReads) {
    typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
    typedef boost::adjacency_matrix<boost::undirectedS, boost::no_property, EdgeWeightProperty> Graph;
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertIter;

    // Don't consider regions that are too large.  If a bucket is more than 2500 reads, then
    // break it up into smaller segments.
    size_t maxBucketSize{2500};
    auto subBucketStart = bucketReads.begin();
    auto endOfReads = bucketReads.end();
    while (subBucketStart != endOfReads) {
        size_t numReads = std::min( static_cast<size_t>(std::distance(subBucketStart, endOfReads)), maxBucketSize);
        size_t offset = std::distance(bucketReads.begin(), subBucketStart);
        Graph g(numReads);

        #pragma omp parallel for
        for (size_t i = 0; i < numReads; ++i) {
            for (size_t j = i+1; j < numReads; ++j) {
                double d = sqEuclideanDistance(std::get<0>(bucketReads[i + offset]), std::get<0>(bucketReads[j + offset]));
                boost::add_edge(i, j, d, g);
            }
        }

        std::vector<boost::graph_traits<Graph>::vertex_descriptor> t;
        boost::metric_tsp_approx_tour_from_vertex(g, 0,
                get(boost::edge_weight, g), std::back_inserter(t));

        reorder_destructive(t.begin(), t.end() - 1, subBucketStart);
        subBucketStart += numReads;
    }
}


int main(int argc, char *argv[]) {
  const int nb_threads      = 20;
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

  if (argv[1][0] == 'd') {
      std::string input(argv[2]);
      std::string output(argv[3]);
      decode(input, output);
      std::exit(0);
  }


  std::string s("AGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAAGGGGGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATT");
  std::string s2("CGGAAAGTGCCACAGAAAATAACCGCCTAAGCGCAACAGCGCCGGTAAGGGCCGGGTTCTGTCGAGGACAGCCATTCCTCTGGGACGCACATCACTNNACGC");
  std::string s3("CTCCCGTCAATTCATTTGAGTTTTAACCTTGCGGCCGTACTCCCCAGGGGGCCGGGTTCTGTCGAGGACAGCCATTCCTCTGGGACGCACATCACTNNACGC");
  auto ad1 = sqEuclideanDistance(s, s);
  auto ad2 = sqEuclideanDistance(s, s2);
  auto ad3 = sqEuclideanDistance(s, s3);
  auto ad4 = sqEuclideanDistance(s2, s3);

  std::cerr << "angular distance self = " << ad1 << "\n";
  std::cerr << "angular distance diff = " << ad2 << "\n";
  std::cerr << "angular distance diff = " << ad3 << "\n";
  std::cerr << "angular distance diff = " << ad4 << "\n";

  stream_manager  streams(argv + 1, argv + argc, concurrent_file);
  sequence_parser parser(4 * nb_threads, max_read_group, concurrent_file, streams);

  std::string outfname(argv[2]);
  std::vector<std::thread> threads;
  std::atomic<uint64_t> total(0);
  std::atomic<uint64_t> totReads{0};
  std::atomic<uint64_t> distinct{0};
  using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

  std::ofstream outfile;
  outfile.open(outfname, std::ios::out | std::ios::binary);
  unsigned int kl{15};
  my_mer::k(15);

  std::cerr << "\n\n\n" << +std::numeric_limits<uint8_t>::max() << "\n";
  //std::vector<std::map<my_mer, std::vector<std::tuple<std::string, size_t>>>> maps(nb_threads);
  CountMap countMap;
  std::vector<MapT> maps(nb_threads);
  for(int i = 0; i < nb_threads; ++i) {
      threads.push_back(std::thread(add_sizes<my_mer>, &parser, &total, std::ref(totReads), std::ref(maps[i]), std::ref(countMap)));
  }

  for(int i = 0; i < nb_threads; ++i)
    threads[i].join();

  std::set<my_mer> keys;
  for (auto& m : maps) {
      for (auto& kv : m ) {
        keys.insert(kv.first);
      }
  }

  uint8_t readLength = std::get<0>(maps[0].begin()->second.front()).length() + kl;
  std::cerr << "\n\nread length = " << +readLength << "\n";
  //outfile << +readLength << "\n";

  outfile.write(reinterpret_cast<char*>(&readLength), sizeof(readLength));

  size_t totOutput{0};
  size_t totBuckets{0};
  size_t bnum{0};
  std::cerr << "Num buckets = " << keys.size() << "\n";
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


      if (reads.size() > 2) {
        approxTSPReorder(reads);
      }
      /*
         std::sort(reads.begin(), reads.end(),
         [&](const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
         int cmp = std::get<1>(a) - std::get<1>(b);
         if (cmp == 0) {
         return std::get<0>(a) < std::get<0>(b);
         } else {
         return cmp < 0;
         }

         });
         */

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
        for (size_t j = 0; j <= sbsize; ++j) {
            auto& r = reads[i*256 + j];
            std::string rstr = std::get<0>(r);
            uint8_t l = static_cast<uint8_t>(std::get<1>(r));

            try {
                //outfile << rstr << l << "\n";

                std::vector<uint8_t> bytes = twoBitEncode(rstr);
                outfile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());

            } catch (std::exception& e ) {
                std::cerr << rstr << "\n";
            }

            outfile.write(reinterpret_cast<char*>(&l), sizeof(l));

            totOutput++;
        }
    }
    totBuckets++;
  }
  outfile.close();


  std::cerr << "Total bases: " << total << "\n";
  std::cerr << "Wrote: " << totOutput << " reads, in " << totBuckets << " buckets\n";
  return 0;
}
