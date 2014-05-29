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
void bucketReads(sequence_parser* parser, std::atomic<uint64_t>* total,
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
          if (totReads++ % 1000000 == 0) {
              std::cerr << "\r\rprocessed " << totReads << " reads";
          }
          std::string s(j->data[i].seq);

          count += s.size();        // Add up the size of the sequence

          cmlen = 0;
          MerT mer(kl);
          MerT minmer(kl);
          MerT rmer(kl);
          minmer.polyT();
          mer.polyT();
          rmer.polyT();

         size_t readLength{s.length()};
          std::map<uint32_t, std::tuple<uint8_t, Direction>> merOffsetMap;
          size_t offset{0};
          size_t firstOffset{offset};
          while (offset < s.size()) {
            int c = jellyfish::mer_dna::code(s[offset]);

            if (jellyfish::mer_dna::not_dna(c)) {
                c = jellyfish::mer_dna::code('A');
                s[offset] = 'A';
            }

            mer.shift_left(c);
            rmer.shift_right(jellyfish::mer_dna::complement(c));

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
            }
            ++offset;
          }

          // Decide the heaviest bucket
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

          // If the heaviest bucket contains the RC mer,
          // then RC the read
          if (maxBucketDirection == Direction::REVERSE) {
              reverseComplement(s);
          }
          // Increment the count of the assigned bucket
          countMap[maxBucketKey]++;

          // Encode the read accounting for the shared string
          minmer.set_bits(0, 2*kl, maxBucketKey);
          // The first appearance of the key in the read
          firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);
          // If the first appearance isn't at the end if the read
          if (firstOffset + kl < s.size()) {
              // orig: x . key . y
              // new : y . x
              std::string reord(s.substr(firstOffset + kl) + s.substr(0, firstOffset));
              // Put the reordered string, offset tuple in the bucket
              buckets[minmer].push_back(make_tuple(reord, firstOffset));
          } else { // if the first appearance is at the end of the string, then the 'y' substring is empty
              buckets[minmer].push_back(make_tuple(s.substr(0, firstOffset), firstOffset));
          }
      } // for all of the reads in this job
  } // for all jobs of this thread
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

    return std::sqrt(sqdiff);
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

    std::vector<std::tuple<std::string, size_t>> newBucketReads;
    newBucketReads.reserve(bucketReads.size());
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
            /*
            StripedSmithWaterman::Aligner aligner;
            StripedSmithWaterman::Filter filter;
            StripedSmithWaterman::Alignment alignment;
            */
            for (size_t j = i+1; j < numReads; ++j) {
                double d = sqEuclideanDistance(std::get<0>(bucketReads[i + offset]), std::get<0>(bucketReads[j + offset]));
                /*
                std::string& query = std::get<0>(bucketReads[i + offset]);
                std::string& reference = std::get<0>(bucketReads[j + offset]);
                aligner.Align(query.c_str(), reference.c_str(), reference.size(), filter, &alignment);
                auto score = alignment.sw_score;
                double d = 2.0 * reference.size() - score;
                */
                boost::add_edge(i, j, d, g);
            }
        }

        std::vector<boost::graph_traits<Graph>::vertex_descriptor> t;
        boost::metric_tsp_approx_tour_from_vertex(g, 0,
                get(boost::edge_weight, g), std::back_inserter(t));

        //reorder_destructive(t.begin(), t.end() - 1, subBucketStart);
        for (auto it = t.begin(); it != t.end() - 1; ++it) {
            newBucketReads.push_back(bucketReads[*it + offset]);
        }
        subBucketStart += numReads;
    }
    std::swap(bucketReads, newBucketReads);
}

void subBucket(std::vector<std::tuple<std::string, size_t>>& bucketReads) {
    size_t kl = 15;
    std::map<uint32_t, std::vector<uint32_t>> subMap;
    for (size_t i=0; i < bucketReads.size(); ++i) {
        size_t cmlen = 0;
        my_mer mer(kl);
        my_mer::k(kl);
        mer.polyT();

        auto& s = std::get<0>(bucketReads[i]);
        size_t readLength{s.length()};
        std::set<uint32_t> mers;
        size_t offset{0};
        size_t firstOffset{offset};
        while (offset < s.size()) {
            int c = jellyfish::mer_dna::code(s[offset]);
            if (jellyfish::mer_dna::not_dna(c)) {
                c = static_cast<int>('A');
                s[offset] = 'A';
            }
            mer.shift_left(c);
            if (jellyfish::mer_dna::not_dna(c)) {
                cmlen = 0;
                ++offset;
                continue;
            }
            if (++cmlen >= kl) {
                cmlen = kl;

                auto key = mer.get_bits(0, 2*kl);
                mers.insert(key);
            }
            ++offset;
          }


          /* decide the bucket */
          uint32_t maxBucketKey{*mers.begin()};
          uint32_t maxBucketValue{0};
          for (auto& k : mers) {
              auto cit = subMap.find(k);
              if (cit != subMap.end() and cit->second.size() > maxBucketValue) {
                  maxBucketValue = cit->second.size();
                  maxBucketKey = k;
              }
          }

          subMap[maxBucketKey].push_back(i);
    }

    std::vector<size_t> newInds;
    newInds.reserve(bucketReads.size());
    std::vector<size_t> sbSizes;
    for (auto& kv : subMap) {
        std::sort(kv.second.begin(), kv.second.end(),
                [&](const size_t& x, const size_t& y) -> bool {
                //const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
                    auto& a = bucketReads[x];
                    auto& b = bucketReads[y];
                    int cmp = std::get<1>(a) - std::get<1>(b);
                    if (cmp == 0) {
                        return std::get<0>(a) < std::get<0>(b);
                    } else {
                        return cmp < 0;
                    }
                });

        newInds.insert(newInds.end(),
                kv.second.begin(),
                kv.second.end());

        sbSizes.push_back(kv.second.size());
    }
    /*
        std::sort(newInds.begin(), newInds.end(),
                [&](const size_t& x, const size_t& y) -> bool {
                //const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
                    auto& a = bucketReads[x];
                    auto& b = bucketReads[y];
                    int cmp = std::get<1>(a) - std::get<1>(b);
                    if (cmp == 0) {
                        return std::get<0>(a) < std::get<0>(b);
                    } else {
                        return cmp < 0;
                    }
                });
    */

    std::vector<std::tuple<std::string, size_t>> newBucketReads;
    newBucketReads.reserve(bucketReads.size());

    for (auto it = newInds.begin(); it != newInds.end(); ++it) {
        newBucketReads.push_back(bucketReads[*it]);
    }
    std::swap(bucketReads, newBucketReads);
    //reorder_destructive(newInds.begin(), newInds.end(), bucketReads.begin());
    //reorder(bucketReads, newInds);

    /*
    std::cerr << "sub bucket size[ ";
    if (orig != std::get<0>(bucketReads[0])) {
        for (auto bs : sbSizes) {
            std::cerr << bs << " ";
        }
        std::cerr << " ]\n";
    }
    */
}

void reassignOnsies(
        std::vector<std::tuple<my_mer, std::string>>& onsies,
        CountMap& countMap,
        MapT& buckets) {

    for (auto& kr : onsies) {
        std::string s = std::get<1>(kr);
        size_t readLength{s.length()};
        size_t offset{0};
        size_t cmlen{0};
        size_t firstOffset{offset};
        size_t kl{15};
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

            if (++cmlen >= kl) {
                cmlen = kl;
                auto key = mer.get_bits(0, 2*kl);
                auto rkey = rmer.get_bits(0, 2*kl);

                auto it = merOffsetMap.find(key);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[key] = std::make_tuple(offset - (kl - 1), Direction::FORWARD);
                    if (s.substr(offset - (kl - 1), kl) != mer.to_str()) {
                        std::cerr << "putting " << mer.to_str() << " into the map, but I really have " << s.substr(offset - (kl-1), kl) << "\n";
                    }
                }
                it = merOffsetMap.find(rkey);
                if (it == merOffsetMap.end()) {
                    merOffsetMap[rkey] = std::make_tuple(readLength - offset - 1, Direction::REVERSE);
                }

            }
            ++offset;
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
            reverseComplement(s);
        }

        countMap[maxBucketKey]++;
        minmer.set_bits(0, 2*kl, maxBucketKey);
        firstOffset = std::get<0>(merOffsetMap[maxBucketKey]);

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

  po::options_description generic("mince options");
  generic.add_options()
      ("version,v", "print version information")
      ("help,h", "print help message")
      ("input,i", po::value<string>(), "input file; FASTA/Q if encoding, MINCE if decoding")
      ("output,o", po::value<string>(), "output file")
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

    const char* input[] = { vm["input"].as<string>().c_str() };

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
    unsigned int kl{15};
    my_mer::k(15);

    CountMap countMap;
    std::vector<MapT> maps(nb_threads);
    // Spawn of the read parsing threads
    for(int i = 0; i < nb_threads; ++i) {
        threads.push_back(std::thread(bucketReads<my_mer>, &parser, &total, std::ref(totReads), std::ref(maps[i]), std::ref(countMap)));
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
            size_t newOffset = readLength - kl - offset;
            // The string on which this read was originally bucketed
            my_mer bucketMer = std::get<0>(combined[0]);
            auto bucketString = bucketMer.to_str();

            // The split / swapped string
            std::string s(std::get<1>(combined[0]));
            std::string recon;

            // Decode the split / swapped string to recover the original
            if (offset == 0) {
                recon = bucketString + s;
            }
            else if (offset < s.size()) {
                recon = s.substr(newOffset) + bucketString + s.substr(0, newOffset);
            } else {
                recon = s + bucketString;
            }
            // Push back the original key and the reconstructed string
            // onsies.push_back(std::forward_as_tuple(k, recon));
        }
    }

    /*

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
    reassignOnsies(onsies, countMap, collatedMap);
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
    */

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


        /*
           if (reads.size() > 2) {
           approxTSPReorder(reads);
           }
           */

        std::sort(reads.begin(), reads.end(),
                [&](const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
                int cmp = std::get<1>(a) - std::get<1>(b);
                if (cmp == 0) {
                return std::get<0>(a) < std::get<0>(b);
                } else {
                return cmp < 0;
                }

                });
        /*
           if (reads.size() > 100) {
           subBucket(reads);
           }
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
            std::cout << k << "\n";
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
                    std::vector<uint8_t> bytes = twoBitEncode(rstr.substr(0, encReadLength - extraBasesPerRead));
                    outfile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
                    if (extraBasesPerRead > 0) {
                        sstream << rstr.substr(encReadLength - extraBasesPerRead);
                    }
                } catch (std::exception& e ) {
                    std::cerr << "EXCEPTION: " << e.what() << rstr << "\n";
                }
            }

            if (extraBasesPerRead > 0) {
                //std::string input = sstream.str().c_str();
                //char* output = new char[input.size()];
                //auto ret = divbwt(input.c_str(), output, NULL, input.size());
                //if (ret == -1 or ret == -2) { std::cerr << "BWT failed\n"; }
                //std::string outStr(output);
                std::vector<uint8_t> bonusBytes = twoBitEncode(sstream.str());
                outfile.write(reinterpret_cast<char*>(&bonusBytes[0]), sizeof(uint8_t) * bonusBytes.size());
                //delete [] output;
            }

            std::vector<uint8_t> offsets;
            offsets.reserve(sbsize);
            for (size_t j = 0; j <= sbsize; ++j) {
                auto& r = reads[i*256 + j];
                uint8_t l = static_cast<uint8_t>(std::get<1>(r));
                offsets.push_back(l);
                //outfile.write(reinterpret_cast<char*>(&l), sizeof(l));
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
            /*        std::cerr << "deltas [ ";
                      for (size_t j = 0; j <= sbsize; ++j) {
                      std::cerr << +deltas[j] << " ";
                      }
                      std::cerr << "]\n";
                      */
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
