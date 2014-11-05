#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <functional>

#include <boost/container/flat_set.hpp>
#include <boost/heap/priority_queue.hpp>
#include <boost/dynamic_bitset.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

struct CountList {
    uint32_t count;
    std::set<uint32_t> rlist;
};

using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
using MapT = std::unordered_map<uint32_t, std::vector<uint32_t>>;

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


// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename MerT>
void add_sizes(sequence_parser* parser, std::atomic<uint64_t>* total,
               std::atomic<uint64_t>& totReads,
               MapT& buckets,
               std::vector<std::tuple<size_t, std::string>>& reads) {
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

      size_t readID{0};
      for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
          readID = totReads++;
          if (readID % 500000 == 0) {
              std::cerr << "\r\rprocessed " << totReads << " reads";
          }
          std::string s(j->data[i].seq);

          cmlen = 0;
          MerT mer(kl);
          MerT minmer(kl);
          MerT rmer(kl);
          minmer.polyT();
          mer.polyT();
          rmer.polyT();

          std::set<MerT> mers;
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
                if (mers.find(mer) == mers.end()) {
                    size_t key{mer.get_bits(0, 2*kl)};
                    auto loc = buckets.find(key);
                    if (loc == buckets.end()) {
                        buckets[key]= {readID};
                    } else {
                        loc->second.push_back(readID);
                    }
                    mers.insert(mer);
                }
            }
            ++offset;
          }
          mers.clear();
          reads.emplace_back(make_tuple(readID, j->data[i].seq));
        /*
          if (firstOffset + kl < j->data[i].seq.size()) {
              std::string reord(j->data[i].seq.substr(firstOffset + kl) + j->data[i].seq.substr(0, firstOffset));
              buckets[minmer].push_back(make_tuple(reord, firstOffset));
              //buckets[minmer].push_back(make_tuple(j->data[i].seq, firstOffset));

              //auto pos = reord.find('*');
              //std::string recon = reord.substr(pos + 1) + minmer.to_str() + reord.substr(0, pos);
              //if (recon != j->data[i].seq) {
                  //std::cerr << "recon: " << recon << "\n";
                  //std::cerr << "orig:  " << j->data[i].seq << "\n";
              //}
          } else {

              //buckets[minmer].push_back(make_tuple(j->data[i].seq, firstOffset));
              buckets[minmer].push_back(make_tuple(j->data[i].seq.substr(0, firstOffset), firstOffset));
          }
          */
          count += j->data[i].seq.size();        // Add up the size of the sequence
      }
  }

  *total += count;
}

template <typename MerT>
void mergeMaps(std::vector<MapT>& maps, MapT& outMap) {
    size_t i{0};
    size_t totOtherKeys{0};
    for (auto& m : maps) {
        totOtherKeys += m.size();
    }
    for (auto& m : maps) {
        for (auto& kv : m ) {
            auto it = outMap.find(kv.first);
            if (it == outMap.end()) {
               outMap[kv.first] = kv.second;
               it->second.clear();
            } else {
               it->second.reserve(it->second.size() + kv.second.size());
               it->second.insert(
                       it->second.end(),
                       kv.second.begin(),
                       kv.second.end());
               it->second.clear();
            }
            ++i;
            if (i % 100000 == 0) {
                std::cerr << "\r\ri = " << i << " / " << totOtherKeys;
            }
        }
        m.clear();
    }

    size_t numKeys{outMap.size()};
    std::cerr << "collated " << numKeys << " keys\n";
    /*
    size_t kn{0};
    for (auto& kv : outMap) {
        auto& k = kv.first;
        size_t tot{0};
        std::vector<std::tuple<MapT&, MapT::iterator>> its;
        for (auto& m : maps) {
            auto kit = m.find(k);
            if (kit != m.end()) {
                its.push_back(make_tuple(std::ref(m), kit));
                tot += kit->second.size();
            }
        }

        for (auto& itTup : its) {
                auto& outV = kv.second;
                outV.reserve(tot);

                auto beg = std::get<1>(itTup)->second.begin();
                auto end = std::get<1>(itTup)->second.end();
                outV.insert(
                        outV.end(),
                        beg,
                        end);
                std::get<0>(itTup).erase(std::get<1>(itTup));
        }
        ++kn;
        if (kn % 1000 == 0) {
            std::cerr << "\r\rmerged " << kn << " / " << numKeys << " keys";
        }
    }
    */
}


void writeBucket(std::ofstream& outfile, my_mer& bucketMer, std::vector<std::tuple<std::string, uint8_t>>& reads) {
    size_t kl{15};
    char sep = '$';
    std::sort(reads.begin(), reads.end(),
            [&](const std::tuple<std::string, size_t>& a, const std::tuple<std::string, size_t>& b) -> bool {
            int cmp = std::get<1>(a) - std::get<1>(b);
            if (cmp == 0) {
            return std::get<0>(a) < std::get<0>(b);
            } else {
            return cmp < 0;
            }

            });

    outfile.write(reinterpret_cast<char*>(&sep), sizeof(sep));
    auto key = bucketMer.get_bits(0, 2*kl);
    outfile.write(reinterpret_cast<char*>(&key), sizeof(key) );

    for (auto& r : reads) {
        std::string rstr = std::get<0>(r);
        uint8_t l = static_cast<uint8_t>(std::get<1>(r));

        try {
            std::vector<uint8_t> bytes = twoBitEncode(rstr);
            outfile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
        } catch (std::exception& e ) {
            std::cerr << rstr << "\n";
        }
        outfile.write(reinterpret_cast<char*>(&l), sizeof(l));
    }
}

int main(int argc, char *argv[]) {
  const int nb_threads      = 20;
  const int concurrent_file = 1;   // Number of files to read simultaneously
  const int max_read_group  = 100; // Number of reads in each "job" group

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

  /*
  bfs::path thashFile = std::string(argv[2]);

  // Read in the Jellyfish hash of the transcripts
  std::ifstream jellyfishDB(thashFile.c_str());
  if (!jellyfishDB.good()) {
      std::cerr << "Couldn't open the Jellyfish hash [" << thashFile << "] quitting\n";
      std::exit(-1);
  }
  jellyfish::file_header header;
  header.read(jellyfishDB);

  std::cerr << "transcript hash size is " << header.size() << "\n";
  size_t nkeys = header.size();
  // Since JF2, key_len() is in terms of bits so merLen = keyLen / 2;
  size_t merLen = header.key_len() / 2;
  my_mer::k(merLen);

  size_t tot{0};

  const jellyfish::RectangularBinaryMatrix m(header.matrix());
  size_t sizeMask{header.size() - 1};

  auto merHasher = [&](const my_mer& mer) -> size_t {
      //auto key = mer.get_bits(0, jellyfish::mer_dna::k());
      return m.times(mer) & sizeMask;
  };

  std::unordered_map<my_mer, std::atomic<size_t>, decltype(merHasher)> merIDMap(10, merHasher);

  if (!header.format().compare(binary_dumper::format)) {
      binary_reader reader(jellyfishDB, &header);
      while ( reader.next() ) {
          tot += reader.val();
          merIDMap[reader.key()] = reader.val();
          ++distinct;
          if (distinct % 1000000 == 0) {
              std::cerr << "\r\rinserted " << distinct << " keys into hash";
          }
      }
  }
  */

  unsigned int kl{15};
  my_mer::k(kl);

  //std::vector<std::map<my_mer, std::vector<std::tuple<std::string, size_t>>>> maps(nb_threads);
  std::vector<MapT> maps(nb_threads);
  std::vector<std::vector<std::tuple<size_t, std::string>>> readList(nb_threads);
  for(int i = 0; i < nb_threads; ++i) {
      threads.push_back(std::thread(add_sizes<my_mer>, &parser, &total, std::ref(totReads), std::ref(maps[i]), std::ref(readList[i])));
  }

  for(int i = 0; i < nb_threads; ++i)
    threads[i].join();

  std::cerr << "\n\ncounted " << totReads << " reads\n";
  std::cerr << "joining all reads . . . ";
  std::vector<std::tuple<size_t, std::string>> reads;
  for (int i = 0; i < nb_threads; ++i) {
      reads.insert(
              reads.end(),
              std::make_move_iterator(readList[i].begin()),
              std::make_move_iterator(readList[i].end())
              );
  }
  std::cerr << "done\n";

  // sort the read sequences so that they are in order
  std::cerr << "sorting reads . . . ";
  std::sort(reads.begin(), reads.end());
  std::cerr << "done\n";

  std::cerr << "joining all maps . . . ";
  MapT outMap;
  mergeMaps<MapT>(maps, outMap);
  std::cerr << "done\n";
  struct MerCount {
    uint32_t mer;
    uint32_t count;
  };

  struct MerCountComparator {
      bool operator() (const MerCount& m1,const MerCount& m2) const {
        return m1.count < m2.count;
      }
  };

  boost::heap::priority_queue<MerCount, boost::heap::compare<MerCountComparator>> pq;
  // populate the queue
  std::cerr << "populating the priority queue . . . ";
  for (auto& kv : outMap) { pq.push({kv.first, kv.second.size()}); }
  std::cerr << "done\n";

  size_t count{totReads};
  while (count > 0) {
    // pick the biggest bucket
    auto tk = pq.top();
    pq.pop();
    size_t i{0};
    auto mit = outMap.find(tk.mer);
    while ( (mit == outMap.end()) or mit->second.size() < pq.top().count) {
        if (mit != outMap.end()) {
            tk.count = mit->second.size();//outMap[tk.mer].size();
            pq.push(tk);
       }
        tk = pq.top();
        pq.pop();
        mit = outMap.find(tk.mer);
    }

    size_t bucketSize = mit->second.size();
    std::vector<std::tuple<std::string, uint8_t>> bucketReads;
    // For each read in this bucket, update all k-mers that overlap it
    for (auto rID : mit->second) {
          auto& s = std::get<1>(reads[rID]);
          size_t cmlen = 0;
          my_mer mer(kl);
          mer.polyT();

          std::set<uint32_t> mers;
          size_t offset{0};
          size_t firstOffset{offset};
          while (offset < s.size()) {
              int c = jellyfish::mer_dna::code(s[offset]);
              mer.shift_left(c);
              if (jellyfish::mer_dna::not_dna(c)) {
                  cmlen = 0;
                  ++offset;
                  continue;
              }
              if (++cmlen >= kl) {
                  cmlen = kl;
                  mers.insert(mer.get_bits(0, 2*kl));
                  if (firstOffset == 0 and mer.get_bits(0, 2*kl) == tk.mer) {
                      firstOffset = offset-(kl - 1);
                  }
              }
              ++offset;
          }
        for (auto& km : mers) {
            if (km == tk.mer) { continue; }
            auto mit = outMap.find(km);
            if (mit != outMap.end()) {
                auto& v = mit->second;
                auto it = std::find(v.begin(), v.end(), rID);
                if (it != v.end()) {
                    v.erase(it);
                }
                if (v.size() == 0) { outMap.erase(mit); }
            }
        }

        if (firstOffset + kl < s.size()) {
            std::string reord(s.substr(firstOffset + kl) + s.substr(0, firstOffset));
            bucketReads.push_back(make_tuple(reord, firstOffset));
        } else {
            bucketReads.push_back(make_tuple(s.substr(0, firstOffset), firstOffset));
        }

    }
    my_mer bucketMer(kl);
    bucketMer.set_bits(0, 2*kl, tk.mer);
    writeBucket(outfile, bucketMer, bucketReads);
    outMap[tk.mer].clear();
    count -= bucketSize;
    std::cerr << "remining reads = " << count << "\n";
  }

  outfile.close();
  /*
  std::set<my_mer> keys;
  for (auto& m : maps) {
      for (auto& kv : m ) {
        keys.insert(kv.first);
      }
  }

  char sep = '$';
  //std::cerr << "Num buckets = " << buckets.size() << "\n";
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
    //outfile << sep << kv.first.to_str();


    outfile.write(reinterpret_cast<char*>(&sep), sizeof(sep));
    auto key = k.get_bits(0, 2*kl);
    outfile.write(reinterpret_cast<char*>(&key), sizeof(key) );

    for (auto& r : reads) {
        std::string rstr = std::get<0>(r);
        uint8_t l = static_cast<uint8_t>(std::get<1>(r));

        try {
            std::vector<uint8_t> bytes = twoBitEncode(rstr);
            outfile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(uint8_t) * bytes.size());
         } catch (std::exception& e ) {
            std::cerr << rstr << "\n";
        }
        outfile.write(reinterpret_cast<char*>(&l), sizeof(l));

        //outfile << std::get<0>(r);
        //outfile << std::get<1>(r);
    }
  }
  outfile.close();
*/

  std::cerr << "Total bases: " << total << "\n";

  return 0;
}
