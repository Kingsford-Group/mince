#ifndef __BUCKET_MODEL__
#define __BUCKET_MODEL__

#include <atomic>
#include <vector>
#include <mutex>
#include <memory>
#include <set>

#include <boost/dynamic_bitset.hpp>

#include <jellyfish/mer_dna.hpp>

extern "C" {
    #include "bloom.h"
}

using Kmer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 3>;
using kmer_t = uint64_t;

class KmerSet {
public:
   KmerSet();
   KmerSet(const KmerSet& o);
   ~KmerSet();
   void add(kmer_t k);
   int contains(kmer_t k);
   void operator=(const KmerSet& o);

protected:
   /*void add_set(kmer_t k);
   bool contains_set(kmer_t k);
   void add_bs(kmer_t k);
   void contains_bs(kmer_t k); */
   void convert_to_bs();
          
private: 
       std::set<uint16_t>* s_;
       boost::dynamic_bitset<>* bs_;
       std::vector<uint16_t>* v_;
   enum { STO_SET, STO_BS, STO_VEC } storage_;
};       

class BucketModel {
    public:
        BucketModel();
        BucketModel(const BucketModel& o);

        void addCount(uint32_t inc=1);
        void subCount(uint32_t inc=1);
        double scoreOfReadRC(std::string& s, uint8_t k);
        double scoreOfRead(std::string& s, uint8_t k, bool rc);
        void addRead(std::string& s, uint8_t k);

    private:
        std::atomic<uint64_t> count_;
        //std::vector<std::atomic<uint32_t>> trimerCount_;
        //boost::dynamic_bitset<> trimerCount_;
        KmerSet kmers_;
        std::mutex bloomMutex_;
        //std::unique_ptr<bloom_t, void(*)(bloom_t*)> bloomFilt_;
};

#endif // __BUCKET_MODEL__
