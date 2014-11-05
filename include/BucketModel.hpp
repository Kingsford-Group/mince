#ifndef __BUCKET_MODEL__
#define __BUCKET_MODEL__

#include <atomic>
#include <vector>
#include <mutex>
#include <memory>

#include <boost/dynamic_bitset.hpp>

#include <jellyfish/mer_dna.hpp>

extern "C" {
    #include "bloom.h"
}

using Kmer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 3>;
using kmer_t = uint64_t;

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
        boost::dynamic_bitset<> trimerCount_;
        std::mutex bloomMutex_;
        //std::unique_ptr<bloom_t, void(*)(bloom_t*)> bloomFilt_;
};

#endif // __BUCKET_MODEL__
