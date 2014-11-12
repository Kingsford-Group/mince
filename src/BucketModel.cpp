#include <cstddef>
#include <iostream>
#include <algorithm>


#include "BucketModel.hpp"



// Create a new kmer set, with a set backing to start
KmerSet::KmerSet() 
    : storage_(STO_SET) 
{
    s_ = new std::set<uint16_t>();
}

// free up the current storage usage
KmerSet::~KmerSet()
{
    switch (storage_) {
        case STO_SET:
            delete s_;
            break;

        case STO_BS: 
            delete bs_; 
            break;
    }
}

// add k to the set
void KmerSet::add(kmer_t k) 
{
    switch (storage_) {
        case STO_SET:
            s_->insert(k);
            if (s_->size() > 4096) {
                convert_to_bs();
            }
            break;

        case STO_BS: 
            (*bs_)[k] = 1;
            break;
    }
}


// return 1 if set contains k, 0 otherwise
int KmerSet::contains(kmer_t k) 
{
    switch (storage_) {
        case STO_SET:
            return (s_->find(k) != s_->end())?1:0;

        case STO_BS: 
            return (*bs_)[k];
    }
    return 0;
}

void KmerSet::convert_to_bs() 
{
    boost::dynamic_bitset<>* b = new boost::dynamic_bitset<>(65536);
    for (auto & k : *s_) {
        (*b)[k] = 1;
    }
    storage_ = STO_BS;
    delete s_;
    s_ = nullptr;
    bs_ = b;
}



BucketModel::BucketModel() : count_(0), kmers_() {
    //bloomFilt_(bloom_alloc(500,8), bloom_free) {
}

BucketModel::BucketModel(const BucketModel& o) : count_(o.count_.load()), kmers_(o.kmers_) { //bloomFilt_(nullptr, bloom_free) {
    /*
    bloomMutex_.lock();
    bloomFilt_.reset(bloom_copy(o.bloomFilt_.get()));
    bloomMutex_.unlock();
    */
    /*
    for (size_t i = 0; i < trimerCount_.size(); ++i) {
        trimerCount_[i] = o.trimerCount_[i].load();
    }
    */
}


void BucketModel::addCount(uint32_t inc) {
    count_ += inc;
}

void BucketModel::subCount(uint32_t inc) {
    count_ -= inc;
}

void BucketModel::addRead(std::string& s, uint8_t k) {
    size_t offset{0};
    size_t kl{k};
    size_t cmlen{0};

    Kmer mer;
    mer.polyT();

    bloomMutex_.lock();
    while (offset < s.size()) {
        // Get the code for the next base
        int c = jellyfish::mer_dna::code(s[offset]);

        // If it's not a valid DNA code
        if (jellyfish::mer_dna::not_dna(c)) {
            // Switch it to an 'A' in the mer
            c = jellyfish::mer_dna::code('A');
            mer.shift_left(c);
        } else { // Otherwise, base is legit
            mer.shift_left(c);
        }

        ++offset;
        // If we've read a full k-mer
        if (++cmlen >= kl) {
            cmlen = kl;
            kmer_t key = mer.get_bits(0, 2*kl);
            /*
            if (!bloom_get(bloomFilt_.get(), key)) {
                bloom_inc(bloomFilt_.get(), key);
            }
            */
            //trimerCount_[key] = 1;
            kmers_.add(key);
        }
    }
    bloomMutex_.unlock();

    count_ += 1;
}

double BucketModel::scoreOfReadRC(std::string& s, uint8_t k) {
    auto start = s.begin();
    auto stop = s.end();
    size_t kl{k};
    size_t cmlen{0};

    Kmer mer;
    mer.polyT();
    double ip{0.0};
    for (auto it = start; it != stop; ++it) {
        // Get the code for the next base
        int c = jellyfish::mer_dna::code(*it);

        // If it's not a valid DNA code
        if (jellyfish::mer_dna::not_dna(c)) {
            // Switch it to an 'A' in the mer
            c = jellyfish::mer_dna::code('A');
            mer.shift_right(c);
        } else { // Otherwise, base is legit
            mer.shift_right(jellyfish::mer_dna::complement(c));
        }

        // If we've read a full k-mer
        if (++cmlen >= kl) {
            cmlen = kl;
            kmer_t key = mer.get_bits(0, 2*kl);
            ip += kmers_.contains(key);
            //ip += trimerCount_[key];//bloom_get(bloomFilt_.get(), key);
        }
    }

    return ip;
}



double BucketModel::scoreOfRead(std::string& s, uint8_t k, bool rc) {
    Kmer::k(k);
    if (rc) {
        return scoreOfReadRC(s, k);
    }

    size_t offset{0};
    size_t kl{k};
    size_t cmlen{0};

    Kmer mer;
    mer.polyT();

    double ip{0.0};
    while (offset < s.size()) {
        // Get the code for the next base
        int c = jellyfish::mer_dna::code(s[offset]);

        // If it's not a valid DNA code
        if (jellyfish::mer_dna::not_dna(c)) {
            // Switch it to an 'A' in the mer
            c = jellyfish::mer_dna::code('A');
            mer.shift_left(c);
        } else { // Otherwise, base is legit
            mer.shift_left(c);
        }

        ++offset;
        // If we've read a full k-mer
        if (++cmlen >= kl) {
            cmlen = kl;
            kmer_t key = mer.get_bits(0, 2*kl);
            ip += kmers_.contains(key);
            //ip += trimerCount_[key];//bloom_get(bloomFilt_.get(), key);
        }

    }
    return ip;
}


/*
double BucketModel::scoreOfRead(std::vector<double>& featureVec) {
    if (count_.load() == 0) { return 0.0; }
    double ip{0.0};
    uint64_t n1{0};
    uint64_t n2{0};
    for (size_t i = 0; i < trimerCount_.size(); ++i) {
        ip += trimerCount_[i] * featureVec[i];
        n1 += trimerCount_[i];
        n2 += featureVec[i];
    }
    if (!(ip >0)) {
        std::cerr << "\ntrimerCount.size() = " << trimerCount_.size() << "\n"
                  << "featureVec.size() = " << featureVec.size() << "\n"
                  << "count = " << count_.load() << "\n"
                  << "n1 = " << n1 << "\nn2 = " << n2 << "\n";
        std::cerr << "trimerCount = [";
        for (size_t i = 0; i < trimerCount_.size(); ++i) {
            if (trimerCount_[i] > 0) {
                std::cerr << i << " : " << trimerCount_[i] << ", ";
            }
        }
        std::cerr << "]\n";
        std::cerr << "featureVec = [";
        for (size_t i = 0; i < featureVec.size(); ++i) {
            if (featureVec[i] > 0) {
                std::cerr << i << " : " << featureVec[i] << ", ";
            }
        }
        std::cerr << "]\n";
        std::cerr << "SHOULD AT LEAST SHARE THE KEY!!!!\n";
        std::exit(1);
    }
    return ip;
    //return count_.load();
    //return count_.load() * ( ip / (std::sqrt(n1*n1) * std::sqrt(n2*n2)));
}

*/