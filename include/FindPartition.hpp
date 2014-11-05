#ifndef FIND_PARTITION_HPP
#define FIND_PARTITION_HPP

#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cmath>

namespace mince{
    class BucketInfo {
        public:
            BucketInfo() = default;
            BucketInfo(int32_t score_, uint32_t startIdx_, uint32_t stopIdx_,
                       int32_t prev_, uint32_t prefLen_) :
                score(score_), startIdx(startIdx_), stopIdx(stopIdx_), prev(prev_), prefLen(prefLen_) {}


            int32_t score;
            uint32_t startIdx, stopIdx;
            int32_t prev;
            uint32_t prefLen;

            uint32_t size() { return stopIdx - startIdx + 1; }
    };

    int32_t bucketScore(uint32_t coreLength, size_t i, size_t j);

    /**
     *  Given a vector of reads, all sharing a core string of length `bucketStringLength and a maximum bucket size,
     *  find the optimal partitioning of the vector into sub-ranges.
     */
    std::vector<BucketInfo> findPartition(std::vector<std::tuple<std::string, size_t>>& reads, uint32_t bucketStringLength,
            size_t readLength,
            size_t maxBucketSize=256);

}

#endif // FIND_PARTITION_HPP
