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

#include "MinceUtils.hpp"
#include "FindPartition.hpp"

namespace mince {
    int32_t bucketScore(uint32_t coreLength, size_t i, size_t j) {
        int32_t coreStringBytes = std::ceil(static_cast<float>(coreLength) / 4.0);
        return (-coreStringBytes - 2) + ((i - j + 1) * coreStringBytes);
    }

    /**
      *  Given a vector of reads, all sharing a core string of length `bucketStringLength and a maximum bucket size,
      *  find the optimal partitioning of the vector into sub-ranges.
      */
    std::vector<BucketInfo> findPartition(std::vector<std::tuple<std::string, size_t>>& reads, uint32_t bucketStringLength,
                       size_t readLength,
                       size_t maxBucketSize) {
        // Sort the reads, or assume they are already sorted lexicographically?
        std::vector<BucketInfo> bucketInfo(reads.size());
        std::vector<int32_t> prevBucketInd(reads.size());
        std::vector<int64_t> opt(reads.size());

        // For each read
        for (size_t i = 0; i < reads.size(); ++i) {
            std::string& ri = std::get<0>(reads[i]);
            auto offseti = std::get<1>(reads[i]);

            bucketInfo[i].score = -2;//std::numeric_limits<int32_t>::max();
            bucketInfo[i].startIdx = i;
            bucketInfo[i].stopIdx = i;
            bucketInfo[i].prev = i - 1;
            bucketInfo[i].prefLen = 0;

            // j starts at either i - bucketSize or 0.
            int64_t jstart = (i > maxBucketSize) ? i - maxBucketSize : 0;

            uint32_t minLenY = std::numeric_limits<uint32_t>::max();

            // For every j less than the current i, evaluate the solution in which
            // [i,j] is a bucket;
            for (int64_t j = i - 1; j >= jstart; --j){ //jstart; j < i; ++j) {
                std::string& rj = std::get<0>(reads[j]);
                auto offsetj = std::get<1>(reads[j]);
                uint32_t lenY = std::min(ri.length() - offseti,
                                         rj.length() - offsetj);

                lenY = utils::lcp(ri, rj, lenY);
                minLenY = std::min(lenY, minLenY);
                lenY = minLenY;

                if (ri.substr(0, lenY) != rj.substr(0, lenY)) {
                    std::cerr << "Found non-matching LCPs: "
                              << ri.substr(0, lenY) << '\n'
                              << rj.substr(0, lenY) << "!\n";
                    std::exit(1);
                }

                if (lenY > ri.length() - offseti or
                    lenY > rj.length() - offsetj) {
                    std::cerr << "ERROR!\n";
                    std::cerr << "found prefix " << ri.substr(0, lenY) << "\n";
                    std::cerr << "that is longer than min("
                              << ri.length() - offseti << ", "
                              << rj.length() - offsetj << ")\n";
                    std::exit(1);
                }
                if (offseti > ri.substr(lenY).length()) {
                    std::cerr << "offseti = " << +offseti << " > "
                              << "read len = " << ri.substr(lenY).length() << "\n";
                    std::exit(1);
                }
                if (offsetj > rj.substr(lenY).length()) {
                    std::cerr << "offsetj = " << +offsetj << " > "
                              << "read len = " << rj.substr(lenY).length() << "\n";
                    std::exit(1);
                }

                /*
                for (size_t k = jstart; k < j; ++k) {
                    auto rk = std::get<0>(reads[k]);
                    auto offsetk = std::get<1>(reads[k]);
                    if (offsetk > rk.substr(lenY).length()) {
                        std::cerr << "offsetk = " << +offsetk << " > "
                            << "read len = " << rk.substr(lenY).length() << "\n";
                        std::exit(1);
                    }
                }
                */

                uint32_t coreLength = bucketStringLength + lenY;
                auto prevOpt = (j > 0) ? bucketInfo[j-1].score : 0;
                int32_t score = bucketScore(coreLength, i, j) + prevOpt;
                if (score >  bucketInfo[i].score) {
                    bucketInfo[i].score = score;
                    bucketInfo[i].startIdx = j; bucketInfo[i].stopIdx = i;
                    bucketInfo[i].prefLen = lenY;
                }
            }
        }


        /*
        std::cerr << "opt = [ ";
        for (auto o : opt) {
            std::cerr << o << ' ';
        }
        std::cerr << "]\n";

        std::cerr << "prevBucketInd = [ ";
        for (auto o : prevBucketInd) {
            std::cerr << o << ' ';
        }
        std::cerr << "]\n";
        std::exit(1);
*/


        // Backtrace
        std::vector<BucketInfo> selectedBuckets;
        selectedBuckets.push_back(bucketInfo.back());
        while (selectedBuckets.back().startIdx > 0) {
            int64_t prevBucket = selectedBuckets.back().startIdx - 1;
            selectedBuckets.push_back(bucketInfo[prevBucket]);
        }
        /*
        bool done{false};
        while(!done) {
            if (prevPtr > 0) {
                uint32_t bucketSize = currEnd - prevPtr;
                selectedBuckets.emplace_back(
                        bucketInfo[currEnd].score,
                        bucketSize,
                        bucketInfo[currEnd].prev,
                        bucketInfo[currEnd].prefLen);
                currEnd = prevPtr;
                prevPtr = bucketInfo[prevPtr].prev;
            } else {
                uint32_t bucketSize = currEnd;
                selectedBuckets.emplace_back(
                        bucketInfo[currEnd].score,
                        bucketSize,
                        bucketInfo[currEnd].prev,
                        bucketInfo[currEnd].prefLen);
                done = true;
            }
        }
        */
        std::reverse(selectedBuckets.begin(), selectedBuckets.end());
        return selectedBuckets;
    }
}
