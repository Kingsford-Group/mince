#ifndef _MINCE_UTILS_HPP_
#define _MINCE_UTILS_HPP_

#include <vector>
#include <fstream>

#include <boost/dynamic_bitset.hpp>

namespace mince {
namespace utils {
    std::string twoBitDecode(uint8_t* read, size_t readLength) {
        boost::dynamic_bitset<uint8_t> bitmap(2*readLength, 0);
        size_t effectiveReadLength = bitmap.num_blocks();
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
}
}

#endif // _MINCE_UTILS_HPP_

