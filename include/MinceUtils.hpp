#ifndef _MINCE_UTILS_HPP_
#define _MINCE_UTILS_HPP_

#include <vector>
#include <fstream>
#include <cstdio>
#include <cctype>
#include <locale>

#include <boost/range/join.hpp>
#include <boost/program_options.hpp>
#include <boost/dynamic_bitset.hpp>

#include "ReadLibrary.hpp"

namespace mince {
namespace utils {

    // Return the length of the shared prefix of strings a and b
    uint32_t lcp(std::string& a, std::string& b, uint32_t maxLen);

    std::string twoBitDecode(const uint8_t* read, size_t readLength);

    std::vector<uint8_t> twoBitEncode(const std::string& str);

    void reverseComplement(std::string& s);

    /**
     * This function parses the library format string that specifies the format in which
     * the reads are to be expected.
     */
    LibraryFormat parseLibraryFormatString(std::string& fmt);

    /**
      * Parses a set of __ordered__ command line options and extracts the relevant
      * read libraries from them.
      */
    std::vector<ReadLibrary> extractReadLibraries(boost::program_options::parsed_options& orderedOptions);

}
}

#endif // _MINCE_UTILS_HPP_
