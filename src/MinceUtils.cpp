#include <exception>
#include <sstream>
#include "MinceUtils.hpp"

namespace mince {
    namespace utils {
        // Return the length of the shared prefix of strings a and b
        uint32_t lcp(std::string& a, std::string& b, uint32_t maxLen) {
            size_t l = a.length();
            uint32_t pl = 0;
            while (a[pl] == b[pl] and pl < l and pl < maxLen) {
                ++pl;
            }
            return pl;
        }

	std::string unpermute(std::string& permS, std::string& key, size_t offset) {
		size_t l = permS.length();
		return permS.substr(l - offset) + key + permS.substr(0, l - offset);
	}

	std::string permute(std::string& s, size_t offset, size_t kl) {
		return s.substr(offset + kl) + s.substr(0, offset);
	}

        std::string twoBitDecode(const uint8_t* read, size_t readLength) {
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
            bitmap.reset();
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

        /**
         * This function parses the library format string that specifies the format in which
         * the reads are to be expected.
         */
        LibraryFormat parseLibraryFormatString(std::string& fmt) {
            using std::vector;
            using std::string;
            using std::map;
            using std::stringstream;

            map<char, ReadOrientation> orientationType = {{'S', ReadOrientation::SAME},
                {'A', ReadOrientation::AWAY},
                {'T', ReadOrientation::TOWARD},
                {'I', ReadOrientation::NONE}};
            map<string, ReadStrandedness> strandType = {{"FF", ReadStrandedness::S},
                {"RR", ReadStrandedness::A},
                {"FR", ReadStrandedness::SA},
                {"RF", ReadStrandedness::AS},
                {"UU", ReadStrandedness::U},
                {"F", ReadStrandedness::S},
                {"R", ReadStrandedness::A},
                {"U", ReadStrandedness::U}};

            // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
            // first convert the string to upper-case
            for (auto& c : fmt) { c = std::toupper(c); }

            auto orientationIt = orientationType.find(fmt.front());
            auto strandednessIt = strandType.find(fmt.substr(1));

            if (orientationIt == orientationType.end()) {
                stringstream errstr;
                errstr << "unknown orientation type: " << fmt.front();
                throw std::invalid_argument(errstr.str());
            }
            if (strandednessIt == strandType.end()) {
                stringstream errstr;
                errstr << "unknown strand type: " << fmt.substr(1);
                throw std::invalid_argument(errstr.str());
            }

            // if we recognize the orientation & strand
            auto orientation = orientationIt->second;
            auto strandedness = strandednessIt->second;

            ReadType rt;
            if (orientation == ReadOrientation::NONE) {
                rt = ReadType::SINGLE_END;
            } else {
                rt = ReadType::PAIRED_END;
            }

            return LibraryFormat(rt, orientation, strandedness);
        }

        /**
         * Parses a set of __ordered__ command line options and extracts the relevant
         * read libraries from them.
         */
        std::vector<ReadLibrary> extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
            // The current (default) format for paired end data
            LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U);
            // The current (default) format for single end data
            LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);

            std::vector<ReadLibrary> peLibs{ReadLibrary(peFormat)};
            std::vector<ReadLibrary> seLibs{ReadLibrary(seFormat)};
            for (auto& opt : orderedOptions.options) {
                // Update the library type
                if (opt.string_key == "libtype") {
                    auto libFmt = parseLibraryFormatString(opt.value[0]);
                    if (libFmt.type == ReadType::PAIRED_END) {
                        peFormat = libFmt;
                        peLibs.emplace_back(libFmt);
                    } else {
                        seFormat = libFmt;
                        seLibs.emplace_back(libFmt);
                    }
                }
                if (opt.string_key == "mates1") {
                    std::cerr << "mates1\n";
                    peLibs.back().addMates1(opt.value);
                }
                if (opt.string_key == "mates2") {
                    std::cerr << "mates2\n";
                    peLibs.back().addMates2(opt.value);
                }
                if (opt.string_key == "unmated_reads") {
                    std::cerr << "unmated\n";
                    seLibs.back().addUnmated(opt.value);
                }
            }

            std::vector<ReadLibrary> libs;
            libs.reserve(peLibs.size() + seLibs.size());
            for (auto& lib : boost::range::join(seLibs, peLibs)) {
                if (lib.format().type == ReadType::SINGLE_END) {
                    if (lib.unmated().size() == 0) {
                        std::cerr << "skipping single-end library w/ no reads\n";
                        continue;
                    }
                } else if (lib.format().type == ReadType::PAIRED_END) {
                    if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
                        std::cerr << "skipping paired-end library w/ no reads\n";
                        continue;
                    }
                }
                libs.push_back(lib);
            }
            std::cerr << "there are " << libs.size() << " libs\n";
            return libs;
        }

    } // namespace utils
} // namespace mince
