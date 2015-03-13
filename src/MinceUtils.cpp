#include <exception>
#include <sstream>

#include "jellyfish/mer_dna.hpp"

#include "MinceUtils.hpp"

namespace mince {
    namespace utils {
        using Trimer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 2>;

        std::vector<double> trimerVectorRC(std::string& s) {
            auto start = s.begin();
            auto stop = s.end();
            size_t kl{6};
            size_t cmlen{0};
            std::vector<double> fvec(4096, 0.0);

            Trimer mer;
            mer.polyT();

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
                    auto key = mer.get_bits(0, 2*kl);
                    fvec[key] += 1.0;
                }
            }

            return fvec;
        }

        // Return a feature vector with the number of occurences of each
        // trimer in the read.
        std::vector<double> trimerVector(std::string& s, bool rc){
            Trimer::k(6);
            if (rc) { return trimerVectorRC(s); }

            size_t offset{0};
            size_t kl{6};
            size_t cmlen{0};
            std::vector<double> fvec(4096, 0.0);

            Trimer mer;
            mer.polyT();

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
                    auto key = mer.get_bits(0, 2*kl);
                    fvec[key] += 1.0;
                }
            }

            return fvec;
        }

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
		//size_t l = permS.length();
        //return permS.substr(0, offset) + key + permS.substr(offset);
        size_t l = permS.length();
		return permS.substr(l - offset) + key + permS.substr(0, l - offset);
	}

	std::string permute(std::string& s, size_t offset, size_t kl) {
        //return s.substr(0, offset) + s.substr(offset+kl);
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
	LibraryFormat parseLibraryFormatStringNew(std::string& fmt) {
		using std::vector;
		using std::string;
		using std::map;
		using std::stringstream;

		map<string, LibraryFormat> formatMap = {
			{"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
			{"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA)},
			{"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS)},
			{"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
			{"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA)},
			{"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS)},
			{"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
			{"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
			{"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)},
			{"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
			{"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
			{"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)}};

		// inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
		// first convert the string to upper-case
		for (auto& c : fmt) { c = std::toupper(c); }


		auto libFmtIt = formatMap.find(fmt);

		if (libFmtIt == formatMap.end()) {
			stringstream errstr;
			errstr << "unknown library format string : " << fmt;
			throw std::invalid_argument(errstr.str());
		}

		return libFmtIt->second;
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
				auto libFmt = parseLibraryFormatStringNew(opt.value[0]);
				if (libFmt.type == ReadType::PAIRED_END) {
					peFormat = libFmt;
					peLibs.emplace_back(libFmt);
				} else {
					seFormat = libFmt;
					seLibs.emplace_back(libFmt);
				}
			}
			if (opt.string_key == "mates1") {
				peLibs.back().addMates1(opt.value);
			}
			if (opt.string_key == "mates2") {
				peLibs.back().addMates2(opt.value);
			}
			if (opt.string_key == "unmated_reads") {
				seLibs.back().addUnmated(opt.value);
			}
		}

		std::vector<ReadLibrary> libs;
		libs.reserve(peLibs.size() + seLibs.size());
		for (auto& lib : boost::range::join(seLibs, peLibs)) {
			if (lib.format().type == ReadType::SINGLE_END) {
				if (lib.unmated().size() == 0) {
					// Didn't use default single end library type
					continue;
				}
			} else if (lib.format().type == ReadType::PAIRED_END) {
				if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
					// Didn't use default paired-end library type
					continue;
				}
			}
			libs.push_back(lib);
		}
		size_t numLibs = libs.size();
		std::cerr << "there " << ((numLibs > 1) ? "are " : "is ") << libs.size() << ((numLibs > 1) ? " libs\n" : " lib\n");
		return libs;
	}


    } // namespace utils
} // namespace mince
