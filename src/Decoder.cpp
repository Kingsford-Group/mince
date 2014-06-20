#include <thread>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <functional>
#include <cmath>
#include <mutex>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>

#include "g2logworker.h"
#include "g2log.h"
#include "bitfile.h"
#include "MinceUtils.hpp"
#include "Decoder.hpp"

void Decoder::decode(std::string& ifname, std::string& ofname) {
    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
    std::ifstream seqFile;
    std::ifstream offFile;
    std::ifstream nsFile;
    bit_file_c flipFile;

    namespace bfs = boost::filesystem;
 
    bfs::path seqPath = ifname+".seqs";
    bfs::path offPath = ifname+".offs";
    bfs::path nsPath = ifname+".nlocs";
    bfs::path flipPath = ifname+".flips";

    if (!bfs::exists(seqPath) or !bfs::exists(offPath) or !bfs::exists(nsPath)) {
	std::string errstr = "Couldn't find mince sequence [.seqs] or offset [.offs] file or n's file [.nlocs]";
	std::cerr << errstr << "\n";
	LOG(FATAL) << errstr;
    }

    seqFile.open(seqPath.string(), std::ios::in | std::ios::binary);
    offFile.open(offPath.string(), std::ios::in | std::ios::binary);
    nsFile.open(nsPath.string(), std::ios::in | std::ios::binary);

    bool haveFlipFile{true};
    if (!bfs::exists(flipPath)) {
	LOG(INFO) << "no flip file found";
	haveFlipFile = false;
    } else {
	LOG(INFO) << "opened flip file";
	flipFile.Open(flipPath.c_str(), BF_READ);
    }

    uint8_t rl{0};
    uint8_t kl{0};
    size_t numOnsies{0};

    seqFile.read(reinterpret_cast<char*>(&kl), sizeof(kl));
    seqFile.read(reinterpret_cast<char*>(&rl), sizeof(rl));
    seqFile.read(reinterpret_cast<char*>(&numOnsies), sizeof(numOnsies));

    std::ofstream ofile;
    ofile.open(ofname);

    size_t readLength = rl;
    uint64_t bucketMer{0};
    std::string bucketString('X', kl);
    my_mer bucketKmer;
    my_mer::k(kl);

    size_t effectiveReadLength = std::ceil((readLength) / 4.0);

    size_t i{0};
    char* onsieRead = new char[effectiveReadLength];

    class NPatch {
	public:
	void getNext(std::ifstream& ifile) {
		ifile.read(reinterpret_cast<char*>(&id), sizeof(id));
		if (ifile.good()) {
			uint8_t nn{0};
			ifile.read(reinterpret_cast<char*>(&nn), sizeof(nn));
			nlocs.resize(nn, 0);
			ifile.read(reinterpret_cast<char*>(&nlocs[0]), nn * sizeof(nlocs[0]));
		} else {
			id = std::numeric_limits<uint32_t>::max();
		}
	}
	void apply(std::string& s) {
		for (auto p : nlocs) { 
			s[p] = 'N';
		}
	}
	uint32_t id;
	std::vector<uint8_t> nlocs;
    };

    std::cerr << "kl = " << +kl << "\n";
    std::cerr << "readLength = " << readLength << "\n";
    std::cerr << "num onsies = " << numOnsies << "\n";

    LOG(INFO) << "bucket string length = " << +kl;
    LOG(INFO) << "effective read length = " << effectiveReadLength;
    LOG(INFO) << "reading " << numOnsies << " onsies"; 

    NPatch np;
    np.getNext(nsFile);
    // decode all of the onsies
    for (size_t j = 0; j < numOnsies; ++j) {
        seqFile.read(reinterpret_cast<char*>(onsieRead), effectiveReadLength);
        std::string recon = mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(onsieRead), readLength);
	if (np.id == j) { 
		np.apply(recon); 
		np.getNext(nsFile);
	}
	bool doFlip = haveFlipFile ? flipFile.GetBit() : 0;
	if (doFlip) { mince::utils::reverseComplement(recon); }
        ofile << ">" << j << "\n" << recon << std::endl;
    }

    std::cerr << "done with onsies\n";
    delete onsieRead;

    uint32_t maxBucketSize{256};

    effectiveReadLength = std::ceil((readLength - kl) / 4.0);
    LOG(INFO) << "reset effective read length to " << effectiveReadLength;

    std::vector<std::string> reads(maxBucketSize, std::string(effectiveReadLength, 'X'));
    std::vector<uint32_t> offsets(maxBucketSize, 0);
    uint8_t offset{0};
    size_t subBucketSize{0};
    uint8_t bucketStringLength{0};

    while (seqFile.good() and offFile.good()) {
        uint8_t bsize{0};
        uint32_t numBytes;

        // Read the length of the bucket string (in characters)
	auto prevBucketStringLength = bucketStringLength;
        seqFile.read(reinterpret_cast<char*>(&bucketStringLength), sizeof(bucketStringLength));
	if (!seqFile.good()) { break; }

	// Read the actual bucket string
	if (bucketStringLength > 0) {
		LOG(INFO) << "bucket has new core of length " << +bucketStringLength;
        	numBytes = std::ceil(bucketStringLength / 4.0);
		std::vector<uint8_t> coreStrBytes(numBytes);
		seqFile.read(reinterpret_cast<char*>(&coreStrBytes[0]), numBytes);
		bucketString = mince::utils::twoBitDecode(&coreStrBytes.front(), bucketStringLength);
		LOG(INFO) << "new core is " << bucketString;
	} else { // 0 indicates that the current bucket has the same core string (and therefore the same length)
		LOG(INFO) << "bucket shares same core as previous bucket";
		LOG(INFO) << "repeated core is " << bucketString;
		bucketStringLength = prevBucketStringLength;
	}	

        // Read the bucket size
        seqFile.read(reinterpret_cast<char*>(&bsize), sizeof(bsize));

        subBucketSize = static_cast<uint32_t>(bsize) + 1;

	LOG(INFO) << "decoding bucket of size " << subBucketSize;
	LOG(INFO) << "bucket string length " << +bucketStringLength;
        //std::cerr << "bucket size = " << subBucketSize << "\n";

        effectiveReadLength = std::ceil((readLength - bucketStringLength) / 4.0);

        size_t subBucketCount{0};
        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            //std::cerr << "next read\n";
            //std::cerr << "effectiveReadLength = " << effectiveReadLength << "\n";
            reads[subBucketCount].resize(effectiveReadLength, 'X');
            seqFile.read(reinterpret_cast<char*>(&reads[subBucketCount][0]), effectiveReadLength);
        }

        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            offFile.read(reinterpret_cast<char*>(&offsets[subBucketCount]), sizeof(uint8_t));
        }
        // Transform offset deltas into absolute offsets
        uint8_t base = offsets.front();
        for (size_t i = 1; i < subBucketSize; ++i) {
            offsets[i] = base + offsets[i];
            base = offsets[i];
        }

        for (subBucketCount = 0; subBucketCount < subBucketSize; ++subBucketCount) {
            const char* read = reads[subBucketCount].c_str();
            uint32_t offset = offsets[subBucketCount];

                /*
                std::cerr << "HERE6\n";
                std::cerr << "bucket size = " << +subBucketSize << "\n";
                std::cerr << "bucket-string = " << bucketString << "\n";
                std::cerr << "bucket-string length = " << +bucketStringLength << "\n";
                */
            std::string s = mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(read), readLength - bucketStringLength);

                /*
                std::cerr << "[" << s << "]( " << offset << ")\n";
                std::cerr << "i = " << i << "\n";
                std::cerr << "HERE7\n";
                */

            if (offset <= s.length()) {
                //std::cerr << "HERE8\n";
              std::string recon = mince::utils::unpermute(s, bucketString, offset);

	      if (np.id == numOnsies + i) { 
		      np.apply(recon); 
		      np.getNext(nsFile);
	      }

	      if (recon.length() < readLength) {
                  std::cerr << "recon = " << recon << "\n";
                  std::cerr << "bucket string = " << bucketString << "\n";
                  std::cerr << "offset = " << +offset << "\n";
                  std::cerr << "s = " << s << "\n";
                  std::exit(1);
              }

	      bool doFlip = haveFlipFile ? flipFile.GetBit() : 0;
	      if (doFlip) { mince::utils::reverseComplement(recon); }

              ofile << ">" << numOnsies + i << "\n" << recon << std::endl;
                //std::cerr << "recon = " << recon << "\n";
            } else {
                ofile << ">" << numOnsies + i << "\nERROR RECONSTRUCTING THIS READ" << std::endl;
                std::cerr << "s.length() = " << s.length() << "\n";
                std::cerr << "offset = " << +offset << "\n";
                std::cerr << "HERE!!!!!!!!!\n";
                std::exit(1);
            }

            ++i;
            if (i % 100000 == 0) {
                std::cerr << "\r\rwrote read " << i;
            }
        }
    }
    if (seqFile.eof()) {
        std::cerr << "\nreached EOF in seq\n";
    }
    if (offFile.eof()) {
        std::cerr << "\nreached EOF in offs\n";
    }

    seqFile.close();
    offFile.close();
    if (haveFlipFile) { flipFile.Close(); }
    ofile.close();
}


