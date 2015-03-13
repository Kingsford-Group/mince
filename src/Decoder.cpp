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

#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>

#include "pstream.h"

#include "bitfile.h"
#include "MinceUtils.hpp"
#include "Decoder.hpp"
#include "spdlog/spdlog.h"
#include "MinceOpts.hpp"


void dumpReadToFile(std::string& r, uint64_t readNum, LibraryFormat& fmt,
		uint32_t readLen, std::ofstream& of1, std::ofstream& of2) {
	if (fmt.type == ReadType::SINGLE_END) {
		of1 << ">" << readNum  << "\n" << r << std::endl;
	} else {
		auto singleReadLen = readLen / 2;
		auto read1 = r.substr(0, singleReadLen);
		auto read2 = r.substr(singleReadLen);
		switch (fmt.orientation) {
			case ReadOrientation::AWAY:
			case ReadOrientation::TOWARD:
				mince::utils::reverseComplement(read2);
				break;
			case ReadOrientation::SAME:
			default:
				break;
		}
		of1 << ">" << readNum << "/1\n" << read1 << "\n";
		of2 << ">" << readNum << "/2\n" << read2 << "\n";
	}

}



void Decoder::decode(MinceOpts& minceOpts, std::string& ifname, std::string& ofname) {
    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
    /*
    std::ifstream seqFile;
    std::ifstream offFile;
    std::ifstream nsFile;
    */
    bit_file_c flipFile;

    namespace bfs = boost::filesystem;

    bfs::path seqPath = ifname+".seqs.lz";
    bfs::path offPath = ifname+".offs.lz";
    bfs::path nsPath = ifname+".nlocs.lz";
    bfs::path flipPath = ifname+".flips.lz";

    auto jointLog = minceOpts.jointLog;

    if (!bfs::exists(seqPath) or !bfs::exists(offPath) or !bfs::exists(nsPath)) {
        std::string errstr = "Couldn't find mince sequence [.seqs] or offset [.offs] file or n's file [.nlocs]";
        jointLog->critical(errstr);
        throw std::logic_error(errstr);
    }

    fmt::MemoryWriter w;
    w.write("plzip -d -c -n {} {}", minceOpts.numThreads - 3, seqPath.string());
    jointLog->info("reading seq buffer with command {}", w.str());
    redi::ipstream seqFile(w.str());
    w.clear();

    w.write("plzip -d -c -n 1 {}", offPath.string());
    jointLog->info("reading offset buffer with command {}", w.str());
    redi::ipstream offFile(w.str());
    w.clear();

    w.write("plzip -d -c -n 1 {}", nsPath.string());
    jointLog->info("reading Ns buffer with command {}", w.str());
    redi::ipstream nsFile(w.str());
    w.clear();


    //seqFile.open(seqPath.string(), std::ios::in | std::ios::binary);
    //offFile.open(offPath.string(), std::ios::in | std::ios::binary);
    //nsFile.open(nsPath.string(), std::ios::in | std::ios::binary);

    bool haveFlipFile{true};
    if (!bfs::exists(flipPath)) {
        jointLog->info("no flip file found");
	    haveFlipFile = false;
    } else {
        // The "flips" file is small enough that it's first written
        // uncompressed.  We call plzip on it here.
        {
            fmt::MemoryWriter unzipCmd;
            unzipCmd.write("plzip -d -k -n 1 {}", flipPath.string());
            redi::opstream zipFlipFile(unzipCmd.str());
            if (!zipFlipFile.good()) {
                jointLog->critical("Error decompressing flip file!");
                throw std::logic_error("Error decompressing flip file!");
            }
        }

        // done compressing flip file
    	flipFile.Open(flipPath.c_str(), BF_READ);
        jointLog->info("opened flip file {}", flipPath);
    }

    uint8_t rltypeID{0};
    uint8_t rl{0};
    uint8_t kl{0};
    size_t numOnsies{0};

    seqFile.read(reinterpret_cast<char*>(&rltypeID), sizeof(rltypeID));
    seqFile.read(reinterpret_cast<char*>(&kl), sizeof(kl));
    seqFile.read(reinterpret_cast<char*>(&rl), sizeof(rl));
    seqFile.read(reinterpret_cast<char*>(&numOnsies), sizeof(numOnsies));

    LibraryFormat libFmt = LibraryFormat::formatFromID(rltypeID);

    std::ofstream ofile1;
    std::ofstream ofile2;

    std::string ofname1 = (libFmt.type == ReadType::PAIRED_END) ?
        ofname + "1.fa" : ofname + ".fa";

    ofile1.open(ofname1);
    if (libFmt.type == ReadType::PAIRED_END) {
        ofile2.open(ofname + "2.fa");
    }

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
	void getNext(redi::ipstream& ifile) {
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

    jointLog->info("bucket string length = {}", +kl);
    jointLog->info("effective read length = {}", readLength);
    jointLog->info("num onsies = ", numOnsies);

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
        dumpReadToFile(recon, j, libFmt, readLength, ofile1, ofile2);
    }

    jointLog->info("done with onsies.");
    delete onsieRead;

    uint32_t maxBucketSize{256};

    effectiveReadLength = std::ceil((readLength - kl) / 4.0);
    jointLog->info("reset effective read length to {}", effectiveReadLength);

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
        jointLog->debug("bucket has new core of length {}", +bucketStringLength);
        numBytes = std::ceil(bucketStringLength / 4.0);

        std::vector<uint8_t> coreStrBytes(numBytes);
        seqFile.read(reinterpret_cast<char*>(&coreStrBytes[0]), numBytes);
        bucketString = mince::utils::twoBitDecode(&coreStrBytes.front(), bucketStringLength);
        jointLog->debug("new core is {}", bucketString);
    } else { // 0 indicates that the current bucket has the same core string (and therefore the same length)
        jointLog->debug("bucket shares same core as previous bucket");
        jointLog->debug("repeated core is {}", bucketString);
        bucketStringLength = prevBucketStringLength;
    }

    // Read the bucket size
    seqFile.read(reinterpret_cast<char*>(&bsize), sizeof(bsize));

    subBucketSize = static_cast<uint32_t>(bsize) + 1;

    jointLog->debug("decoding bucket of size {}", subBucketSize);
    jointLog->debug("bucket string length {}", +bucketStringLength);

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

        std::string s = mince::utils::twoBitDecode(reinterpret_cast<const uint8_t*>(read), readLength - bucketStringLength);

        if (offset <= s.length()) {
            std::string recon = mince::utils::unpermute(s, bucketString, offset);

            if (np.id == numOnsies + i) {
                np.apply(recon);
                np.getNext(nsFile);
            }

            if (recon.length() < readLength) {
                fmt::MemoryWriter w;
                w.write("recon = {}", recon);
                w.write("bucket string = {}", bucketString);
                w.write("offset = {}", +offset);
                w.write("s = {}", s);
                jointLog->critical(w.str());
                throw std::logic_error(w.str());
            }

            bool doFlip = haveFlipFile ? flipFile.GetBit() : 0;
            if (doFlip) { mince::utils::reverseComplement(recon); }

            dumpReadToFile(recon, numOnsies + i, libFmt, readLength,
                    ofile1, ofile2);

        } else {
            fmt::MemoryWriter w;
            w.write("ERROR RECONSTRUCTING A READ");
            w.write("> {}", numOnsies + i);
            w.write("s.length() = {}", s.length());
            w.write("offset = {}", +offset);
            jointLog->critical(w.str());
            throw std::logic_error(w.str());
        }

        ++i;
        if (i % 100000 == 0) {
            std::cerr << "\r\rwrote read " << i;
        }
    }
    }
    if (seqFile.eof()) {
        jointLog->info("reached EOF in seq");
    }
    if (offFile.eof()) {
        jointLog->info("reached EOF in offs");
    }

    seqFile.close();
    offFile.close();
    if (haveFlipFile) {
        flipFile.Close();
        // Delete the temporarily decompressed flip file
        if (bfs::exists(ifname+".flips") and
            bfs::exists(ifname+".flips.lz")) {
            bfs::remove(ifname+".flips");
        }
    }
    ofile1.close();
    if (libFmt.type == ReadType::PAIRED_END) {
        ofile2.close();
    }


}


