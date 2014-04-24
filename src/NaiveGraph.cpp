#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "jellyfish/jellyfish.hpp"
#include "jellyfish/mer_dna.hpp"

#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>

template <typename VertDesc, typename Graph>
void collapsePath(VertDesc startVert, Graph& g, size_t& numNodes, size_t& numEdges) {


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> vertsToRemove;
    auto neighbor = *(boost::adjacent_vertices(startVert, g).first);
    size_t pathLength{1};
    while (boost::out_degree(neighbor, g) == 1) {
        vertsToRemove.push_back(neighbor);
        neighbor = *(boost::adjacent_vertices(neighbor, g).first);
        ++pathLength;
    }
    end = std::chrono::system_clock::now();

     auto outNeigh = boost::adjacent_vertices(neighbor, g);
    for (auto it = outNeigh.first; it != outNeigh.second; ++it) {
        boost::add_edge(startVert, *it, g);
    }

    for (auto v : vertsToRemove) {
        boost::clear_vertex(v, g);
    }

    std::chrono::duration<double> elapsed_seconds = end-start;

    if (pathLength > 0) {
        numNodes -= pathLength;
        numEdges -= pathLength;

        std::cerr << "\r\rcollapsed path of length " << pathLength
            << " V[" << numNodes << "], E[" << numEdges << "]"
            << " : (" << elapsed_seconds.count() << "s)";
    }
}


int main(int argc, char* argv[]) {

    namespace bfs = boost::filesystem;

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;

    bfs::path thashFile = std::string(argv[1]);

    // Read in the Jellyfish hash of the transcripts
    std::ifstream jellyfishDB(thashFile.c_str());
    if (!jellyfishDB.good()) {
        std::cerr << "Couldn't open the Jellyfish hash [" << thashFile << "] quitting\n";
        std::exit(-1);
    }
    jellyfish::file_header header;
    header.read(jellyfishDB);

    std::cerr << "transcript hash size is " << header.size() << "\n";
    size_t nkeys = header.size();
    // Since JF2, key_len() is in terms of bits so merLen = keyLen / 2;
    size_t merLen = header.key_len() / 2;
    jellyfish::mer_dna::k(merLen);

    size_t tot{0};
    size_t distinct{0};

    const jellyfish::RectangularBinaryMatrix m(header.matrix());
    size_t sizeMask{header.size() - 1};

    auto merHasher = [&](const jellyfish::mer_dna& mer) -> size_t {
            //auto key = mer.get_bits(0, jellyfish::mer_dna::k());
            return m.times(mer) & sizeMask;
        };

    std::unordered_map<jellyfish::mer_dna, size_t, decltype(merHasher)> merIDMap(10, merHasher);

    if (!header.format().compare(binary_dumper::format)) {
        binary_reader reader(jellyfishDB, &header);
        while ( reader.next() ) {
            tot += reader.val();
            merIDMap[reader.key()] = distinct;
            ++distinct;
            if (distinct % 1000000 == 0) {
                std::cerr << "\r\rinserted " << distinct << " keys into hash";
            }
        }
    }
    std::cerr << "\n";
    std::cerr << "There were " << tot << " total k-mers and " << distinct << " distinct k-mers\n";
    Graph g(distinct);

    size_t numEdge{0};
    size_t nv{0};
    for (auto& kv : merIDMap) {
        ++nv;
        auto key = kv.first;

        size_t numNeighbor{numEdge};

        jellyfish::mer_dna extA(key);
        extA.shift_left('A');
        auto it = merIDMap.find(extA);
        if (it != merIDMap.end()) {
            boost::add_edge(kv.second, it->second, g);
            ++numEdge;
        }

        jellyfish::mer_dna extC(key);
        extC.shift_left('C');
        it = merIDMap.find(extC);
        if (it != merIDMap.end()) {
            boost::add_edge(kv.second, it->second, g);
            ++numEdge;
        }

        jellyfish::mer_dna extG(key);
        extG.shift_left('G');
        it = merIDMap.find(extG);
        if (it != merIDMap.end()) {
            boost::add_edge(kv.second, it->second, g);
            ++numEdge;
        }

        jellyfish::mer_dna extT(key);
        extT.shift_left('T');
        it = merIDMap.find(extT);
        if (it != merIDMap.end()) {
            boost::add_edge(kv.second, it->second, g);
            ++numEdge;
        }

        if (nv % 1000000 == 0) {
            std::cerr << "\r\rprocessed " << nv << " graph verts; have " << numEdge << " edges\n";
        }
    }

    size_t numProc{0};
    boost::graph_traits<Graph>::vertex_iterator vi, vi_end, next;
    std::tie(vi, vi_end) = vertices(g);
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pathStarts;
    for (next = vi; vi != vi_end; vi = next) {
        ++next;
        if (boost::in_degree(*vi, g) != 1 and boost::out_degree(*vi, g) == 1) {
            pathStarts.push_back(*vi);
        }
        ++numProc;
        if (numProc % 1000000 == 0) {
            std::cerr << "\r\rprocessed " << numProc << " vertices, have " << pathStarts.size() << " paths";
        }
    }


    size_t numNodes = boost::num_vertices(g);
    size_t numEdges = boost::num_edges(g);
    for (auto v : pathStarts) {
        collapsePath(v, g, numNodes, numEdges);
    }

    std::cerr << "\n\nfinal # of nodes = " << boost::num_vertices(g) << ", final # edges = "
              << boost::num_edges(g) << "\n";
}
