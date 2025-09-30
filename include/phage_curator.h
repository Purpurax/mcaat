#ifndef PHAGE_CURATOR_H
#define PHAGE_CURATOR_H

#include <vector>
#include <cstdint>
#include <sdbg/sdbg.h>
#include "graph_generic_func.h"
#include <string>
#include <stack>
#include <functional>
#include <map>
#include <set>
using namespace std;

struct BeamPathInfo {
    std::vector<uint64_t> path;
    double total_mult;
};

class PhageCurator {

private:
    // Variables
    SDBG& sdbg;
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> grouped_paths;
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> cycles;
    std::set<uint64_t> cycle_nodes;
    
    string _FetchNodeLastBase(size_t node);
    string _FetchFirstNode(size_t node);
public:
    //ctor
    PhageCurator(SDBG& sdbg);
    PhageCurator(SDBG& sdbg, const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles);
    std::vector<std::string> reconstructed_sequences;
    std::vector<std::vector<uint64_t>> DepthLimitedPaths(uint64_t start, int lower, int higher, std::function<void(const std::vector<uint64_t>&)> path_callback = nullptr);
    std::vector<std::vector<uint64_t>> DepthLimitedPathsAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden);
    std::vector<std::vector<uint64_t>> ExtendFromGroupedPaths(int min_depth, int max_depth);
    std::vector<std::vector<uint64_t>> ProcessPathsFromFile(const std::string& file_path, size_t min_group_size);
    void ReconstructPaths(std::vector<std::vector<uint64_t>> paths);
    void WriteSequencesToFasta(const std::string& filename);
    std::vector<BeamPathInfo> BeamSearchPaths(uint64_t start, int length, int beam_width);
};

#endif // PHAGE_CURATOR_H