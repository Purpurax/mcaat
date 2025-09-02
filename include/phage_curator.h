#ifndef PHAGE_CURATOR_H
#define PHAGE_CURATOR_H

#include <vector>
#include <cstdint>
#include <sdbg/sdbg.h>
#include "graph_generic_func.h"
#include <string>
#include <stack>
using namespace std;

struct BeamPathInfo {
    std::vector<uint64_t> path;
    double total_mult;
};

class PhageCurator {

private:
    // Variables
    SDBG& sdbg;
    
    string _FetchNodeLastBase(size_t node);
    string _FetchFirstNode(size_t node);
public:
    //ctor
    PhageCurator(SDBG& sdbg);
    std::vector<std::string> reconstructed_sequences;
    std::vector<std::vector<uint64_t>> DepthLimitedPaths(uint64_t start, int lower,int higher);
    void ReconstructPaths(std::vector<std::vector<uint64_t>> paths);
    std::vector<BeamPathInfo> BeamSearchPaths(uint64_t start, int length, int beam_width);
};

#endif // PHAGE_CURATOR_H