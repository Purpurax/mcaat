#ifndef PHAGE_CURATOR_H
#define PHAGE_CURATOR_H

#include <vector>
#include <cstdint>
#include <sdbg/sdbg.h>
#include "graph_generic_func.h"

class PhageCurator {

private:
    // Variables
    SDBG& sdbg;
public:
    //ctor
    PhageCurator(SDBG& sdbg);
    std::vector<std::vector<uint64_t>> results;
    void DepthLimitedPaths(uint64_t start, int limits[]);

};

#endif // PHAGE_CURATOR_H