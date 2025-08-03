#ifndef JUMPS_H
#define JUMPS_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <string>
#include <filesystem>
#include "sdbg/sdbg.h"
#include "settings.h"
#include "sequence/io/sequence_lib.h"
#include "idba/sequence.h"

struct Jump {
    uint64_t start_k_mer_id;
    uint64_t end_k_mer_id;
    uint64_t nodes_in_between;
};

std::vector<Jump> get_jumps_from_reads(SDBG& sdbg, Settings& settings);
#endif
