#ifndef INCLUDE_CYCLE_FILTER_H_
#define INCLUDE_CYCLE_FILTER_H_

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdint>

#include "core/cft.hpp"
#include "core/Instance.hpp"
#include "algorithms/Refinement.hpp"

using namespace std;

int get_cycle_count(unordered_map<uint64_t, vector<vector<uint64_t>>>& cycle_map);

void keep_relevant_cycles(unordered_map<uint64_t, vector<vector<uint64_t>>>& cycle_map);

#endif
