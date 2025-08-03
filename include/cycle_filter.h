#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include "core/cft.hpp"
#include "core/Instance.hpp"
#include "algorithms/Refinement.hpp"

int get_cycle_count(std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycle_map);
void keep_relevant_cycles(std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycle_map);
