#include "cycle_filter.h"

int get_cycle_count(std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles_map) {
	int counter = 0;
    for(const auto& [_, cycles] : cycles_map) {
        counter += cycles.size();
    }
    return counter;
}

void keep_relevant_cycles(std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles_map) {
    auto set_cover_instance = cft::Instance();
    std::vector<std::vector<uint64_t>> sets;
    for (const auto& [_, cycles] : cycles_map) {
        for (const auto& cycle : cycles) {
            sets.push_back(cycle);
        }
    }

    for (const auto& set : sets) {
        set_cover_instance.cols.push_back(set);
        set_cover_instance.costs.push_back(1);
    }
    auto env = cft::Environment();

    // Use cft to solve minimal set cover
    auto solution = cft::run(env, set_cover_instance);

    // Remove cycles not in the solution
    size_t idx = 0;
    for (auto& [_, cycles] : cycles_map) {
        std::vector<std::vector<uint64_t>> filtered;
        for (const auto& cycle : cycles) {
            if (solution.sol.idxs[idx]) {
                filtered.push_back(cycle);
            }
            ++idx;
        }
        cycles = std::move(filtered);
    }
}
