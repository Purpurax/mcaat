#include "cycle_filter.h"

int get_cycle_count(unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map) {
	int counter = 0;
    for(const auto& [_, cycles] : cycles_map) {
        counter += cycles.size();
    }
    return counter;
}

void keep_relevant_cycles(unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map) {
    if (get_cycle_count(cycles_map) < 3) {
        return;
    }
    
    auto set_cover_instance = cft::Instance();
    vector<vector<uint64_t>> sets;
    unordered_set<uint64_t> universe;
    for (const auto& [_, cycles] : cycles_map) {
        for (const auto& cycle : cycles) {
            if (!cycle.empty()) {
                sets.push_back(cycle);
                for (const auto& node : cycle) {
                    universe.insert(node);
                }
            }
        }
    }

    // for (const auto& element : universe) {
    //     vector<int32_t> new_row;
    //     for (size_t idx = 0; idx < sets.size(); ++idx) {
    //         if (std::find(sets[idx].begin(), sets[idx].end(), element) != sets[idx].end()) {
    //         new_row.push_back(idx);
    //         }
    //     }
    //     set_cover_instance.rows.push_back(new_row);
    // }
    
    for (const auto& set : sets) {
        set_cover_instance.cols.push_back(set);
        set_cover_instance.costs.push_back(1.0);
    }

    cft::fill_rows_from_cols(set_cover_instance.cols, universe.size(), set_cover_instance.rows);
    CFT_IF_DEBUG(cft::col_and_rows_check(set_cover_instance.cols, set_cover_instance.rows));
    
    auto env = cft::Environment();
    env.time_limit = 10.0;
    env.verbose = 5;
    env.timer.restart();

    // Use cft to solve minimal set cover
    if (!sets.empty() && !universe.empty()) {
        auto solution = cft::run(env, set_cover_instance);
    }

    // unordered_set<size_t> chosen_idxs(solution.sol.idxs.begin(), solution.sol.idxs.end());

    // // Remove cycles not in the solution
    // unordered_map<uint64_t, vector<vector<uint64_t>>> filtered_map;
    // for (size_t idx : solution.sol.idxs) {
    //     const auto& cycle = sets[idx];
    //     if (!cycle.empty()) {
    //         uint64_t start_node = cycle.front();
    //         filtered_map[start_node].push_back(cycle);
    //     }
    // }
    // cycles_map = std::move(filtered_map);
}
