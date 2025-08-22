#include "phage_curator.h"

PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {}

// Add this new function to CycleFinder:
void PhageCurator::DepthLimitedPaths(uint64_t start, int limits[]) {
    std::vector<std::vector<uint64_t>> all_paths;
    std::stack<std::pair<uint64_t, std::vector<uint64_t>>> dls_stack;
    dls_stack.emplace(start, std::vector<uint64_t>{start});

    while (!dls_stack.empty()) {
        auto [v, path] = dls_stack.top();
        dls_stack.pop();

        if ((int)path.size() - 1 >= limits[0] && (int)path.size() - 1 <= limits[1]) {
            all_paths.push_back(path);
            continue;
        }

        std::unordered_set<uint64_t> adj;
        adj.reserve(kAlphabetSize + 1);
        this->_GetOutgoings(v, adj);

        for (uint64_t neighbor : adj) {
            // Prevent cycles in path
            if (std::find(path.begin(), path.end(), neighbor) == path.end()) {
                auto new_path = path;
                new_path.push_back(neighbor);
                dls_stack.emplace(neighbor, std::move(new_path));
            }
        }
    }
    this->results = std::move(all_paths);
    return;
}

