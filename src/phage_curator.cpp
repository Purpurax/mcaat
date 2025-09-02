#include "phage_curator.h"
#include <fstream>
#include <deque>
#include <cstdio> // for std::remove



//constructor
PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {
 
}

string PhageCurator::_FetchFirstNode(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

string PhageCurator::_FetchNodeLastBase(size_t node) {
    uint8_t seq[sdbg.k()];
    sdbg.GetLabel(node, seq);
    char base = "ACGT"[seq[0] - 1];
    return string(1, base);
}

void PhageCurator::ReconstructPaths(std::vector<std::vector<uint64_t>> paths) {
    
    for (const auto& path : paths) {
        string result_path = _FetchFirstNode(path.front());
        // start loop from second element
        for (size_t i = 1; i < path.size(); ++i) {
            size_t node = path[i];
            string last_base = _FetchNodeLastBase(node);
            result_path += last_base;
        }
        reconstructed_sequences.push_back(result_path);
    }
    return;
}

// Add this new function to CycleFinder:
std::vector<std::vector<uint64_t>> PhageCurator::DepthLimitedPaths(uint64_t start, int lower, int higher) {
    std::vector<std::vector<uint64_t>> all_paths;
    std::deque<std::pair<uint64_t, std::vector<uint64_t>>> dls_stack;
    dls_stack.emplace_back(start, std::vector<uint64_t>{start});

    // Spill-to-disk parameters
    const size_t max_stack_nodes = 10000; // keep last 10000 in memory
    const std::string temp_file = "dls_stack_spill.tmp";

    auto spill_stack_to_file = [&](std::deque<std::pair<uint64_t, std::vector<uint64_t>>>& stack) {
        std::ofstream ofs(temp_file, std::ios::app | std::ios::binary);
        // Spill oldest entries, keep last max_stack_nodes
        while (stack.size() > max_stack_nodes) {
            auto& entry = stack.front();
            uint64_t node = entry.first;
            const auto& path = entry.second;
            ofs.write(reinterpret_cast<const char*>(&node), sizeof(node));
            size_t path_size = path.size();
            ofs.write(reinterpret_cast<const char*>(&path_size), sizeof(path_size));
            ofs.write(reinterpret_cast<const char*>(path.data()), path_size * sizeof(uint64_t));
            stack.pop_front();
        }
        ofs.close();
    };

    auto load_stack_from_file = [&](std::deque<std::pair<uint64_t, std::vector<uint64_t>>>& stack) {
        std::ifstream ifs(temp_file, std::ios::binary);
        if (!ifs) return;
        while (ifs.peek() != EOF) {
            uint64_t node;
            size_t path_size;
            ifs.read(reinterpret_cast<char*>(&node), sizeof(node));
            ifs.read(reinterpret_cast<char*>(&path_size), sizeof(path_size));
            std::vector<uint64_t> path(path_size);
            ifs.read(reinterpret_cast<char*>(path.data()), path_size * sizeof(uint64_t));
            stack.emplace_back(node, std::move(path));
        }
        ifs.close();
        std::remove(temp_file.c_str());
    };

    while (!dls_stack.empty() || std::ifstream(temp_file).good()) {
        // Spill stack if too large (keep only last max_stack_nodes)
        if (dls_stack.size() > max_stack_nodes) {
            std::cout << "Triggering spill, stack size: " << dls_stack.size() << std::endl;
            std::cout << "Spilling " << dls_stack.size() - max_stack_nodes << " entries to " << temp_file << std::endl;
            spill_stack_to_file(dls_stack);
        }
        // Reload stack from file only if empty
        if (dls_stack.empty() && std::ifstream(temp_file).good()) {
            load_stack_from_file(dls_stack);
        }
        if (dls_stack.empty()) break;

        auto [v, path] = dls_stack.back();
        dls_stack.pop_back();

        if ((int)path.size() - 1 >= lower && (int)path.size() - 1 <= higher) {
            all_paths.push_back(path);
            continue;
        }

        std::unordered_set<uint64_t> adj;
        adj.reserve(kAlphabetSize + 1);
        graph_generic_func::_GetOutgoings(v, adj, sdbg);

        for (uint64_t neighbor : adj) {
            // Prevent cycles in path
            if (std::find(path.begin(), path.end(), neighbor) == path.end()) {
                auto new_path = path;
                new_path.push_back(neighbor);
                dls_stack.emplace_back(neighbor, std::move(new_path));
            }
        }
    }
    std::remove(temp_file.c_str());
    return all_paths;
}


std::vector<BeamPathInfo> PhageCurator::BeamSearchPaths(uint64_t start, int length, int beam_width) {
    std::vector<BeamPathInfo> beam = {{ {start}, 0.0 }};
    for (int step = 0; step < length; ++step) {
        std::vector<BeamPathInfo> candidates;
        for (const auto& info : beam) {
            uint64_t current = info.path.back();
            std::unordered_set<uint64_t> adj;
            graph_generic_func::_GetOutgoings(current, adj, sdbg);
            for (uint64_t neighbor : adj) {
                // Prevent cycles
                if (std::find(info.path.begin(), info.path.end(), neighbor) == info.path.end()) {
                    double mult = sdbg.EdgeMultiplicity(neighbor);
                    auto new_path = info.path;
                    new_path.push_back(neighbor);
                    candidates.push_back({new_path, info.total_mult + mult});
                }
            }
        }
        // Keep only top beam_width candidates (max multiplicity)
        std::sort(candidates.begin(), candidates.end(),
                  [](const BeamPathInfo& a, const BeamPathInfo& b) { return a.total_mult > b.total_mult; });
        if ((int)candidates.size() > beam_width) candidates.resize(beam_width);
        beam = std::move(candidates);
        if (beam.empty()) break;
    }
    return beam;
}