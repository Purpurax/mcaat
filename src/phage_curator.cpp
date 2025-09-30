#include "phage_curator.h"
#include <fstream>
#include <deque>
#include <cstdio> // for std::remove
#include <parallel_hashmap/phmap.h>  // For phmap optimizations
#include <functional>  // For std::function callback



//constructor
PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {
 
}

PhageCurator::PhageCurator(SDBG& sdbg, const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), grouped_paths(grouped_paths), cycles(cycles) {
    for (const auto& [id, cycle] : cycles) {
        for (const auto& path : cycle) {
            for (uint64_t node : path) {
                cycle_nodes.insert(node);
            }
        }
    }
}

string PhageCurator::_FetchFirstNode(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    sdbg.GetLabel(node, seq);
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

void PhageCurator::WriteSequencesToFasta(const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    int count = 1;
    for (const auto& seq : reconstructed_sequences) {
        out << ">sequence_" << count << "\n" << seq << "\n";
        count++;
    }
}

// Adapted DepthLimitedPaths with optimizations from DepthLevelSearch (no logic changes)
std::vector<std::vector<uint64_t>> PhageCurator::DepthLimitedPaths(uint64_t start, int lower, int higher, std::function<void(const std::vector<uint64_t>&)> path_callback) {
    std::vector<std::vector<uint64_t>> all_paths;
    if (!path_callback) {
        all_paths.reserve(1024);  // Pre-allocate if collecting
    }
    std::vector<std::pair<uint64_t, std::vector<uint64_t>>> dls_stack;
    dls_stack.emplace_back(start, std::vector<uint64_t>{start});

    // Thread-local pools (extended for paths)
    static thread_local std::vector<std::pair<uint64_t, std::vector<uint64_t>>> stack_pool;
    static thread_local phmap::flat_hash_set<uint64_t> visited_pool;
    static thread_local std::vector<std::vector<uint64_t>> path_pool;

    // Clear but keep capacity (reuse pattern)
    stack_pool.clear();
    visited_pool.clear();
    path_pool.clear();

    // Reserve initial capacities
    if (stack_pool.capacity() == 0) {
        stack_pool.reserve(64);
    }
    if (visited_pool.capacity() == 0) {
        visited_pool.reserve(1024);
    }
    if (path_pool.capacity() == 0) {
        path_pool.reserve(64);
    }

    // Use pooled structures
    auto& pooled_visited = visited_pool;
    auto& pooled_paths = path_pool;

    // Cache values
    const uint64_t start_node = start;
    const int min_depth = lower;
    const int max_depth = higher;

    // Helper to get or create a path vector from pool (with shrinking)
    auto get_path_from_pool = [&]() -> std::vector<uint64_t>& {
        if (pooled_paths.empty()) {
            pooled_paths.emplace_back();
            pooled_paths.back().reserve(max_depth + 1);  // Pre-reserve based on max depth
        }
        auto& path = pooled_paths.back();
        path.clear();
        return path;
    };

    while (!dls_stack.empty()) {
        auto [v, path] = dls_stack.back();
        dls_stack.pop_back();

        int current_depth = (int)path.size() - 1;
        if (__builtin_expect(current_depth >= min_depth && current_depth <= max_depth, 1)) {
            if (path_callback) {
                path_callback(path);
            } else {
                all_paths.push_back(path);
            }
            continue;
        }

        // Early check for zero outdegree (from DepthLevelSearch)
        if (__builtin_expect(sdbg.EdgeOutdegreeZero(v), 0)) {
            continue;
        }

        int outdegree = sdbg.EdgeOutdegree(v);
        if (__builtin_expect(outdegree == 0 || !sdbg.IsValidEdge(v), 0)) {
            continue;
        }

        // Fixed-size array for neighbors (from DepthLevelSearch)
        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT];
        int flag = sdbg.OutgoingEdges(v, neighbors);
        if (__builtin_expect(flag == -1, 0)) {
            continue;
        }

        // Prefetch neighbors (from DepthLevelSearch)
        __builtin_prefetch(&neighbors[0], 0, 1);

        for (int i = 0; i < outdegree; ++i) {
            uint64_t neighbor = neighbors[i];
            // Prevent cycles (logic unchanged)
            if (__builtin_expect(std::find(path.begin(), path.end(), neighbor) == path.end(), 1)) {
                // Global visited check (phmap from DepthLevelSearch)
                auto visited_it = pooled_visited.find(neighbor);
                bool not_visited = (visited_it == pooled_visited.end());
                bool is_start_revisit = (neighbor == start_node && current_depth > 0);

                if (__builtin_expect(not_visited || is_start_revisit, 1)) {
                    pooled_visited.insert(neighbor);
                    auto new_path = get_path_from_pool();
                    new_path = path;
                    new_path.push_back(neighbor);
                    // Extend along simple path using Megahit's NextSimplePathEdge
                    uint64_t current = neighbor;
                    while (true) {
                        uint64_t next = sdbg.NextSimplePathEdge(current);
                        if (next == SDBG::kNullID || (int)new_path.size() - 1 >= max_depth) {
                            break;
                        }
                        new_path.push_back(next);
                        current = next;
                    }
                    dls_stack.emplace_back(current, std::move(new_path));
                    // After moving, shrink the pooled vector if it's oversized
                    if (!pooled_paths.empty()) {
                        auto& last_path = pooled_paths.back();
                        if (last_path.capacity() > static_cast<size_t>((max_depth + 1) * 2)) {  // If capacity is >2x needed
                            last_path.shrink_to_fit();  // Shrink to actual size
                        }
                    }
                }
            }
        }
    }
    return all_paths;
}


std::vector<std::vector<uint64_t>> PhageCurator::DepthLimitedPathsAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden) {
    std::vector<std::vector<uint64_t>> all_paths;
    all_paths.reserve(1024);  // Pre-allocate
    std::vector<std::pair<uint64_t, std::vector<uint64_t>>> dls_stack;
    dls_stack.emplace_back(start, std::vector<uint64_t>{start});

    // Thread-local pools
    static thread_local std::vector<std::pair<uint64_t, std::vector<uint64_t>>> stack_pool;
    static thread_local phmap::flat_hash_set<uint64_t> visited_pool;
    static thread_local std::vector<std::vector<uint64_t>> path_pool;

    stack_pool.clear();
    visited_pool.clear();
    path_pool.clear();

    if (stack_pool.capacity() == 0) stack_pool.reserve(64);
    if (visited_pool.capacity() == 0) visited_pool.reserve(1024);
    if (path_pool.capacity() == 0) path_pool.reserve(64);

    auto& pooled_visited = visited_pool;
    auto& pooled_paths = path_pool;

    const uint64_t start_node = start;
    const int min_depth = lower;
    const int max_depth = higher;

    auto get_path_from_pool = [&]() -> std::vector<uint64_t>& {
        if (pooled_paths.empty()) {
            pooled_paths.emplace_back();
            pooled_paths.back().reserve(max_depth + 1);
        }
        auto& path = pooled_paths.back();
        path.clear();
        return path;
    };

    while (!dls_stack.empty()) {
        auto [v, path] = dls_stack.back();
        dls_stack.pop_back();

        int current_depth = (int)path.size() - 1;
        if (__builtin_expect(current_depth >= min_depth && current_depth <= max_depth, 1)) {
            all_paths.push_back(path);
            continue;
        }

        if (__builtin_expect(sdbg.EdgeOutdegreeZero(v), 0)) continue;

        int outdegree = sdbg.EdgeOutdegree(v);
        if (__builtin_expect(outdegree == 0 || !sdbg.IsValidEdge(v), 0)) continue;

        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT];
        int flag = sdbg.OutgoingEdges(v, neighbors);
        if (__builtin_expect(flag == -1, 0)) continue;

        __builtin_prefetch(&neighbors[0], 0, 1);

        for (int i = 0; i < outdegree; ++i) {
            uint64_t neighbor = neighbors[i];
            if (__builtin_expect(std::find(path.begin(), path.end(), neighbor) == path.end(), 1)) {
                // Check forbidden, but allow start
                if (forbidden.find(neighbor) != forbidden.end() && neighbor != start_node) continue;
                auto visited_it = pooled_visited.find(neighbor);
                bool not_visited = (visited_it == pooled_visited.end());
                bool is_start_revisit = (neighbor == start_node && current_depth > 0);

                if (__builtin_expect(not_visited || is_start_revisit, 1)) {
                    pooled_visited.insert(neighbor);
                    auto new_path = get_path_from_pool();
                    new_path = path;
                    new_path.push_back(neighbor);
                    uint64_t current = neighbor;
                    while (true) {
                        uint64_t next = sdbg.NextSimplePathEdge(current);
                        if (next == SDBG::kNullID || (int)new_path.size() - 1 >= max_depth) break;
                        // Also check forbidden for extended
                        if (forbidden.find(next) != forbidden.end() && next != start_node) break;
                        new_path.push_back(next);
                        current = next;
                    }
                    dls_stack.emplace_back(current, std::move(new_path));
                    if (!pooled_paths.empty()) {
                        auto& last_path = pooled_paths.back();
                        if (last_path.capacity() > static_cast<size_t>((max_depth + 1) * 2)) {
                            last_path.shrink_to_fit();
                        }
                    }
                }
            }
        }
    }
    return all_paths;
}


std::vector<std::vector<uint64_t>> PhageCurator::ExtendFromGroupedPaths(int min_depth, int max_depth) {
    std::vector<std::vector<uint64_t>> all_extended_paths;
    for (const auto& [outer_key, inner_map] : grouped_paths) {
        for (const auto& [inner_key, paths_vec] : inner_map) {
            for (const auto& path : paths_vec) {
                if (!path.empty()) {
                    uint64_t start = path.back();
                    auto paths = DepthLimitedPathsAvoiding(start, min_depth, max_depth, cycle_nodes);
                    all_extended_paths.insert(all_extended_paths.end(), paths.begin(), paths.end());
                }
            }
        }
    }
    return all_extended_paths;
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