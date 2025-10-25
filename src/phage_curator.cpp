#include "phage_curator.h"
#include <fstream>
#include <deque>
#include <cstdio> // for std::remove
#include <parallel_hashmap/phmap.h>  // For phmap optimizations
#include <functional>  // For std::function callback



//constructor
PhageCurator::PhageCurator(SDBG& sdbg) : sdbg(sdbg) {
    bool validation = RevalidateAllNodesButSingleton();
    if(validation){
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }
}

PhageCurator::PhageCurator(SDBG& sdbg, const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), grouped_paths(grouped_paths), cycles(cycles) {
    bool validation = RevalidateAllNodesButSingleton();
    if(validation){
        std::cout << "Graph nodes have successfully been revalidated." << std::endl;
    }
    for (const auto& [id, cycle] : cycles) {
        for (const auto& path : cycle) {
            for (uint64_t node : path) {
                cycle_nodes.insert(node);
            }
        }
        // Compute avg_spacers for this cycle
        phmap::flat_hash_set<uint64_t> unique_nodes;
        for (const auto& path : cycle) {
            for (uint64_t node : path) {
                unique_nodes.insert(node);
            }
        }
        double sum_mult = 0.0;
        for (uint64_t node : unique_nodes) {
            sum_mult += sdbg.EdgeMultiplicity(node);
        }
        avg_spacers[id] = sum_mult / unique_nodes.size();
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


bool PhageCurator::RevalidateAllNodesButSingleton() {
    
    #pragma omp parallel for
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        
        if (sdbg.EdgeMultiplicity(node) > 1 && !sdbg.IsValidEdge(node)) {
            sdbg.SetValidEdge(node);
        }
    }
    return true;
}

std::vector<std::vector<uint64_t>> PhageCurator::BeamSearchPathsAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden, int beam_width, double min_mult, double max_mult, std::function<void(const std::vector<uint64_t>&)> path_callback) {
    std::vector<std::vector<uint64_t>> all_paths;
    if (!path_callback) {
        all_paths.reserve(beam_width);  // Pre-allocate based on beam width
    }

    // Thread-local pools for paths, scores, and currents
    static thread_local std::vector<std::vector<uint64_t>> path_pool;
    static thread_local std::vector<double> score_pool;
    static thread_local std::vector<uint64_t> current_pool;

    // Clear pools
    path_pool.clear();
    score_pool.clear();
    current_pool.clear();

    // Reserve capacities based on beam width
    if (path_pool.capacity() == 0) {
        path_pool.reserve(beam_width + 64);
    }
    if (score_pool.capacity() == 0) {
        score_pool.reserve(beam_width + 64);
    }
    if (current_pool.capacity() == 0) {
        current_pool.reserve(beam_width + 64);
    }

    // Set for beam: stores {score, id}, ordered descending (begin() is highest score)
    auto comp = std::greater<std::pair<double, size_t>>();
    std::set<std::pair<double, size_t>, decltype(comp)> beam_set(comp);

    size_t unique_id = 0;

    // Initial path
    double initial_mult = sdbg.EdgeMultiplicity(start);
    if (initial_mult <= 1 || initial_mult < min_mult || initial_mult > max_mult) {
        return all_paths;
    }

    path_pool.emplace_back(std::vector<uint64_t>{start});
    path_pool.back().reserve(higher + 1);  // Pre-reserve for max depth
    score_pool.push_back(initial_mult);
    current_pool.push_back(start);

    beam_set.insert({initial_mult, unique_id++});

    while (!beam_set.empty()) {
        // Pop the best (highest score)
        auto it = beam_set.begin();
        double score = it->first;
        size_t id = it->second;
        beam_set.erase(it);

        auto& path = path_pool[id];
        uint64_t v = current_pool[id];
        int current_depth = static_cast<int>(path.size()) - 1;

        if (current_depth >= lower && current_depth <= higher) {
            if (path_callback) {
                path_callback(path);
            } else {
                all_paths.push_back(path);
            }
            continue;  // Preserve original behavior: skip expansion after collection
        }

        // Early check for zero outdegree
        if (sdbg.EdgeOutdegreeZero(v)) {
            continue;
        }

        int outdegree = sdbg.EdgeOutdegree(v);
        if (outdegree == 0 || !sdbg.IsValidEdge(v)) {
            continue;
        }

        const int MAX_EDGE_COUNT = 4;
        uint64_t neighbors[MAX_EDGE_COUNT];
        int flag = sdbg.OutgoingEdges(v, neighbors);
        if (flag == -1) {
            continue;
        }

        __builtin_prefetch(&neighbors[0], 0, 1);

        for (int i = 0; i < outdegree; ++i) {
            uint64_t neighbor = neighbors[i];

            // Cycle check in current path
            if (std::find(path.begin(), path.end(), neighbor) != path.end()) {
                continue;
            }

            // Forbidden check
            if (forbidden.find(neighbor) != forbidden.end() && neighbor != start) {
                continue;
            }

            // Multiplicity check for quality
            double mult = sdbg.EdgeMultiplicity(neighbor);
            if (mult <= 1 || mult < min_mult || mult > max_mult) {
                continue;
            }

            // Create new entry in pools
            size_t new_id = unique_id++;
            path_pool.emplace_back(path);
            path_pool.back().push_back(neighbor);
            double new_score = (score * current_depth + mult) / (current_depth + 1);
            uint64_t current = neighbor;
            current_pool.push_back(current);
            score_pool.push_back(new_score);

            // Extend along simple path, preserving original logic
            while (true) {
                uint64_t next = sdbg.NextSimplePathEdge(current);
                if (next == SDBG::kNullID || static_cast<int>(path_pool.back().size()) - 1 >= higher) {
                    break;
                }

                // Cycle check
                if (std::find(path_pool.back().begin(), path_pool.back().end(), next) != path_pool.back().end()) {
                    break;
                }

                // Forbidden check
                if (forbidden.find(next) != forbidden.end() && next != start) {
                    break;
                }

                // Multiplicity check
                double next_mult = sdbg.EdgeMultiplicity(next);
                if (next_mult <= 1 || next_mult < min_mult || next_mult > max_mult) {
                    break;
                }

                path_pool.back().push_back(next);
                int new_depth = static_cast<int>(path_pool.back().size()) - 1;
                new_score = (new_score * (new_depth - 1) + next_mult) / new_depth;
                current = next;
            }

            // Update score and current
            score_pool.back() = new_score;
            current_pool.back() = current;

            // Add to beam set
            beam_set.insert({new_score, new_id});

            // Prune the lowest score if over beam width
            if (beam_set.size() > static_cast<size_t>(beam_width)) {
                auto prune_it = beam_set.end();
                --prune_it;
                beam_set.erase(prune_it);
            }

            // Shrink if oversized
            if (path_pool.back().capacity() > static_cast<size_t>((higher + 1) * 2)) {
                path_pool.back().shrink_to_fit();
            }
        }
    }

    return all_paths;
}

std::string PhageCurator::_ReconstructPath(const std::vector<uint64_t>& path) {
    if (path.empty()) return "";
    std::string result = _FetchFirstNode(path.front());
    for (size_t i = 1; i < path.size(); ++i) {
        result += _FetchNodeLastBase(path[i]);
    }
    return result;
}

                    

// New function using beam search, mirroring FindQualityPathsDLSFromGroupedPaths
std::map<std::string,vector<string>> PhageCurator::FindQualityPathsBeamSearchFromGroupedPaths(int min_length, int max_length, const std::string& filename, int beam_width) {
    std::ofstream out(filename, std::ios::app);  // Append mode to add to file
    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }
    std::map<std::string,vector<string>> consensus_map;
    int path_count = 0;
    struct counts {
        int first;
        int second;
        int third;
    }
    counts_by_nodes = {0,0,0};
    for (const auto& [group_id, cycle_map] : grouped_paths) {
        std::vector<std::string> quality_paths;
        counts_by_nodes.first += 1;
        for (const auto& [cycle_id, paths_vec] : cycle_map) {
            counts_by_nodes.second +=1;
            string final_path;
            double avg = avg_spacers[cycle_id];
            double min_mult = 0.5 * avg;
            double max_mult = 2 * avg;
            string local_quality_path;
            
            for (const auto& path : paths_vec) {
                counts_by_nodes.third += 1;
                vector<string> local_potential_paths;
                string result_path;
                if (path.empty()) continue;
                uint64_t start = path.back();
                auto extended_paths = BeamSearchPathsAvoiding(start, min_length, max_length, cycle_nodes, beam_width, min_mult, max_mult, nullptr);
                for (const auto& ext_path : extended_paths) {
                   result_path = _FetchFirstNode(ext_path.front());
                    for (size_t i = 1; i < ext_path.size(); ++i) {
                        size_t node = ext_path[i];
                        string last_base = _FetchNodeLastBase(node);
                        result_path += last_base;
                        
                    }
                    //std::cout << "Found quality path " << path_count << " with length " << ext_path.size() << std::endl;
                    local_potential_paths.push_back(result_path);
                }
                local_quality_path = ComputeConsensusForCurrentGroup(local_potential_paths);
                if (local_quality_path.empty()) continue;
                out << ">quality_path_" << local_quality_path.substr(0, 30) << "\n" << local_quality_path << "\n";
                quality_paths.push_back(local_quality_path);
            }

        }
        std::cout<<"Processed counts:"<<counts_by_nodes.first<<","<<counts_by_nodes.second<<","<<counts_by_nodes.third<<";";
        string group_id_str = _FetchFirstNode(group_id);
        consensus_map[group_id_str] = quality_paths;
    }
    std::cout << "\nSaved in " << filename << std::endl;
    out.close();
    return consensus_map;
}

std::string PhageCurator::ComputeConsensusForCurrentGroup(vector<string> sequences) {
    // use spoa to compute consensus
    if (sequences.empty()) return "";
    // Create a spoa::Graph object with the desired alignment parameters
    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }
    auto consensus = graph.GenerateConsensus();
    return consensus;
}