#include "spacer_ordering.h"

class TarjanSCC {
private:
    SDBG& sdbg;
    std::unordered_map<uint64_t, int> index_map;
    std::unordered_map<uint64_t, int> lowlink_map;
    std::unordered_set<uint64_t> on_stack;
    std::stack<uint64_t> stack;
    std::vector<std::vector<uint64_t>> components;
    int index_counter;

    void strongconnect(uint64_t node) {
        index_map[node] = index_counter;
        lowlink_map[node] = index_counter;
        index_counter++;
        stack.push(node);
        on_stack.insert(node);

        // Get outgoing edges
        int outdegree = sdbg.EdgeOutdegree(node);
        if (outdegree > 0) {
            uint64_t* outgoings = new uint64_t[outdegree];
            if (sdbg.OutgoingEdges(node, outgoings) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    uint64_t neighbor = outgoings[i];
                    if (!sdbg.IsValidEdge(neighbor)) continue;
                    
                    if (index_map.find(neighbor) == index_map.end()) {
                        strongconnect(neighbor);
                        lowlink_map[node] = std::min(lowlink_map[node], lowlink_map[neighbor]);
                    } else if (on_stack.find(neighbor) != on_stack.end()) {
                        lowlink_map[node] = std::min(lowlink_map[node], index_map[neighbor]);
                    }
                }
            }
            delete[] outgoings;
        }

        // If node is a root node, pop the stack and create an SCC
        if (lowlink_map[node] == index_map[node]) {
            std::vector<uint64_t> component;
            uint64_t w;
            do {
                w = stack.top();
                stack.pop();
                on_stack.erase(w);
                component.push_back(w);
            } while (w != node);
            
            if (component.size() > 1) { // Only keep non-trivial components
                components.push_back(component);
            }
        }
    }

public:
    TarjanSCC(SDBG& graph) : sdbg(graph), index_counter(0) {}

    std::vector<std::vector<uint64_t>> find_components() {
        // Find all valid nodes first
        std::vector<uint64_t> valid_nodes;
        for (uint64_t node = 0; node < sdbg.size(); ++node) {
            if (sdbg.IsValidEdge(node)) {
                valid_nodes.push_back(node);
            }
        }

        // Run Tarjan's algorithm on all unvisited valid nodes
        for (uint64_t node : valid_nodes) {
            if (index_map.find(node) == index_map.end()) {
                strongconnect(node);
            }
        }

        return components;
    }
};

void keep_crispr_regions(
    SDBG& sdbg,
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles
) {
    std::unordered_set<uint64_t> cycle_nodes;
    for (const auto& [start_node, cycle_list] : cycles) {
        // Add the start node
        cycle_nodes.insert(start_node);
        
        // Add all nodes from all cycles
        for (const auto& cycle : cycle_list) {
            for (uint64_t node : cycle) {
                cycle_nodes.insert(node);
            }
        }
    }

    #pragma omp parallel for
    for (uint64_t node_id = 0; node_id < sdbg.size(); ++node_id) {
        if (sdbg.IsValidEdge(node_id)) {
            // If this node is not in any cycle, invalidate it
            if (cycle_nodes.find(node_id) == cycle_nodes.end()) {
                sdbg.SetInvalidEdge(node_id);
            }
        }
    }

    sdbg.FreeMultiplicity();
}



std::vector<Graph> divide_graph_into_subgraphs(SDBG& sdbg) {
    std::vector<Graph> subgraphs;
    
    // Find strongly connected components using Tarjan's algorithm
    TarjanSCC tarjan(sdbg);
    auto components = tarjan.find_components();
    
    std::cout << "Found " << components.size() << " strongly connected components" << std::endl;
    
    // Convert each component to a SubGraph
    for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
        const auto& component = components[comp_idx];
        Graph subgraph;
        
        // Add all edges within this component
        for (uint64_t node : component) {
            if (!sdbg.IsValidEdge(node)) continue;
            
            int outdegree = sdbg.EdgeOutdegree(node);
            if (outdegree > 0) {
                uint64_t* outgoings = new uint64_t[outdegree];
                if (sdbg.OutgoingEdges(node, outgoings) != -1) {
                    for (int i = 0; i < outdegree; ++i) {
                        uint64_t neighbor = outgoings[i];
                        // Only add edge if neighbor is in the same component
                        if (std::find(component.begin(), component.end(), neighbor) != component.end()) {
                            subgraph.add_edge(node, neighbor);
                        }
                    }
                }
                delete[] outgoings;
            }
        }
        
        if (subgraph.nodes.size() > 0) {
            subgraphs.push_back(std::move(subgraph));
        }
    }

    return subgraphs;
}

std::vector<Graph> get_crispr_regions(
    SDBG& sdbg,
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles
) {
    keep_crispr_regions(sdbg, cycles);
    return divide_graph_into_subgraphs(sdbg);
}

std::vector<Jump> get_relevant_jumps(const Graph& graph, std::vector<Jump>& all_jumps) {
    std::vector<Jump> relevant_jumps;

    for (const auto& jump : all_jumps) {
        if (graph.nodes.find(jump.start_k_mer_id) != graph.nodes.end() &&
            graph.nodes.find(jump.end_k_mer_id) != graph.nodes.end()) {
            relevant_jumps.push_back(jump);
        }
    }

    return relevant_jumps;
}

std::vector<std::vector<uint64_t>> get_relevant_cycles(
    const Graph& graph,
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& all_cycles_map
) {
    std::vector<std::vector<uint64_t>> relevant_cycles;

    for (const auto& [_, cycles] : all_cycles_map) {
        for (const auto& cycle : cycles) {
            bool all_nodes_in_graph = true;
            for (uint64_t node : cycle) {
                if (graph.nodes.find(node) == graph.nodes.end()) {
                    all_nodes_in_graph = false;
                    break;
                }
            }
            if (all_nodes_in_graph) {
                relevant_cycles.push_back(cycle);
            }
        }
    }

    return relevant_cycles;
}

std::unordered_map<uint64_t, int32_t> get_node_to_cycle_map(std::vector<std::vector<uint64_t>>& cycles) {
    std::vector<std::unordered_set<uint64_t>> cycle_node_sets;
    for (const auto& cycle : cycles) {
        cycle_node_sets.emplace_back(cycle.begin(), cycle.end());
    }

    std::unordered_map<uint64_t, int32_t> node_to_cycle_map;
    for (size_t i = 0; i < cycle_node_sets.size(); ++i) {
        // Build union of all other sets
        std::unordered_set<uint64_t> other_nodes;
        for (size_t j = 0; j < cycle_node_sets.size(); ++j) {
            if (j == i) continue;
            other_nodes.insert(cycle_node_sets[j].begin(), cycle_node_sets[j].end());
        }
        // Find unique nodes in cycle i
        for (uint64_t node : cycle_node_sets[i]) {
            if (other_nodes.find(node) == other_nodes.end()) {
                node_to_cycle_map[node] = static_cast<int32_t>(i);
            }
        }
    }
    return node_to_cycle_map;
}

std::vector<int32_t> get_all_cycle_indices(
    std::unordered_map<uint64_t, int32_t> node_to_cycle_map
) {
    std::vector<int32_t> cycle_indices;
    for (const auto& [node, cycle_idx] : node_to_cycle_map) {
        if (std::find(cycle_indices.begin(), cycle_indices.end(), cycle_idx) == cycle_indices.end()) {
            cycle_indices.push_back(cycle_idx);
        }
    }
    return cycle_indices;
}

std::vector<std::vector<uint64_t>> find_all_possible_paths(
    const Graph& graph,
    uint64_t start,
    uint64_t end,
    uint64_t nodes_left
) {
    std::vector<std::vector<uint64_t>> paths;

    if (nodes_left == 0 && start == end) {
        paths.push_back({end});
        return paths;
    } else if (nodes_left == 0) {
        return paths;
    }

    if (graph.adjacency_list.find(start) == graph.adjacency_list.end()) {
        return paths;
    }

    const auto& neighbors = graph.adjacency_list.find(start)->second;
    for (uint64_t next_node : neighbors) {
        auto sub_paths = find_all_possible_paths(graph, next_node, end, nodes_left - 1);
        for (auto& sub_path : sub_paths) {
            sub_path.insert(sub_path.begin(), start);
            paths.push_back(sub_path);
        }
    }

    return paths;
}

std::vector<std::tuple<uint32_t, uint32_t>> generate_constraints(
    const Graph& graph,
    std::vector<Jump>& jumps,
    std::unordered_map<uint64_t, int32_t>& node_to_cycle_map
) {
    std::vector<std::tuple<uint32_t, uint32_t>> constraints;

    for (const auto& jump : jumps) {
        std::vector<std::vector<uint64_t>> all_paths = find_all_possible_paths(
            graph,
            jump.start_k_mer_id,
            jump.end_k_mer_id,
            jump.nodes_in_between + 1
        );

        if (all_paths.empty()) {
            continue;
        }

        std::vector<std::unordered_set<uint64_t>> path_with_cycle_indices;
        for (int i = 0; i < all_paths[0].size(); ++i) {
            std::unordered_set<uint64_t> cycle_indices;
            for (const auto& path : all_paths) {
                if (i < path.size()) {
                    uint64_t node = path[i];
                    auto it = node_to_cycle_map.find(node);
                    if (it != node_to_cycle_map.end()) {
                        cycle_indices.insert(it->second);
                    }
                }
            }
            path_with_cycle_indices.push_back(cycle_indices);
        }

        std::vector<uint64_t> path_with_unique_cycle_indices;
        for (const auto& cycle_indices : path_with_cycle_indices) {
            if (cycle_indices.size() == 1) {
                path_with_unique_cycle_indices.push_back(*cycle_indices.begin());
            }
        }

        if (path_with_unique_cycle_indices.size() < 2) {
            continue;
        }

        for (size_t i = 0; i < path_with_unique_cycle_indices.size(); ++i) {
            for (size_t j = 0; j < path_with_unique_cycle_indices.size(); ++j) {
                uint32_t element_i = path_with_unique_cycle_indices[i];
                uint32_t element_j = path_with_unique_cycle_indices[j];

                if (i != j && element_i != element_j) {
                    constraints.push_back(std::make_tuple(
                        static_cast<uint32_t>(element_i),
                        static_cast<uint32_t>(element_j)
                    ));
                }
            }
        }
    }

    return constraints;
}

struct TupleHash {
    std::size_t operator()(const std::tuple<uint32_t, uint32_t>& t) const {
        return std::hash<uint32_t>()(std::get<0>(t)) ^ (std::hash<uint32_t>()(std::get<1>(t)) << 1);
    }
};

bool has_cycle(std::unordered_map<std::tuple<uint32_t, uint32_t>, int, TupleHash>& edges_with_weights) {
    // Detect if the directed graph has a cycle using DFS
    std::unordered_map<uint32_t, std::vector<uint32_t>> adj;
    for (const auto& [edge, _] : edges_with_weights) {
        adj[std::get<0>(edge)].push_back(std::get<1>(edge));
    }

    std::unordered_set<uint32_t> visited;
    std::unordered_set<uint32_t> rec_stack;

    std::function<bool(uint32_t)> dfs = [&](uint32_t node) {
        if (rec_stack.count(node)) return true;
        if (visited.count(node)) return false;
        visited.insert(node);
        rec_stack.insert(node);
        for (uint32_t neighbor : adj[node]) {
            if (dfs(neighbor)) return true;
        }
        rec_stack.erase(node);
        return false;
    };

    for (const auto& [node, _] : adj) {
        if (dfs(node)) return true;
    }
    return false;
}

// Returns the edge with the smalles value = {edges_occurence_count_in_cycles} * {edge_weight}
std::tuple<uint32_t, uint32_t> resolve_cycles_greedy_best_edge(
    std::unordered_map<std::tuple<uint32_t, uint32_t>, int, TupleHash>& edges_with_weights,
    std::vector<std::vector<uint64_t>>& cycles
) {
    // Count edge occurrences in cycles
    std::unordered_map<std::tuple<uint32_t, uint32_t>, int, TupleHash> edge_occurence_count_in_cycles;
    for (const auto& cycle : cycles) {
        for (size_t i = 0; i + 1 < cycle.size(); ++i) {
            uint32_t from = static_cast<uint32_t>(cycle[i]);
            uint32_t to = static_cast<uint32_t>(cycle[i + 1]);
            std::tuple<uint32_t, uint32_t> edge(from, to);
            edge_occurence_count_in_cycles[edge]++;
        }
        // If cycle is closed, add edge from last to first
        if (!cycle.empty()) {
            uint32_t from = static_cast<uint32_t>(cycle.back());
            uint32_t to = static_cast<uint32_t>(cycle.front());
            std::tuple<uint32_t, uint32_t> edge(from, to);
            edge_occurence_count_in_cycles[edge]++;
        }
    }

    // Find edge with minimal value = edge_weight * edge_occurence_count_in_cycles
    std::tuple<uint32_t, uint32_t> best_edge;
    int min_value = std::numeric_limits<int>::max();
    for (const auto& [edge, weight] : edges_with_weights) {
        int occur = edge_occurence_count_in_cycles[edge];
        int value = weight * occur;
        if (value < min_value) {
            min_value = value;
            best_edge = edge;
        }
    }
    return best_edge;
}

void resolve_cycles_greedy(
    std::vector<std::tuple<uint32_t, uint32_t>>& constraints,
    std::vector<std::vector<uint64_t>>& cycles
) {
    std::unordered_map<std::tuple<uint32_t, uint32_t>, int, TupleHash> edges_with_weights;
    for (const auto& constraint : constraints) {
        edges_with_weights[constraint]++;
    }

    std::vector<std::tuple<uint32_t, uint32_t>> removed_edges;

    while (has_cycle(edges_with_weights)) {
        auto edge_to_remove = resolve_cycles_greedy_best_edge(edges_with_weights, cycles);
        removed_edges.push_back(edge_to_remove);
        edges_with_weights.erase(edge_to_remove);
        std::cout << "Removed edge: (" << std::get<0>(edge_to_remove) << ", " << std::get<1>(edge_to_remove) << ")" << std::endl;
    }

    // Remove constraints that were removed
    std::vector<std::tuple<uint32_t, uint32_t>> filtered_constraints;
    std::unordered_set<std::tuple<uint32_t, uint32_t>, TupleHash> removed_set(removed_edges.begin(), removed_edges.end());
    for (const auto& constraint : constraints) {
        if (removed_set.find(constraint) == removed_set.end()) {
            filtered_constraints.push_back(constraint);
        }
    }
    constraints = std::move(filtered_constraints);
}

// Helper function to perform all possible topological sorts
std::vector<std::vector<int32_t>> apply_topological_sort(
    std::vector<int32_t>& possible_start_nodes,
    std::unordered_map<std::pair<int32_t, int32_t>, int, TupleHash>& edges
) {
    std::vector<std::vector<int32_t>> result;

    if (possible_start_nodes.empty()) {
        std::vector<int32_t> empty_vec;
        result.push_back(empty_vec);
        return result;
    }

    for (size_t idx = 0; idx < possible_start_nodes.size(); ++idx) {
        int32_t start_node = possible_start_nodes[idx];

        // Prepare new possible_start_nodes
        std::vector<int32_t> new_possible_start_nodes;
        for (size_t i = 0; i < possible_start_nodes.size(); ++i) {
            if (possible_start_nodes[i] != start_node) {
                new_possible_start_nodes.push_back(possible_start_nodes[i]);
            }
        }
        // Add nodes that become available after removing start_node
        for (const auto& edge : edges) {
            if (edge.first.first == start_node) {
                int32_t dest = edge.first.second;
                // Check if all incoming edges to dest are removed
                bool has_other_incoming = false;
                for (const auto& e : edges) {
                    if (e.first.second == dest && e.first.first != start_node) {
                        has_other_incoming = true;
                        break;
                    }
                }
                if (!has_other_incoming) {
                    if (std::find(new_possible_start_nodes.begin(), new_possible_start_nodes.end(), dest) == new_possible_start_nodes.end()) {
                        new_possible_start_nodes.push_back(dest);
                    }
                }
            }
        }

        // Prepare new edges map
        std::unordered_map<std::pair<int32_t, int32_t>, int, TupleHash> new_edges;
        for (const auto& edge : edges) {
            if (edge.first.first != start_node && edge.first.second != start_node) {
                new_edges[edge.first] = edge.second;
            }
        }

        auto solutions_rec = apply_topological_sort(new_possible_start_nodes, new_edges);
        for (const auto& solution_rec : solutions_rec) {
            std::vector<int32_t> solution = solution_rec;
            solution.insert(solution.begin(), start_node);
            result.push_back(solution);
        }
    }

    return result;
}

std::vector<std::vector<int32_t>> solve_constraints_with_topological_sort(
    std::vector<std::tuple<uint32_t, uint32_t>> constraints,
    std::vector<int32_t> nodes
) {
    // Build edges map
    std::unordered_map<std::pair<int32_t, int32_t>, int, TupleHash> edges;
    // std::unordered_set<int32_t> node_set;
    for (const auto& constraint : constraints) {
        int32_t source = static_cast<int32_t>(std::get<0>(constraint));
        int32_t dest = static_cast<int32_t>(std::get<1>(constraint));
        edges[{source, dest}]++;
        // node_set.insert(source);
        // node_set.insert(dest);
    }
    // std::vector<int32_t> nodes(node_set.begin(), node_set.end());

    // Find possible start nodes (nodes with no incoming edges)
    std::vector<int32_t> possible_start_nodes;
    for (int32_t node : nodes) {
        bool has_incoming = false;
        for (const auto& constraint : constraints) {
            if (std::get<1>(constraint) == node) {
                has_incoming = true;
                break;
            }
        }
        if (!has_incoming) {
            possible_start_nodes.push_back(node);
        }
    }
    std::cout << "Node set: ";
    for (const auto& node : nodes) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
    auto all_orders = apply_topological_sort(possible_start_nodes, edges);
    std::cout << "All possible topological orders:" << std::endl;
    for (const auto& order : all_orders) {
        std::cout << "[ ";
        for (size_t i = 0; i < order.size(); ++i) {
            std::cout << order[i];
            if (i + 1 < order.size()) std::cout << ", ";
        }
        std::cout << " ]" << std::endl;
    }
    return all_orders;
}

std::vector<std::vector<int32_t>> order_cycles(
    const Graph& graph,
    std::vector<Jump>& jumps,
    std::vector<std::vector<uint64_t>>& cycles
) {
    auto node_to_cycle_map = get_node_to_cycle_map(cycles);
    auto all_cycle_indices = get_all_cycle_indices(node_to_cycle_map);
    auto constraints = generate_constraints(graph, jumps, node_to_cycle_map);
    std::cout << "Constraints before resolving cycles: " << constraints.size() << std::endl;
    resolve_cycles_greedy(constraints, cycles);
    std::cout << "Constraints after resolving cycles: " << constraints.size() << std::endl;
    return solve_constraints_with_topological_sort(constraints, all_cycle_indices);
}
