#include "spacer_ordering.h"

class TarjanSCC {
private:
    const SDBG& sdbg;
    unordered_map<uint64_t, int> index_map;
    unordered_map<uint64_t, int> lowlink_map;
    unordered_set<uint64_t> on_stack;
    std::stack<uint64_t> stack;
    vector<vector<uint64_t>> components;
    int index_counter;

    void strongconnect(const uint64_t node) {
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
                    const uint64_t neighbor = outgoings[i];
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
            vector<uint64_t> component;
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
    TarjanSCC(const SDBG& graph) : sdbg(graph), index_counter(0) {}

    vector<vector<uint64_t>> find_components() {
        // Find all valid nodes first
        vector<uint64_t> valid_nodes;
        for (uint64_t node = 0; node < sdbg.size(); ++node) {
            if (sdbg.IsValidEdge(node)) {
                valid_nodes.push_back(node);
            }
        }

        // Run Tarjan's algorithm on all unvisited valid nodes
        for (const uint64_t node : valid_nodes) {
            if (index_map.find(node) == index_map.end()) {
                strongconnect(node);
            }
        }

        return components;
    }
};

void keep_crispr_regions(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& cycles
) {
    unordered_set<uint64_t> cycle_nodes;
    for (const auto& cycle : cycles) {
        for (uint64_t node : cycle) {
            cycle_nodes.insert(node);
        }
    }

    for (uint64_t node_id = 0; node_id < sdbg.size(); ++node_id) {
        if (sdbg.IsValidEdge(node_id)) {
            // If this node is not in any cycle, invalidate it
            if (cycle_nodes.find(node_id) == cycle_nodes.end()) {
                sdbg.SetInvalidEdge(node_id);
            }
        }
    }

    // std::cout << "Set set for crispr regions is ";
    // std::cout << cycle_nodes.size() << " nodes" << endl;

    // size_t valid_edge_count = 0;
    // for (uint64_t node_id = 0; node_id < sdbg.size(); ++node_id) {
    //     if (sdbg.IsValidEdge(node_id)) {
    //         valid_edge_count++;
    //     }
    // }
    // std::cout << "keep_crispr_regions sdbg has " << valid_edge_count;
    // std::cout << " valid edges remaining" << endl;
}

vector<Graph> divide_graph_into_subgraphs(const SDBG& sdbg) {
    vector<Graph> subgraphs;
    
    // Find strongly connected components using Tarjan's algorithm
    TarjanSCC tarjan(sdbg);
    auto components = tarjan.find_components();

    // std::cout << "Components: [";
    // for (size_t i = 0; i < components.size(); ++i) {
    //     vector<uint64_t> sorted_component = components[i];
    //     std::sort(sorted_component.begin(), sorted_component.end());
    //     std::cout << "[";
    //     for (size_t j = 0; j < sorted_component.size(); ++j) {
    //         std::cout << sorted_component[j];
    //         if (j + 1 < sorted_component.size()) std::cout << ",";
    //     }
    //     std::cout << "]";
    //     if (i + 1 < components.size()) std::cout << ", ";
    // }
    // std::cout << "]" << endl;
    
    // std::cout << "Found " << components.size() << " strongly connected components" << endl;
    
    // Convert each component to a SubGraph
    for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
        const auto& component = components[comp_idx];
        Graph subgraph;
        
        // Add all edges within this component
        for (const uint64_t node : component) {
            if (!sdbg.IsValidEdge(node)) continue;
            
            int outdegree = sdbg.EdgeOutdegree(node);
            if (outdegree > 0) {
                uint64_t* outgoings = new uint64_t[outdegree];
                if (sdbg.OutgoingEdges(node, outgoings) != -1) {
                    for (int i = 0; i < outdegree; ++i) {
                        const uint64_t neighbor = outgoings[i];
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

vector<Graph> get_crispr_regions(
    SDBG& sdbg,
    unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    vector<vector<uint64_t>> cycles;
    for (const auto& [_start_node, cycles_of_start_node] : cycles_map) {
        for (const auto& cycle : cycles_of_start_node) {
            cycles.push_back(cycle);
        }
    }

    keep_crispr_regions(sdbg, cycles);
    return divide_graph_into_subgraphs(sdbg);
}

vector<Jump> get_relevant_jumps(const Graph& graph, const vector<Jump>& all_jumps) {
    vector<Jump> relevant_jumps;

    for (const auto& jump : all_jumps) {
        if (graph.nodes.find(jump.start_k_mer_id) != graph.nodes.end() &&
            graph.nodes.find(jump.end_k_mer_id) != graph.nodes.end()) {
            relevant_jumps.push_back(jump);
        }
    }

    return relevant_jumps;
}

vector<vector<uint64_t>> get_relevant_cycles(
    const Graph& graph,
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& all_cycles_map
) {
    vector<vector<uint64_t>> relevant_cycles;

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

unordered_map<uint64_t, uint32_t> get_node_to_unique_cycle_map(
    const vector<vector<uint64_t>>& cycles
) {
    vector<unordered_set<uint64_t>> cycle_node_sets;
    for (const auto& cycle : cycles) {
        cycle_node_sets.emplace_back(cycle.begin(), cycle.end());
    }

    unordered_map<uint64_t, uint32_t> node_to_cycle_map;
    for (uint32_t i = 0; i < cycle_node_sets.size(); ++i) {
        // Build union of all other sets
        unordered_set<uint64_t> other_nodes;
        for (uint32_t j = 0; j < cycle_node_sets.size(); ++j) {
            if (j == i) continue;
            other_nodes.insert(cycle_node_sets[j].begin(), cycle_node_sets[j].end());
        }
        // Find unique nodes in cycle i
        for (const uint64_t node : cycle_node_sets[i]) {
            if (other_nodes.find(node) == other_nodes.end()) {
                node_to_cycle_map[node] = i;
            }
        }
    }
    return node_to_cycle_map;
}

vector<uint32_t> get_all_cycle_indices(
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<uint32_t> cycle_indices;

    for (const auto& [node, cycle_idx] : node_to_cycle_map) {
        if (std::find(cycle_indices.begin(), cycle_indices.end(), cycle_idx) == cycle_indices.end()) {
            cycle_indices.push_back(cycle_idx);
        }
    }

    return cycle_indices;
}

vector<vector<uint64_t>> find_all_possible_paths(
    const Graph& graph,
    uint64_t start,
    uint64_t end,
    uint64_t nodes_left
) {
    vector<vector<uint64_t>> paths;

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

vector<tuple<uint32_t, uint32_t>> every_possible_combination(const vector<uint32_t>& v) {
    vector<tuple<uint32_t, uint32_t>> possible_combination;
    
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = i + 1; j < v.size(); ++j) {
            const uint32_t element_i = v[i];
            const uint32_t element_j = v[j];

            if (element_i != element_j) {
                possible_combination.push_back(std::make_tuple(
                    static_cast<uint32_t>(element_i),
                    static_cast<uint32_t>(element_j)
                ));
            }
        }
    }

    return possible_combination;
}

vector<tuple<uint32_t, uint32_t>> generate_constraints_from_jump(
    const Graph& graph,
    const Jump& jump,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<tuple<uint32_t, uint32_t>> constraints;

    const vector<vector<uint64_t>> all_paths = find_all_possible_paths(
        graph,
        jump.start_k_mer_id,
        jump.end_k_mer_id,
        jump.nodes_in_between + 1 // + 1, to count the end node as well
    );

    if (all_paths.empty()) {
        return constraints;
    }

    vector<uint32_t> cycle_indices_in_order;
    // + 2 as start and end nodes are included in the path
    for (size_t i = 0; i < static_cast<size_t>(jump.nodes_in_between) + 2; ++i) {
        bool all_nodes_can_be_mapped = true;
        for (const auto& path : all_paths) {
            const uint64_t node = path[i];
            if (node_to_cycle_map.find(node) == node_to_cycle_map.end()) {
                all_nodes_can_be_mapped = false;
                break;
            }
        }

        if (!all_nodes_can_be_mapped) {
            continue;
        }

        bool cycle_index_differs = false;
        optional<uint32_t> common_cycle_index = std::nullopt;
        for (const auto& path : all_paths) {
            const uint64_t node = path[i];
            const uint32_t cycle_index = node_to_cycle_map.find(node)->second;

            if (common_cycle_index.has_value() && common_cycle_index.value() != cycle_index) {
                cycle_index_differs = true;
                break;
            } else if (!common_cycle_index.has_value()) {
                common_cycle_index = std::make_optional(cycle_index);
            }
        }

        if (cycle_index_differs || !common_cycle_index.has_value()) {
            continue;
        }

        cycle_indices_in_order.push_back(common_cycle_index.value());
    }

    return every_possible_combination(cycle_indices_in_order);
}

vector<tuple<uint32_t, uint32_t>> generate_constraints(
    const Graph& graph,
    const vector<Jump>& jumps,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
) {
    vector<tuple<uint32_t, uint32_t>> constraints;

    for (const auto& jump : jumps) {
        const auto jump_constraints = generate_constraints_from_jump(
            graph,
            jump,
            node_to_cycle_map
        );

        for (const auto& constraint : jump_constraints) {
            constraints.push_back(constraint);
        }
    }

    return constraints;
}

bool has_cycle(const unordered_map<uint32_t, vector<uint32_t>>& edges) {
    unordered_set<uint32_t> visited;
    unordered_set<uint32_t> rec_stack;

    std::function<bool(uint32_t)> dfs = [&](uint32_t node) {
        if (rec_stack.count(node)) return true;
        if (visited.count(node)) return false;
        visited.insert(node);
        rec_stack.insert(node);
        auto it = edges.find(node);
        if (it != edges.end()) {
            for (uint32_t neighbor : it->second) {
                if (dfs(neighbor)) return true;
            }
        }
        rec_stack.erase(node);
        return false;
    };

    for (const auto& [node, _] : edges) {
        if (dfs(node)) return true;
    }
    return false;
}

tuple<uint32_t, uint32_t> resolve_cycles_greedy_best_edge(
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash>& edges_with_weights,
    const vector<vector<uint64_t>>& cycles
) {
    // Count edge occurrences in cycles
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash> edge_occurence_count_in_cycles;
    for (const auto& cycle : cycles) {
        for (size_t i = 0; i + 1 < cycle.size(); ++i) {
            uint32_t from = static_cast<uint32_t>(cycle[i]);
            uint32_t to = static_cast<uint32_t>(cycle[i + 1]);
            tuple<uint32_t, uint32_t> edge(from, to);
            edge_occurence_count_in_cycles[edge]++;
        }
        // If cycle is closed, add edge from last to first
        if (!cycle.empty()) {
            uint32_t from = static_cast<uint32_t>(cycle.back());
            uint32_t to = static_cast<uint32_t>(cycle.front());
            tuple<uint32_t, uint32_t> edge(from, to);
            edge_occurence_count_in_cycles[edge]++;
        }
    }

    // Find edge with minimal value = edge_weight * edge_occurence_count_in_cycles
    tuple<uint32_t, uint32_t> best_edge;
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
    vector<tuple<uint32_t, uint32_t>>& constraints,
    const vector<vector<uint64_t>>& cycles
) {
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash> edges_with_weights;
    unordered_map<uint32_t, vector<uint32_t>> edges;
    for (const auto& constraint : constraints) {
        edges_with_weights[constraint]++;

        uint32_t from = std::get<0>(constraint);
        uint32_t to = std::get<1>(constraint);
        edges[from].push_back(to);
    }

    vector<tuple<uint32_t, uint32_t>> removed_edges;

    while (has_cycle(edges)) {
        auto edge_to_remove = resolve_cycles_greedy_best_edge(edges_with_weights, cycles);
        removed_edges.push_back(edge_to_remove);
        edges_with_weights.erase(edge_to_remove);
        
        uint32_t from = std::get<0>(edge_to_remove);
        uint32_t to = std::get<1>(edge_to_remove);
        auto it = edges.find(from);
        if (it != edges.end()) {
            auto& vec = it->second;
            vec.erase(std::remove(vec.begin(), vec.end(), to), vec.end());
            if (vec.empty()) {
                edges.erase(it);
            }
        }
    }

    const unordered_set<tuple<uint32_t, uint32_t>, TupleHash> removed_set(
        removed_edges.begin(),
        removed_edges.end()
    );
    
    vector<tuple<uint32_t, uint32_t>> filtered_constraints;
    for (const auto& constraint : constraints) {
        if (removed_set.find(constraint) == removed_set.end()) {
            filtered_constraints.push_back(constraint);
        }
    }

    constraints = std::move(filtered_constraints);
}

void apply_topological_sort(
    const vector<uint32_t>& possible_start_nodes,
    const unordered_map<uint32_t, vector<uint32_t>>& edges,
    uint32_t& possible_choices,
    vector<uint32_t>& total_order
) {
    if (possible_start_nodes.empty()) {
        return;
    }

    vector<uint32_t> new_possible_start_nodes = possible_start_nodes;
    unordered_map<uint32_t, vector<uint32_t>> new_edges = edges;
    
    // Choose some start_node
    uint32_t start_node = new_possible_start_nodes[0];
    possible_choices *= new_possible_start_nodes.size();
    total_order.push_back(start_node);

    // Remove the start_node as choosable
    new_possible_start_nodes.erase(
        std::remove(new_possible_start_nodes.begin(), new_possible_start_nodes.end(), start_node),
        new_possible_start_nodes.end()
    );

    // Removes every edge containing the start_node
    vector<uint32_t> start_node_candidates;
    auto it = new_edges.find(start_node);
    if (it != new_edges.end()) {
        for (const auto& node : it->second) {
            start_node_candidates.push_back(node);
        }
    }
    for (auto it = new_edges.begin(); it != new_edges.end(); ) {
        if (it->first == start_node) {
            it = new_edges.erase(it);
        } else {
            // Remove any occurrence of start_node in the value vector
            auto& vec = it->second;
            vec.erase(std::remove(vec.begin(), vec.end(), start_node), vec.end());
            ++it;
        }
    }

    // Finding new possible start nodes
    for (const auto& start_node_candidate : start_node_candidates) {
        bool has_other_incoming = false;
        for (const auto& edge : new_edges) {
            bool destinations_contain_candidate = std::find(edge.second.begin(), edge.second.end(), start_node_candidate) != edge.second.end();
            if (destinations_contain_candidate) {
                has_other_incoming = true;
                break;
            }
        }
        if (!has_other_incoming) {
            new_possible_start_nodes.push_back(start_node_candidate);
        }
    }

    apply_topological_sort(new_possible_start_nodes, new_edges, possible_choices, total_order);
}

vector<uint32_t> solve_constraints_with_topological_sort(
    const vector<tuple<uint32_t, uint32_t>>& constraints,
    const vector<uint32_t>& nodes,
    float& confidence
) {
    unordered_map<uint32_t, vector<uint32_t>> edges;
    for (const auto& constraint : constraints) {
        uint32_t source = static_cast<uint32_t>(std::get<0>(constraint));
        uint32_t dest = static_cast<uint32_t>(std::get<1>(constraint));
        auto it = edges.find(source);
        if (it != edges.end()) {
            if (std::find(it->second.begin(), it->second.end(), source) == it->second.end()) {
            it->second.push_back(dest);
            }
        } else {
            edges[source] = {dest};
        }
    }

    // Find possible start nodes (nodes with no incoming edges)
    vector<uint32_t> possible_start_nodes;
    for (uint32_t node : nodes) {
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

    // std::cout << "Node set: ";
    // for (const auto& node : nodes) {
    //     std::cout << node << " ";
    // }
    // std::cout << endl;

    vector<uint32_t> total_order;
    uint32_t possible_choices = 1;
    apply_topological_sort(possible_start_nodes, edges, possible_choices, total_order);

    confidence = 1.0f / possible_choices;

    // std::cout << "All possible topological orders:" << endl;
    // for (const auto& order : all_orders) {
    //     std::cout << "[ ";
    //     for (size_t i = 0; i < order.size(); ++i) {
    //         std::cout << order[i];
    //         if (i + 1 < order.size()) std::cout << ", ";
    //     }
    //     std::cout << " ]" << endl;
    // }

    return total_order;
}

vector<uint32_t> order_cycles(
    const Graph& graph,
    const vector<Jump>& jumps,
    const vector<vector<uint64_t>>& cycles,
    float& confidence
) {
    const auto node_to_cycle_map = get_node_to_unique_cycle_map(cycles);
    const auto all_cycle_indices = get_all_cycle_indices(node_to_cycle_map);
    auto constraints = generate_constraints(graph, jumps, node_to_cycle_map);

    std::cout << "      ▸ " << constraints.size() << " constraints derived" << std::endl;

    resolve_cycles_greedy(constraints, cycles);

    std::cout << "      ▸ " << constraints.size();
    std::cout << " constraints remain after resolving cycles" << endl;

    return solve_constraints_with_topological_sort(constraints, all_cycle_indices, confidence);
}

vector<uint64_t> turn_cycle_order_into_node_order(
    const vector<uint32_t>& cycle_order,
    const vector<vector<uint64_t>>& cycles
) {
    vector<uint64_t> node_order;

    for (const auto& cycle_index : cycle_order) {
        if (cycle_index < cycles.size()) {
            for (const auto& node_id : cycles[cycle_index]) {
                node_order.push_back(node_id);
            }
        }
    }

    return node_order;
}
