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
        if (graph.nodes.find(jump.start_k_mer_id) != graph.nodes.end()
        || graph.nodes.find(jump.end_k_mer_id) != graph.nodes.end()) {
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
            sub_path.push_back(start);
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

vector<tuple<uint32_t, uint32_t>> generate_out_of_cycles_constraints_from_jump(
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
        bool cycle_index_differs = false;
        optional<uint32_t> common_cycle_index = std::nullopt;
        for (const auto& path : all_paths) {
            const uint64_t node = path[i];

            uint32_t cycle_index;
            auto it = node_to_cycle_map.find(node);
            if (it != node_to_cycle_map.end()) {
                cycle_index = it->second;
            } else {
                cycle_index = NOT_IN_ANY_CYCLE_INDEX;
            }

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

    if (cycle_indices_in_order.size() < 2) {
        return {};
    }

    // Be z the index for the "not being in any cycle" and X a sequence of cycle indices
    // The following cases are possible for cycle_indices_in_order:
    //  - z X => Indicates the indices X are close to the start, i.e. are the start
    //  - X z => Indicates the indices X are close to the end, i.e. are the end
    //  - z X z => Doesn't provide any value
    //  - X z X => Doesn't provide any value
    if (cycle_indices_in_order[0] == NOT_IN_ANY_CYCLE_INDEX
    && cycle_indices_in_order[cycle_indices_in_order.size() - 1] != NOT_IN_ANY_CYCLE_INDEX) {
        bool is_valid_start = true;
        bool expecting_z = true;
        for (const auto& cycle_index : cycle_indices_in_order) {
            if (expecting_z && cycle_index != NOT_IN_ANY_CYCLE_INDEX) {
                expecting_z = false;
            } else if (!expecting_z && cycle_index == NOT_IN_ANY_CYCLE_INDEX) {
                is_valid_start = false;
                break;
            }
        }

        if (is_valid_start) {
            return every_possible_combination(cycle_indices_in_order);
        }
    } else if (cycle_indices_in_order[0] != NOT_IN_ANY_CYCLE_INDEX
    && cycle_indices_in_order[cycle_indices_in_order.size() - 1] == NOT_IN_ANY_CYCLE_INDEX) {
        bool is_valid_end = true;
        bool expecting_z = false;
        for (const auto& cycle_index : cycle_indices_in_order) {
            if (expecting_z && cycle_index != NOT_IN_ANY_CYCLE_INDEX) {
                is_valid_end = false;
                break;
            } else if (!expecting_z && cycle_index == NOT_IN_ANY_CYCLE_INDEX) {
                expecting_z = true;
            }
        }

        if (is_valid_end) {
            return every_possible_combination(cycle_indices_in_order);
        }
    }
    
    return {};
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
        const auto extra_constraints = generate_out_of_cycles_constraints_from_jump(
            graph,
            jump,
            node_to_cycle_map
        );

        for (const auto& constraint : jump_constraints) {
            constraints.push_back(constraint);
        }
        for (const auto& constraint : extra_constraints) {
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

pair<tuple<uint32_t, uint32_t>, float> resolve_cycles_greedy_best_edge(
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash>& edges_with_weights
) {
    tuple<uint32_t, uint32_t> best_edge;
    int total_weight = 0;
    int min_weight = std::numeric_limits<int>::max();
    for (const auto& [edge, weight] : edges_with_weights) {
        total_weight += weight;
        if (weight < min_weight) {
            min_weight = weight;
            best_edge = edge;
        }
    }

    float confidence = static_cast<float>(total_weight - min_weight);
    confidence /= static_cast<float>(total_weight);

    return std::make_pair(best_edge, confidence);
}

void resolve_cycles_greedy(
    vector<tuple<uint32_t, uint32_t>>& constraints,
    float& confidence
) {
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash> edges_with_weights;
    unordered_map<uint32_t, vector<uint32_t>> edges;
    for (const auto& constraint : constraints) {
        uint32_t from = std::get<0>(constraint);
        uint32_t to = std::get<1>(constraint);
        
        if (from == NOT_IN_ANY_CYCLE_INDEX || to == NOT_IN_ANY_CYCLE_INDEX) {
            continue;
        }

        edges_with_weights[constraint]++;
        edges[from].push_back(to);
    }

    vector<tuple<uint32_t, uint32_t>> removed_edges;

    while (has_cycle(edges)) {
        auto result = resolve_cycles_greedy_best_edge(edges_with_weights);
        auto edge_to_remove = std::get<0>(result);
        auto edge_remove_confidence = std::get<1>(result);

        confidence *= edge_remove_confidence;

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
    const unordered_map<uint32_t, int>& node_affection_to_start,
    const unordered_map<uint32_t, vector<uint32_t>>& edges,
    float& confidence,
    vector<uint32_t>& total_order
) {
    if (possible_start_nodes.empty()) {
        return;
    }

    vector<uint32_t> new_possible_start_nodes = possible_start_nodes;
    unordered_map<uint32_t, vector<uint32_t>> new_edges = edges;
    
    // Choose start_node by the affection to the start
    int best_start_node = 0;
    int best_affection = node_affection_to_start.at(new_possible_start_nodes[best_start_node]);
    
    for (int i = 0; i < new_possible_start_nodes.size(); ++i) {
        int current_affection = node_affection_to_start.at(new_possible_start_nodes[i]);
        if (current_affection >= best_affection) {
            best_affection = current_affection;
            best_start_node = i;
        }
    }

    int second_best_affection = node_affection_to_start.at(new_possible_start_nodes[
        (best_start_node + 1) % new_possible_start_nodes.size()
    ]);
    for (int i = 0; i < new_possible_start_nodes.size(); ++i) {
        int current_affection = node_affection_to_start.at(new_possible_start_nodes[i]);
        if (i != best_start_node && current_affection >= second_best_affection) {
            second_best_affection = current_affection;
        }
    }

    // Using some function f: [1.0;inf] -> [0.5; 1.0]
    // with f(1.0) = 0.5 and for greater x it gets closer to 1.0
    // e.g. f(x) = -1 / (x**6 + 1) + 1
    int shift = 0; // affection can be negative
    if (second_best_affection <= 0) {
        shift = -second_best_affection + 1;
    }

    float ratio = static_cast<float>(best_affection + shift) / static_cast<float>(second_best_affection + shift);
    float confidence_for_best = (-1.0 / (pow(ratio, 6) + 1)) + 1;
    if (new_possible_start_nodes.size() == 1) {
        confidence_for_best = 1.0;
    }
    confidence *= confidence_for_best;


    uint32_t start_node = new_possible_start_nodes[best_start_node];
    total_order.push_back(start_node);

    // Remove the start_node as choosable (by index, not value)
    new_possible_start_nodes.erase(new_possible_start_nodes.begin() + best_start_node);

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

    apply_topological_sort(
        new_possible_start_nodes,
        node_affection_to_start,
        new_edges,
        confidence,
        total_order
    );
}

vector<uint32_t> solve_constraints_with_topological_sort(
    const vector<tuple<uint32_t, uint32_t>>& constraints,
    const vector<uint32_t>& nodes,
    float& confidence
) {
    unordered_map<uint32_t, vector<uint32_t>> edges;
    for (const auto& constraint : constraints) {
        const uint32_t source = std::get<0>(constraint);
        const uint32_t dest = std::get<1>(constraint);

        if (source == NOT_IN_ANY_CYCLE_INDEX || dest == NOT_IN_ANY_CYCLE_INDEX) {
            continue;
        }

        auto edge_search = edges.find(source);
        if (edge_search != edges.end()) {
            if (std::find(edge_search->second.begin(), edge_search->second.end(), dest)
            == edge_search->second.end()) {
                edge_search->second.push_back(dest);
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
            const uint32_t source = std::get<0>(constraint);
            const uint32_t dest = std::get<1>(constraint);

            if (source != NOT_IN_ANY_CYCLE_INDEX && dest == node) {
                has_incoming = true;
                break;
            }
        }
        if (!has_incoming) {
            possible_start_nodes.push_back(node);
        }
    }

    // Used for the heuristic, which start node to choose from
    // Big positive number indicates a strong affection towards being the first choosen start node
    // Big negative, affection towards the end
    unordered_map<uint32_t, int> node_affection_to_start;
    for (uint32_t node : nodes) {
        node_affection_to_start[node] = 0;
    }

    for (const auto& constraint : constraints) {
        const uint32_t source = std::get<0>(constraint);
        const uint32_t dest = std::get<1>(constraint);

        if (source != NOT_IN_ANY_CYCLE_INDEX && dest != NOT_IN_ANY_CYCLE_INDEX) {
            continue;
        }

        if (source == NOT_IN_ANY_CYCLE_INDEX) {
            node_affection_to_start[dest]++;
        } else {
            node_affection_to_start[source]--;
        }
    }

    // Sort possible_start_nodes by affection descending before printing
    std::vector<std::pair<uint32_t, int>> sorted_affection;
    for (const auto& [node, affection] : node_affection_to_start) {
        sorted_affection.emplace_back(node, affection);
    }
    std::sort(sorted_affection.begin(), sorted_affection.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

    vector<uint32_t> total_order;
    apply_topological_sort(
        possible_start_nodes,
        node_affection_to_start,
        edges,
        confidence,
        total_order
    );

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

    float resolve_cycles_confidence = 1.0;
    resolve_cycles_greedy(constraints, resolve_cycles_confidence);
    confidence *= resolve_cycles_confidence;

    std::cout << "      ▸ " << constraints.size();
    std::cout << " constraints remain after resolving cycles.";
    std::cout << "The confidence of the resolution is " << std::fixed << std::setprecision(2);
    std::cout << (resolve_cycles_confidence * 100) << "%" << endl;

    return solve_constraints_with_topological_sort(constraints, all_cycle_indices, confidence);
}

vector<vector<uint64_t>> get_ordered_cycles(
    const vector<uint32_t>& cycle_order,
    const vector<vector<uint64_t>>& cycles
) {
    vector<vector<uint64_t>> ordered_cycles;

    for (const auto& cycle_index : cycle_order) {
        if (cycle_index < cycles.size()) {
            ordered_cycles.push_back(cycles[cycle_index]);
        }
    }

    return ordered_cycles;
}
