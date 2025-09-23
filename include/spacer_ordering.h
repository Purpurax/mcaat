/**
 * @file spacer_ordering.h
 * @brief Functions for ordering cycles and regions in Succinct De Bruijn Graphs (SDBG).
 *
 * Contains:
 * - Graph structure for basic graph operations.
 * - Functions to extract and order cycles and regions from SDBG.
 * - Functions to get relevant jumps and cycles.
 * - Functions to combine cycles into ordered node sequences.
 *
 * Main public functions:
 * - get_crispr_regions_extended_by_k
 * - get_relevant_jumps
 * - get_relevant_cycles
 * - order_cycles
 * - turn_cycle_order_into_node_order
 */
#ifndef INCLUDE_SPACER_ORDERING_H_
#define INCLUDE_SPACER_ORDERING_H_

#include <unordered_set>
#include <stack>
#include <algorithm>

#include "sdbg/sdbg.h"
#include "jumps.h"

#include "core/cft.hpp"
#include "core/Instance.hpp"
#include "algorithms/Refinement.hpp"

using namespace std;

/**
 * @brief A basic graph with an adjacency_list and nodes
 */
struct Graph {
    unordered_map<uint64_t, vector<uint64_t>> adjacency_list;
    unordered_set<uint64_t> nodes;

    /**
     * @brief Adds a new edge to the graph
     * 
     * @param from
     * @param to
     */
    void add_edge(uint64_t from, uint64_t to) {
        adjacency_list[from].push_back(to);
        nodes.insert(from);
        nodes.insert(to);
    }

    /**
     * @brief Returns the total number of edges in the graph
     *
     * @return The edge count (size_t)
     */
    size_t edge_count() const {
        size_t count = 0;
        for (const auto& kv : adjacency_list) {
            count += kv.second.size();
        }
        return count;
    }
};

const uint32_t NOT_IN_ANY_CYCLE_INDEX = std::numeric_limits<uint32_t>::max();

/**
 * @internal
 * @brief Used to get a hash from a tuple
 */
struct TupleHash {
    size_t operator()(const tuple<uint32_t, uint32_t>& t) const {
        return std::hash<uint32_t>()(std::get<0>(t)) ^ (std::hash<uint32_t>()(std::get<1>(t)) << 1);
    }
};

/**
 * @internal
 * @brief Reduces the succinct de bruijn graphs edges to only have edges
 * that can be found in the cycles (and others with a distance of k)
 * 
 * @param sdbg The Succinct De Bruijn Graph that should be modified
 * @param k How many edges it should extend from that region
 * @param cycles The cycles that determine which edge is kept
 */
void keep_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const size_t& k,
    const vector<vector<uint64_t>>& cycles
);

/**
 * @internal
 * @brief Finds the strongly connected components and splits the graph
 * 
 * @todo Simplify the strongly connected components algorithm.
 * @todo Instead of strongly connected components you might want to judge by
 * grouping the cycles together by repeats and extracting the subgraphs by repeats
 * 
 * @param sdbg The Succinct De Bruijn Graph used for the subgraphs
 * 
 * @return The subgraphs (vector<Graph>)
 */
vector<Graph> divide_graph_into_subgraphs(const SDBG& sdbg);

/**
 * @brief Creates the subgraphs from the total graph reduced to the cycles extended by k
 * 
 * @param sdbg The Succinct De Bruijn Graph
 * @param k Extend the region by k edges in all directions
 * @param cycles The cycles outlining the subgraphs
 * 
 * @return The subgraphs (vector<Graph>)
 */
vector<Graph> get_crispr_regions_extended_by_k(
    SDBG& sdbg,
    const size_t& k,
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles
);

/**
 * @internal
 * @brief Gets all cycle indices that occur in the map
 * 
 * @param node_to_cycle_map The map used
 * 
 * @return All Cycle indices (vector<uint32_t>)
 */
vector<uint32_t> get_all_cycle_indices(
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
);

/**
 * @brief Get the jumps where the end- and start-kmers are part of the graphs nodes
 * 
 * @param graph The graph used to check the nodes
 * @param jumps The jumps from which the relevant ones are copied
 * 
 * @return The relevant jumps (vector<Jump>)
 */
vector<Jump> get_relevant_jumps(const Graph& graph, const vector<Jump>& jumps);

/**
 * @brief Get the relevant cycles
 * 
 * Cycles are relevant, if they all of their nodes are part of the graph
 * 
 * @todo change parameter all_cycles_map to a more meaningful type
 * 
 * @param graph The graph used to check the nodes
 * @param all_cycles_map The cycles from which the relevant ones are copied
 * 
 * @return The relevant cycles (vector<vector<uint64_t>>)
 */
vector<vector<uint64_t>> get_relevant_cycles(
    const Graph& graph,
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& all_cycles_map
);

/**
 * @brief Discards the cycles that do not contribute new nodes to the total cover
 * 
 * Uses min-set-cover to find the relevant cycles and remove the rest
 * 
 * @param cycles
 */
void get_minimum_cycles_for_full_coverage(vector<vector<uint64_t>>& cycles);

/**
 * @internal
 * @brief Solves the min cover problem
 * 
 * Given a universe and sets, find the minimum amount of sets,
 * such that all elements of the universe are in one of those sets.
 * 
 * @note Keep in mind that the elements of the universe must be from 0..universe.size()
 * And the sets contain which columns they satisfy, so from 0..universe.size() again.
 * 
 * @param universe All elements
 * @param sets Sets that are chosen to cover
 * @return The indices of the sets that are chosen as a solution (vector<size_t>)
 */
vector<size_t> solve_min_cover_problem(
    const unordered_set<uint32_t>& universe,
    const vector<vector<uint32_t>>& sets
);

/**
 * @brief Get the minimum cycles to have a full coverage of nodes
 * 
 * Every cycle must have a unique node and therefore contributes to the coverage.
 * By computing the minimum set cover, only those cycles remain.
 * 
 * @post cycles.size() will be the same or smaller
 * 
 * @param cycles
 */
void get_minimum_cycles_for_full_coverage(vector<vector<uint64_t>>& cycles);

/**
 * @internal
 * @brief Creates a map from nodes to cycle indices
 * 
 * The unique nodes from a cycle compared to other cycles, so the nodes that only exist in one cycle,
 * are mapped to the cycle index where they are unique.
 * If no unqiue nodes exist, the map will be empty.
 * 
 * @param cycles The cycles used to create the map
 * 
 * @return The node_id to cycle_index map (unordered_map<uint64_t, uint32_t>)
 */
unordered_map<uint64_t, uint32_t> get_node_to_unique_cycle_map(
    const vector<vector<uint64_t>>& cycles
);

/**
 * @internal
 * @brief Get every possible ORDERED combination
 * 
 * Gets every possible combination without loosing the order in the process.
 * 
 * @param v Vector to generate possible combinations 
 * @return The possible combinations (vector<tuple<uint32_t, uint32_t>>)
 */
vector<tuple<uint32_t, uint32_t>> every_possible_combination(const vector<uint32_t>& v);

/**
 * @internal
 * @brief Uses a single jump to derive constraints
 * 
 * Uses the information provided by the jump to derive as many constraints as possible.
 * 
 * **Constraint**: is composed of 2 indices, standing for cycle indices.
 * 
 * Internally it tries to reconstruct the path from the start and end node using the jump.
 * With that information specific nodes are mapped to unique cycle indices.
 * This vector of cycle indices can then be used for constraint construction,
 * as every two-value-combination of that vector is a valid constraint.
 * 
 * @param graph The graph used to find the jumps path
 * @param jump The jump for which constraints have to be found
 * @param node_to_cycle_map The map used to turn node_ids into usable cycle indices
 * 
 * @return All valid constraints (vector<tuple<uint32_t, uint32_t>>)
 */
vector<tuple<uint32_t, uint32_t>> generate_constraints_from_jump(
    const Graph& graph,
    const Jump& jump,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
);

/**
 * @internal
 * @brief Uses a single jump to derive constraints involving nodes outside the cycles
 * 
 * This function creates constraints using node ids that cannot be mapped to any cycle index,
 * treating these nodes as representing the "outside" of the CRISPR array.
 * 
 * Assumes that the CRISPR array must have at least one jump entering from outside and one jump leading outside.
 * Constraints are constructed to reflect transitions between the outside and the cycles.
 * 
 * **Constraint**: is composed of 2 indices, standing for cycle indices or a special value representing "outside".
 * 
 * Internally, the function reconstructs the path from the start and end node using the jump.
 * Nodes that cannot be mapped to any cycle index are considered outside.
 * Constraints are generated for transitions from outside to inside and from inside to outside.
 * 
 * @param graph The graph used to find the jump's path
 * @param jump The jump for which constraints have to be found
 * @param node_to_cycle_map The map used to turn node_ids into usable cycle indices
 * 
 * @return All valid constraints involving outside nodes (vector<tuple<uint32_t, uint32_t>>)
 */
vector<tuple<uint32_t, uint32_t>> generate_out_of_cycles_constraints_from_jump(
    const Graph& graph,
    const Jump& jump,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
);

/**
 * @internal
 * @brief Uses the jumps to derive constraints
 * 
 * A Constraint is composed of 2 indices, standing for cycle indices.
 * Uses the information provided by the jumps to derive as many constraints as possible.
 * Also extra constraints that use the uint32_t max as an index to indicate outside of cycles
 * 
 * For more technical information look into:
 * @see generate_constraints_from_jump
 * 
 * @param graph The graph used to find the jumps path
 * @param jumps All jumps for which constraints shall be derived
 * @param node_to_cycle_map The map used to turn node_ids into usable cycle indices
 * 
 * @return All valid constraints (vector<tuple<uint32_t, uint32_t>>)
 */
vector<tuple<uint32_t, uint32_t>> generate_constraints(
    const Graph& graph,
    const vector<Jump>& jumps,
    const unordered_map<uint64_t, uint32_t>& node_to_cycle_map
);

/**
 * @internal
 * @brief Checks whether the edges have some cycle
 * 
 * @param edges as an adjecency list
 * 
 * @return true If a cycle exists
 * @return false If no cycle was found
 */
bool has_cycle(const unordered_map<uint32_t, vector<uint32_t>>& edges);

/**
 * @internal
 * @brief Returns the edge with the minimum weight and its confidence that it will resolve the cycles
 * 
 * @todo outdated doc
 * 
 * The edge with the smallest weight will be returned.
 * The confidence is just the relative weight amount compared to the total weight
 * 
 * @post 0.0 <= confidence <= 1.0
 * 
 * @param edges_with_weights Used for the edge_weight
 * 
 * @return Best edge (tuple<uint32_t, uint32_t>) and confidence (float)
 */
pair<tuple<uint32_t, uint32_t>, int> resolve_cycles_greedy_best_edge(
    const unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash>& edges_with_weights
);

/**
 * @internal
 * @brief Resolves any cycles caused when viewing constraints as edges
 * 
 * @todo outdated doc
 * 
 * @pre 0.0 < confidence <= 1.0
 * 
 * @post Removes constraints that lead to cycles
 * @post confidence will be reduced according to the resolving confidence
 * @post 0.0 < confidence <= 1.0
 * 
 * @param constraints Viewed as edges of a graph
 */
void resolve_cycles_greedy(
    vector<tuple<uint32_t, uint32_t>>& constraints,
    unordered_map<uint32_t, int>& heuristic_node_values
);

/**
 * @internal
 * @brief Recursive algorithm to solve topological sort
 * 
 * 1. Chooses the first possible start node from the list
 * 2. Removes any influence of that start_node
 * 3. Finds new possible start nodes
 * 4. Repeat 1,2,3 until no possible start node is left
 * 
 * @post possible_choices is a number estimating the amount of choices the topological
 * sort could have taken during the algorithm
 * @post total_order will contain the end result
 * 
 * @param possible_start_nodes Choices for the next node
 * @param node_affection_to_start Used as a heuristic for what start node to choose next
 * @param edges Used to increase the start nodes list for the next rec call
 * @param heuristic_node_values The other part of the heuristic to help choose the next start node
 * @param total_order The resulting order
 */
void apply_topological_sort(
    vector<uint32_t>& possible_start_nodes,
    const unordered_map<uint32_t, int>& node_affection_to_start,
    unordered_map<uint32_t, int>& heuristic_node_values,
    unordered_map<tuple<uint32_t, uint32_t>, int, TupleHash>& edges,
    vector<uint32_t>& total_order
);

/**
 * @internal
 * @brief Finds an order that complies with every constraint
 * 
 * Constructs edges from the constraints and uses topological sort to get an answer
 * 
 * @pre The constraints cannot have any cycles (ignoring constraints that contain uint32_t max)
 * 
 * @post confidence is in [0.0, 1.0]
 * @post The result is a permutation of nodes
 * 
 * @param constraints Used to get an order
 * @param heuristic_node_values Used for the recursive toposort call
 * @param nodes Result will contain all of these elements
 * @param confidence Shows an estimate of how deterministic the result is
 * 
 * @return Order layed out by constraints (vector<uint32_t>)
 */
vector<uint32_t> solve_constraints_with_topological_sort(
    const vector<tuple<uint32_t, uint32_t>>& constraints,
    unordered_map<uint32_t, int>& heuristic_node_values,
    const vector<uint32_t>& nodes
);

/**
 * @brief Reconstructs the order of cycles using the jumps
 * 
 * Takes the jumps to derive constraints, which are used to find a total order
 * 
 * @todo Don't use the get_all_cycle_indices, instead rely on cycles.size()
 * 
 * @post confidence is in [0.0, 1.0]
 * @post The result is a permutation of nodes from the graph
 * 
 * @param graph 
 * @param jumps 
 * @param cycles 
 * 
 * @return vector<uint32_t> 
 */
vector<uint32_t> order_cycles(
    const Graph& graph,
    const vector<Jump>& jumps,
    const vector<vector<uint64_t>>& cycles
);

/**
 * @brief Uses the permutation of cycle_order to order the cycles
 * 
 * @note Any invalid cycle index in the cycle_order will be ignored
 * 
 * @param cycle_order Represents the order of the cycles
 * @param cycles The cycles that need to be ordered
 * 
 * @return The ordered cycles (vector<uint64_t>)
 */
vector<vector<uint64_t>> get_ordered_cycles(
    const vector<uint32_t>& cycle_order,
    const vector<vector<uint64_t>>& cycles
);

#endif
