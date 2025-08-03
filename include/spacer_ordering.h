#include <sdbg/sdbg.h>
#include <unordered_set>
#include <stack>
#include <algorithm>
#include "jumps.h"

// #include "definitions.h"
// #include "sorting/kmer_counter.h"
// #include "sorting/read_to_sdbg.h"
// #include "sorting/seq_to_sdbg.h"
// #include "utils/options_description.h"
// #include "utils/utils.h"

struct Graph {
    std::unordered_map<uint64_t, std::vector<uint64_t>> adjacency_list;
    std::unordered_set<uint64_t> nodes;

    void add_edge(uint64_t from, uint64_t to) {
        adjacency_list[from].push_back(to);
        nodes.insert(from);
        nodes.insert(to);
    }
};

std::vector<Graph> divide_graph_into_subgraphs(SDBG& sdbg);
std::vector<Graph> get_crispr_regions(
    SDBG& sdbg,
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles
);

std::vector<Jump> get_relevant_jumps(const Graph& graph, std::vector<Jump>& jumps);
std::vector<std::vector<uint64_t>> get_relevant_cycles(
    const Graph& graph,
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& all_cycles_map
);

std::vector<std::vector<int32_t>> order_cycles(
    const Graph& graph,
    std::vector<Jump>& jumps,
    std::vector<std::vector<uint64_t>>& cycles
);
