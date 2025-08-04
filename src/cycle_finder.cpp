#include "cycle_finder.h"
#include "filters.h"

// Parallel hashmap for better performance in DLS
#include <parallel_hashmap/phmap.h>

/**
 * @file cycle_finder.cpp
 * @brief Implementation of functions for cycle detection and analysis in a sequence graph.
 * - Checking if any incoming edge of a node is not equal to the node itself.
 * - Performing a background check on a neighbor node to determine if it meets certain criteria.
 * - Getting the outgoing edges of a node that pass the background check.
 * - Retrieving the incoming edges of a node that pass the background check.
 * - Finding cycles in the graph.
 * - Performing a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * - Collecting tips (nodes with no outgoing edges) in the graph.
 * - Performing recursive reduction to invalidate edges.
 * - Invalidating nodes with edge multiplicity of one.
 * - Chunking start nodes based on their multiplicity for parallel processing.
 * - Finding all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */


/**
 * @brief Checks if any incoming edge of a node is not equal to the node itself.
 */
bool CycleFinder::_IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree) {
    uint64_t incomings[edge_indegree];
    this->sdbg.IncomingEdges(node, incomings);
    for (const auto& incoming : incomings)
        if (node==incoming)
            return true;
    return false;
}
/**
 * @brief Performs a background check on a neighbor node to determine if it meets certain criteria.
 */
bool CycleFinder::_BackgroundCheck(uint64_t original_node, size_t repeat_multiplicity, uint64_t neighbor_node) {
    auto neighbor_node_multiplicity = sdbg.EdgeMultiplicity(neighbor_node);
    if(this->visited[neighbor_node]) {
        return false;
    }
    if (repeat_multiplicity / neighbor_node_multiplicity > 500) {
        return false;
    }
    if(original_node == neighbor_node) {
        return false;
    }
    return true;
}


/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity) {
   
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
     uint64_t outgoings[edge_outdegree];
    int flag =sdbg.OutgoingEdges(node, outgoings);
    if(flag!=-1)    
        for (const auto& outgoing : outgoings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, outgoing))
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity) {
  
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag =sdbg.IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            if (this->_BackgroundCheck(node, repeat_multiplicity, incoming))
                incomings_set.insert(incoming);
}
// ## START: HELPER FUNCTIONS FOR DLS ##
/**
 * @brief Gets the outgoing edges of a node that pass the background check.
 */
void CycleFinder::_GetOutgoings(uint64_t node, unordered_set<uint64_t>& outgoings_set) {
   
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t outgoings[edge_outdegree];
    int flag = sdbg.OutgoingEdges(node, outgoings);
    if(flag!=-1)
        for (const auto& outgoing : outgoings)
            //if (sdbg.EdgeMultiplicity(outgoing) > 1)
                outgoings_set.insert(outgoing);
    
    
}
/**
 * @brief Retrieves the incoming edges of a node that pass the background check.
 */
void CycleFinder::_GetIncomings(uint64_t node, unordered_set<uint64_t>& incomings_set) {
  
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !this->sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = sdbg.IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            //if (sdbg.EdgeMultiplicity(incoming) > 1)
                incomings_set.insert(incoming);
}

// ## END: HELPER FUNCTIONS FOR DLS ##

/**
 * @brief Finds the cycles in the graph
 * 
 * @param sdbg The succinct de Bruijn graph.
 * @param length_bound The maximum length of the path to be considered.
 * @param minimal_length The minimum length of the path to be considered.
 * @param genome_name The name of the genome.
 */
CycleFinder::CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name,int threads_count)
    : maximal_length(length_bound), minimal_length(minimal_length), sdbg(sdbg), genome_name(genome_name), cluster_bounds(500), threads_count(threads_count) {
    this->FindApproximateCRISPRArrays();
}

vector<vector<uint64_t>> CycleFinder::FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, 
                                        vector<int> backtrack_lengths) {
    int counter = 0;
    uint64_t current_node = start_node;
    vector<vector<uint64_t>> cycles;
    int steps_counter = 0;

    while (!stack.empty()) {
        steps_counter += 1;
        if (steps_counter > 10000000) {
            break;
        }
        
        unordered_set<uint64_t> neighbors = stack.back();
        bool flag = true;
        for (auto neighbor : neighbors) {
            current_node = neighbor;
            if (current_node == start_node ) {
                backtrack_lengths.back() = 1;
                if (path.size() > this->minimal_length) {
                    cycles.push_back(path);
                    counter += 1;
                    if (counter >= this->cluster_bounds) {
                        cycles.clear();
                        flag=false;
                    }
                }
            } 
            else if (static_cast<int>(path.size()) < lock.try_emplace(neighbor, this->maximal_length).first->second) {
                neighbors.erase(neighbor);
                path.push_back(neighbor);
                backtrack_lengths.push_back(this->maximal_length);
                lock[neighbor] = path.size();
                stack.back().erase(neighbor);
                unordered_set<uint64_t> outgoings;
                this->_GetOutgoings(neighbor, outgoings, sdbg.EdgeMultiplicity(start_node));
                stack.push_back(outgoings);
                flag = false;
                break;
            }
        }
        if (flag) {
            stack.pop_back();
            uint64_t v = path.back();
            path.pop_back();
            int backtrack_length = backtrack_lengths.back();
            backtrack_lengths.pop_back();

            if (!backtrack_lengths.empty()) {
                backtrack_lengths.back() = min(backtrack_lengths.back(), backtrack_length);
            }
            if (backtrack_length < this->maximal_length) {
                vector<pair<int, int>> relax_stack;
                relax_stack.push_back(make_pair(backtrack_length, v));

                unordered_set<uint64_t> path_set(path.begin(), path.end());

                while (!relax_stack.empty()) {
                    int bl = relax_stack.back().first;
                    int u = relax_stack.back().second;
                    relax_stack.pop_back();
                    if (lock.try_emplace(u, this->maximal_length).first->second < this->maximal_length - bl + 1) {
                        lock[u] = this->maximal_length - bl + 1;
                        unordered_set<uint64_t> incomings;
                        this->_GetIncomings(u, incomings, sdbg.EdgeMultiplicity(start_node));
                        for (auto w : incomings)
                            if (path_set.find(w) == path_set.end())
                                relax_stack.push_back(make_pair(bl + 1, w));
                    }
                }
            }
        }
    }
    
    
    if(cycles.empty()) return {};
    
    #pragma omp critical
    {
        for (const auto& cycle : cycles) 
            for (const auto& node : cycle) 
                this->visited[node] = true;
    }
    
    return cycles;
}


/**
 * @brief Utility function to initialize and start the cycle finding process from a given node.
 */
vector<vector<uint64_t>> CycleFinder::FindCycleUtil(uint64_t start_node) {
    vector<uint64_t> path;
    map<uint64_t, int> lock;
    vector<unordered_set<uint64_t>> stack;
    vector<int> backtrack_lengths;
    path.push_back(start_node);
    lock[start_node] = 0;
    unordered_set<uint64_t> outgoings;
    this->_GetOutgoings(start_node, outgoings, sdbg.EdgeMultiplicity(start_node));
    stack.push_back(outgoings);
    backtrack_lengths.push_back(maximal_length);
    return FindCycle(start_node, path, lock, stack, backtrack_lengths);
}
/**
 * @brief Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.
 * Optimized using megahit's graph traversal patterns.
 */
bool CycleFinder::DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth) {
    // Megahit-style memory pool: Dynamic thread-local pools for reuse (no heuristic sizing)
    struct StackEntry {
        uint64_t node;
        int depth;
    };
    
    static thread_local std::vector<StackEntry> dls_stack_pool;
    static thread_local phmap::flat_hash_set<uint64_t> dls_visited_pool;
    
    // Clear but keep capacity (megahit memory reuse pattern) - no fixed reserve
    dls_stack_pool.clear();
    dls_visited_pool.clear();
    
    // Reserve reasonable initial capacity to avoid reallocations
    if (dls_stack_pool.capacity() == 0) {
        dls_stack_pool.reserve(64);  // Small initial reserve
    }

    // Use pooled structures
    auto& dls_stack = dls_stack_pool;
    auto& dls_visited = dls_visited_pool;

    // Cache values for faster comparison
    const uint64_t target_node = target;
    const uint64_t start_node = start;

    dls_stack.push_back({start_node, 0});
    reached_depth = 0;

    while (!dls_stack.empty()) {
        StackEntry current = dls_stack.back();
        dls_stack.pop_back();
        uint64_t v = current.node;
        int depth = current.depth;

        // Update reached depth
        reached_depth = depth;

        // Get neighbors using SDBG API directly (megahit's pattern) - moved up for better branch prediction
        // Use megahit's efficient zero-degree check for early termination
        if (__builtin_expect(sdbg.EdgeOutdegreeZero(v), 0)) {
            continue;
        }
        
        int outdegree = sdbg.EdgeOutdegree(v);

        // Use fixed-size array for neighbors (megahit's pattern)
        uint64_t neighbors[MAX_EDGE_COUNT];
        int flag = sdbg.OutgoingEdges(v, neighbors);

        if (__builtin_expect(flag == -1, 0)) {
            continue;
        }

        // Prefetch neighbors for better cache performance
        __builtin_prefetch(&neighbors[0], 0, 1);

        // Exceeded depth limit - check after we know we have neighbors
        if (__builtin_expect(depth >= limit, 0)) {
            continue;
        }

        // Process all neighbors to maintain correctness (removed faulty simple path optimization)
        // Process neighbors in forward order (SIMD-friendly pattern from megahit)
        // Unroll loop for small outdegrees to reduce overhead
        // Unrolled loop for common case (de Bruijn graph max degree = 4)
        for (int i = 0; i < outdegree; ++i) {
            uint64_t neighbor = neighbors[i];
            auto visited_it = dls_visited.find(neighbor);
            bool not_visited = (visited_it == dls_visited.end());
            bool is_start_revisit = (neighbor == start_node && depth > 0);
                
            if (__builtin_expect(not_visited || is_start_revisit, 1)) {
                dls_visited.insert(neighbor);
                dls_stack.push_back({neighbor, depth + 1});
            }
        }
        

        // Megahit-style aggressive early cycle detection
        if (__builtin_expect(v == target_node && depth > 1, 0)) {
            return true;  // Found cycle - exit immediately
        }
    }

    return false;
}

vector<uint64_t> CycleFinder::CollectTips() {
    
    unordered_set<uint64_t> tips;

    #pragma omp parallel for
    for (uint64_t node = 0; node < this->sdbg.size(); node++) 
        if(this->sdbg.EdgeOutdegree(node) == 0) {
            #pragma omp critical
            tips.insert(node);
        }
    return vector<uint64_t>(tips.begin(), tips.end());
}

void CycleFinder::RecursiveReduction(uint64_t tip) {
    if (this->sdbg.EdgeOutdegree(tip)> 0) 
        return;
    unordered_set<uint64_t> parents;
    this->_GetIncomings(tip, parents);
    this->sdbg.SetInvalidEdge(tip);
    for (uint64_t parent : parents) 
        if(this->sdbg.IsValidEdge(parent)) 
            this->RecursiveReduction(parent);
        else 
            continue;
    return;
}
void CycleFinder::InvalidateMultiplicityOneNodes() {
    #pragma omp parallel for
    for (uint64_t node = 0; node < this->sdbg.size(); node++) {
        if (this->sdbg.EdgeMultiplicity(node) == 1) {
            this->sdbg.SetInvalidEdge(node);
        }
    }
}
/**
 * @brief Chunks the start nodes based on their multiplicity for parallel processing.
 */
size_t CycleFinder::ChunkStartNodes(map<int, vector<uint64_t>, greater<int>>& start_nodes_chunked) {
    uint64_t loaded = 0;
    const int chunk_size = 20000;
    //this->InvalidateMultiplicityOneNodes();
    #pragma omp parallel num_threads(this->threads_count)
    {
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t node = 0; node < this->sdbg.size(); node++) {
                size_t edge_indegree = this->sdbg.EdgeIndegree(node);
                // size_t edge_outdegree = this->sdbg.EdgeOutdegree(node); // unused
                loaded+=1; 
                if(loaded % 10000000 == 0) std::cout << "Loaded " << loaded << " nodes\n";
                if (edge_indegree >= 2 && this->sdbg.EdgeMultiplicity(node) > 20)
                {
                    
                    if(this->_IncomingNotEqualToCurrentNode(node,edge_indegree)) continue;
                    int reached_depth = 0;
                    
                    bool dls = this->DepthLevelSearch(node, node, this->maximal_length, reached_depth);
                    if(!dls) continue; //|-> the last version!
            
                    double log2_mult = ceil(log2(double(this->sdbg.EdgeMultiplicity(node))));
                    #pragma omp critical
                    start_nodes_chunked[log2_mult].push_back(node);
            }
        }
    }
   //writeStartNodesToFile(start_nodes_chunked, "start_nodes.txt");
    size_t sum_of_all_quantities_in_all_chunks = 0;
    for (const auto& [key, value] : start_nodes_chunked) {
        std::cout << "log2_mult: " << key << " Number of nodes: " << value.size() << endl;
        sum_of_all_quantities_in_all_chunks += value.size();
    }
    return sum_of_all_quantities_in_all_chunks;
}


/**
 * @brief Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.
 */
int CycleFinder::FindApproximateCRISPRArrays()
 {
    /*
    vector<uint64_t> tips = this->CollectTips();
    std::cout<<"BEFORE number of nodes: " << this->sdbg.size() << endl;
    std::cout<< "Number of tips: " << tips.size() << endl;

    for (uint64_t tip : tips) {
        this->RecursiveReduction(tip);
    } 
    int invalid_edges = 0;
    #pragma omp parallel for reduction(+:invalid_edges)
    for (uint64_t node = 0; node < this->sdbg.size(); node++) {
        if (node==-1) {
            invalid_edges += 1;
        }
    }
    */
    
    // struct mallinfo mem_info = mallinfo(); // deprecated
    // size_t graph_mem_info = mem_info.uordblks; // unused
    int cumulative = 0;
    printf("Number of nodes in a graph: %lu\n", this->sdbg.size());
    string mode = "fastq";
    this->visited.resize(this->sdbg.size(), false);
    std::cout << "Starting the cycle search:\n\n";
        std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> all_cycles;

    map<int, vector<uint64_t>, greater<int>> start_nodes_chunked;
    size_t start_nodes_amount=this->ChunkStartNodes(start_nodes_chunked);
    std::cout << "Number of start_nodes: " << start_nodes_amount << endl;
    size_t counter = 0;
    size_t n_th_counter = 0;
    for (auto nodes_iterator = start_nodes_chunked.begin(); nodes_iterator != start_nodes_chunked.end(); nodes_iterator++) {
        auto thread_count = this->threads_count;
        if (static_cast<int>(nodes_iterator->second.size()) < thread_count)
            thread_count = nodes_iterator->second.size();
        // #pragma omp parallel for num_threads(thread_count) reduction(+:cumulative) shared(nodes_iterator, sdbg, visited)
        for (uint64_t start_node_index = 0; start_node_index < nodes_iterator->second.size(); start_node_index++) {
            uint64_t start_node = nodes_iterator->second[start_node_index];
            if (this->visited[start_node]) continue;
            vector<vector<uint64_t>> cycles = this->FindCycleUtil(start_node);
            cumulative += cycles.size();
            this->results[start_node] = cycles;
            ++counter;
            if (counter % 100000 == 0) {
                n_th_counter += 1;
                std::cout << n_th_counter << " 100k nodes\n";
            }
        }
        malloc_trim(0);
        std::cout << "processed " << nodes_iterator->second.size() << " with ";
        printf("Cycles: %d\n", cumulative);
    }
    std::cout << "Number of cycles: " << cumulative << endl;
    std::cout << "Number of nodes in results: " << this->results.size() << endl;

    std::cout << endl;
    std::cout << "Number of Cycles: " << cumulative << endl;
    return cumulative;
}

CycleFinder::~CycleFinder() {}