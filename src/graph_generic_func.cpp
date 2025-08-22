#include "graph_generic_func.h"


// generic helper functions for graph operations

void graph_generic_func::_GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set, SDBG& sdbg) {
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !sdbg.IsValidEdge(node)) {
        return;
    }
     uint64_t outgoings[edge_outdegree];
    int flag =sdbg.OutgoingEdges(node, outgoings);
    if(flag!=-1)    
        for (const auto& outgoing : outgoings)
            outgoings_set.insert(outgoing);
    
    
}
void graph_generic_func::_GetIncomings(auto node, unordered_set<uint64_t>& incomings_set, SDBG& sdbg) {
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !sdbg.IsValidEdge(node)) {
        return;
    }
    uint64_t incomings[edge_indegree];
    int flag = sdbg.IncomingEdges(node, incomings);
    if (flag!=-1)
        for (const auto& incoming : incomings)
            incomings_set.insert(incoming);
}