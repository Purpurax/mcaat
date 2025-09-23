#include "isolate_protospacers.h"
#include <vector>
#include <sstream>

IsolateProtospacers::IsolateProtospacers(const SDBG& g, const std::string& cyclesFilePath) : graph(g) {
    this->cycleNodes= ReadCycles(cyclesFilePath);
    protospacerNodes = IsolateProtospacersMethod();
}

std::map<uint64_t, std::set<uint64_t>> IsolateProtospacers::ReadCycles(const std::string& filePath) {
    std::map<uint64_t, std::set<uint64_t>> cyclesMap;
    std::set<uint64_t> allNodes;  // Big set for all nodes, as per user request

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return cyclesMap;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        uint64_t firstNode;
        ss >> firstNode;
        std::set<uint64_t> nodesSet;
        uint64_t node;
        while (ss >> node) {
            if (node != firstNode) {  // Ensure first node not in values
                nodesSet.insert(node);
            }
            allNodes.insert(node);  // Add to big set
        }
        allNodes.insert(firstNode);  // Add first node to big set
        if (!nodesSet.empty()) {
            cyclesMap[firstNode] = nodesSet;
        }
    }

    file.close();
    return cyclesMap;
}

IN_OUT_PAIR_MAP_SET IsolateProtospacers::IsolateProtospacersMethod() {
    IN_OUT_PAIR_MAP_SET possibleProtospacerNodes;
    std::map<uint64_t, std::set<uint64_t>> possibleOutProtospacerNodes;
    std::map<uint64_t, std::set<uint64_t>> possibleInProtospacerNodes;
    // create a new map that saves only the groups that have nodes that occur only once in that group using Counts like mechanism in python
    std::map<uint64_t, std::set<uint64_t>> cycleNodesUnique;
    for (const auto& pair : cycleNodes) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;
        std::map<uint64_t, int> nodeCount;
        for (uint64_t v : values) {
            nodeCount[v]++;
        }
        std::set<uint64_t> uniqueNodes;
        for (const auto& countPair : nodeCount) {
            if (countPair.second == 1) {
                uniqueNodes.insert(countPair.first);
            }
        }
        if (!uniqueNodes.empty()) {
            cycleNodesUnique[k] = uniqueNodes;
        }
    }


    for (const auto& pair : cycleNodesUnique) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;

        for (uint64_t v : values) {
            // Check for outgoing external neighbor
            int outdegree = graph.EdgeOutdegree(v);
            if (outdegree > 0) {
                std::vector<uint64_t> outgoings(outdegree);
                graph.OutgoingEdges(v, outgoings.data());
                for (uint64_t neighbor : outgoings) {
                    if (values.find(neighbor) == values.end()) {
                        possibleOutProtospacerNodes[k].insert(v);
                        break;
                    }
                }
            }
            // Check for incoming external neighbor
            int indegree = graph.EdgeIndegree(v);
            if (indegree > 0) {
            std::vector<uint64_t> incomings(indegree);
            graph.IncomingEdges(v, incomings.data());
            for (uint64_t neighbor : incomings) {
                if (values.find(neighbor) == values.end()) {
                    possibleInProtospacerNodes[k].insert(v);
                    break;
                }
            }
        }
    }
    }

return {possibleOutProtospacerNodes, possibleInProtospacerNodes};

}

void IsolateProtospacers::PrintProtospacerNodesToConsole(const std::map<uint64_t, std::set<uint64_t>>& possibleProtospacerNodes, const std::string& type_edge) {
    std::cout << "Identified " << possibleProtospacerNodes.size() << " potential protospacer nodes." << std::endl;

    for (const auto& pair : possibleProtospacerNodes) {
        std::cout << "Cycle starting at node " << pair.first << " has " << pair.second.size() << " possible protospacer nodes." << std::endl;
        for (uint64_t node : pair.second) {
            std::cout << "  " << node << " : " << graph.EdgeMultiplicity(node) <<"; [";
            //list all the neighbors of the protospacer node
            if(type_edge == "outgoing" ){
                int outdegree = graph.EdgeOutdegree(node);
                if (outdegree > 0) {
                    std::vector<uint64_t> outgoings(outdegree);
                    graph.OutgoingEdges(node, outgoings.data());
                    std::cout << "  Outgoing neighbors: ";
                    for (uint64_t neighbor : outgoings) {
                        std::cout << neighbor << " : "<< graph.EdgeMultiplicity(neighbor) << ",";
                    }
                    std::cout << "]" << std::endl;
                }
                else{
                    //incomings
                    int indegree = graph.EdgeIndegree(node);
                    if (indegree > 0) {
                        std::vector<uint64_t> incomings(indegree);
                        graph.IncomingEdges(node, incomings.data());
                        std::cout << "  Incoming neighbors: ";
                        for (uint64_t neighbor : incomings) {
                            std::cout << neighbor << " : "<< graph.EdgeMultiplicity(neighbor) << ",";
                        }
                        std::cout << "]" << std::endl;
                    }
                }
            }
        }
    }
}

// do a dls from each incoming protospacer node
// and see if in depth of maximum depth we can reach the outgoing protospacer nodes
std::vector<std::vector<uint64_t>> IsolateProtospacers::DepthLimitedPathsFromInToOut(
    const std::map<uint64_t, std::set<uint64_t>>& inGroup,
    const std::map<uint64_t, std::set<uint64_t>>& outGroup,
    int maxDepth,
    int minDepth
) {
    
    //keep only those that share the same group
    std::map<uint64_t, std::set<uint64_t>> possibleInProtospacerNodes, possibleOutProtospacerNodes;
    for (const auto& pair : inGroup) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;
        if (outGroup.find(k) != outGroup.end()) {
            possibleInProtospacerNodes[k] = values;
        }
    }
    for (const auto& pair : outGroup) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& values = pair.second;
        if (inGroup.find(k) != inGroup.end()) {
            possibleOutProtospacerNodes[k] = values;
        }
    }
    //print first values of the maps
    std::cout << "After filtering, we have " << possibleOutProtospacerNodes.size() << " cycles with both incoming and outgoing protospacer nodes." << std::endl;
    
    
    std::vector<std::vector<uint64_t>> successfulPaths;
    std::function<void(uint64_t, int, std::vector<uint64_t>&, std::set<uint64_t>&, const std::set<uint64_t>&, const std::set<uint64_t>&)> dls = [&](uint64_t currentNode, int depth, std::vector<uint64_t>& path, std::set<uint64_t>& visited, const std::set<uint64_t>& outNodes, const std::set<uint64_t>& cycleNodeSet) {
        if (depth > maxDepth) return;
        visited.insert(currentNode);
        path.push_back(currentNode);

        if (outNodes.find(currentNode) != outNodes.end() && depth >= minDepth) {
            successfulPaths.push_back(path);
        } else {
            int outdegree = graph.EdgeOutdegree(currentNode);
            if (outdegree > 0) {
                std::vector<uint64_t> outgoings(outdegree);
                graph.OutgoingEdges(currentNode, outgoings.data());
                for (uint64_t neighbor : outgoings) {
                    if (visited.find(neighbor) == visited.end() && cycleNodeSet.find(neighbor) != cycleNodeSet.end()) {
                        dls(neighbor, depth + 1, path, visited, outNodes, cycleNodeSet);
                    }
                }
            }
        }

        path.pop_back();
        visited.erase(currentNode);
    };


    for (const auto& inPair : possibleInProtospacerNodes) {
        uint64_t cycleStartNode = inPair.first;
        const std::set<uint64_t>& inNodes = inPair.second;

        auto outIt = possibleOutProtospacerNodes.find(cycleStartNode);
        if (outIt == possibleOutProtospacerNodes.end()) {
            continue; // No corresponding outgoing protospacer nodes
        }
        const std::set<uint64_t>& outNodes = outIt->second;

        // Get the full cycle node set
        auto cycleIt = this->cycleNodes.find(cycleStartNode);
        if (cycleIt == this->cycleNodes.end()) {
            continue; // No cycle nodes found
        }
        const std::set<uint64_t>& cycleNodeSet = cycleIt->second;

        for (uint64_t startNode : inNodes) {
            std::vector<uint64_t> path;
            std::set<uint64_t> visited;
            dls(startNode, 0, path, visited, outNodes, cycleNodeSet);
        }
    }
    return successfulPaths;
}