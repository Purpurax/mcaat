#include "isolate_protospacers.h"
#include <vector>
#include <sstream>

IsolateProtospacers::IsolateProtospacers(const SDBG& g, const std::string& cyclesFilePath) : graph(g) {
    this->cycleNodes= ReadCycles(cyclesFilePath);
    protospacerNodes = IsolateProtospacersMethod();
}

IsolateProtospacers::IsolateProtospacers(const SDBG& g, const std::map<uint64_t, std::vector<std::vector<uint64_t>>>& repeatToSpacerNodes) : graph(g) {
    for (const auto& pair : repeatToSpacerNodes) {
        uint64_t groupId = pair.first;
        const auto& vecVec = pair.second;
        for (const auto& vec : vecVec) {
            if (!vec.empty()) {
                uint64_t cycleId = vec[0];
                std::set<uint64_t> nodes(vec.begin(), vec.end());
                cycleNodes[cycleId] = nodes;
                cycleToGroup[cycleId] = groupId;
            }
        }
    }
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
    std::map<uint64_t, std::set<uint64_t>> incomingOutersMap, outgoingOutersMap;

    for (const auto& pair : cycleNodes) {
        uint64_t k = pair.first;
        const std::set<uint64_t>& cycleNodesSet = pair.second;

        // Find outer incomings: nodes outside that point into the cycle
        std::set<uint64_t> incomingOuters;
        for (uint64_t cycleNode : cycleNodesSet) {
            int indegree = graph.EdgeIndegree(cycleNode);
            if (indegree > 0) {
                std::vector<uint64_t> incomings(indegree);
                graph.IncomingEdges(cycleNode, incomings.data());
                for (uint64_t neighbor : incomings) {
                    if (cycleNodesSet.find(neighbor) == cycleNodesSet.end() && cycleNodes.find(neighbor) == cycleNodes.end()) {
                        incomingOuters.insert(neighbor);
                    }
                }
            }
        }

        // Find outer outgoings: nodes outside that the cycle points to
        std::set<uint64_t> outgoingOuters;
        for (uint64_t cycleNode : cycleNodesSet) {
            int outdegree = graph.EdgeOutdegree(cycleNode);
            if (outdegree > 0) {
                std::vector<uint64_t> outgoings(outdegree);
                graph.OutgoingEdges(cycleNode, outgoings.data());
                for (uint64_t neighbor : outgoings) {
                    if (cycleNodesSet.find(neighbor) == cycleNodesSet.end() && cycleNodes.find(neighbor) == cycleNodes.end()) {
                        outgoingOuters.insert(neighbor);
                    }
                }
            }
        }

        // Only add if both exist
        if (!incomingOuters.empty() && !outgoingOuters.empty()) {
            incomingOutersMap[k] = incomingOuters;
            outgoingOutersMap[k] = outgoingOuters;
        }
    }
    //print all the incomingOutersMap and outgoingOutersMap
    std::cout << "Identified " << incomingOutersMap.size() << " cycles with incoming protospacer nodes." << std::endl;
    std::cout << "Identified " << outgoingOutersMap.size() << " cycles with outgoing protospacer nodes." << std::endl;
    for (const auto& [cycleId, nodes] : incomingOutersMap) {
        std::cout << "Cycle " << cycleId << " incoming protospacer nodes: ";
        for (uint64_t node : nodes) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
    for (const auto& [cycleId, nodes] : outgoingOutersMap) {
        std::cout << "Cycle " << cycleId << " outgoing protospacer nodes: ";
        for (uint64_t node : nodes) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    return {incomingOutersMap, outgoingOutersMap};
}

std::pair<std::map<uint64_t, std::set<uint64_t>>, std::map<uint64_t, std::set<uint64_t>>> IsolateProtospacers::FilterSharedGroups(const std::map<uint64_t, std::set<uint64_t>>& inGroup, const std::map<uint64_t, std::set<uint64_t>>& outGroup) {
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
    return {possibleInProtospacerNodes, possibleOutProtospacerNodes};
}

void IsolateProtospacers::DepthLimitedSearch(uint64_t currentNode, int depth, std::vector<uint64_t>& path, std::set<uint64_t>& visited, const std::set<uint64_t>& outNodes, const std::set<uint64_t>& cycleNodeSet, int maxDepth, int minDepth, std::vector<std::vector<uint64_t>>& successfulPaths) {
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
                if (visited.find(neighbor) == visited.end() && (cycleNodeSet.find(neighbor) != cycleNodeSet.end() || outNodes.find(neighbor) != outNodes.end()) && cycleNodes.find(neighbor) == cycleNodes.end()) {
                    DepthLimitedSearch(neighbor, depth + 1, path, visited, outNodes, cycleNodeSet, maxDepth, minDepth, successfulPaths);
                }
            }
        }
    }

    path.pop_back();
    visited.erase(currentNode);
}

// do a dls from each incoming protospacer node
std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> IsolateProtospacers::DepthLimitedPathsFromInToOut(
    const std::map<uint64_t, std::set<uint64_t>>& inGroup,
    const std::map<uint64_t, std::set<uint64_t>>& outGroup,
    int maxDepth,
    int minDepth
) {
    
    auto [possibleInProtospacerNodes, possibleOutProtospacerNodes] = FilterSharedGroups(inGroup, outGroup);
    //print first values of the maps
    std::cout << "After filtering, we have " << possibleOutProtospacerNodes.size() << " cycles with both incoming and outgoing protospacer nodes." << std::endl;
    
    
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> groupedPaths;
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

        std::vector<std::vector<uint64_t>> cyclePaths;
        for (uint64_t startNode : inNodes) {
            std::vector<uint64_t> path;
            std::set<uint64_t> visited;
            int cycleMaxDepth = cycleNodeSet.size();
            DepthLimitedSearch(startNode, 0, path, visited, outNodes, cycleNodeSet, cycleMaxDepth, minDepth, cyclePaths);
        }

        // Filter out subpaths for this cycle
        // Sort by length descending
        std::sort(cyclePaths.begin(), cyclePaths.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
            return a.size() > b.size();
        });
        std::vector<std::vector<uint64_t>> filtered;
        for (const auto& p : cyclePaths) {
            bool isSubpath = false;
            for (const auto& longer : filtered) {
                if (longer.size() > p.size() && std::search(longer.begin(), longer.end(), p.begin(), p.end()) != longer.end()) {
                    isSubpath = true;
                    break;
                }
            }
            if (!isSubpath) {
                filtered.push_back(p);
            }
        }

        // Select node-disjoint paths (greedy: longest first)
        std::vector<std::vector<uint64_t>> disjointPaths;
        std::set<uint64_t> usedNodes;
        for (const auto& p : filtered) {
            bool canAdd = true;
            for (uint64_t node : p) {
                if (usedNodes.count(node)) {
                    canAdd = false;
                    break;
                }
            }
            if (canAdd) {
                disjointPaths.push_back(p);
                for (uint64_t node : p) {
                    usedNodes.insert(node);
                }
            }
        }
        filtered = disjointPaths;

        // Trim each path to [1:-1] for protospacers
        std::vector<std::vector<uint64_t>> trimmedPaths;
        for (const auto& p : filtered) {
            if (p.size() > 2) {
                std::vector<uint64_t> trimmed(p.begin() + 1, p.end() - 1);
                trimmedPaths.push_back(trimmed);
            }
        }

        uint64_t groupId = cycleToGroup[cycleStartNode];
        groupedPaths[groupId][cycleStartNode] = trimmedPaths;
    }

    // Global subpath filtering across all paths
    std::vector<std::vector<uint64_t>> allPaths;
    std::map<std::vector<uint64_t>, std::set<uint64_t>> pathToCycles;
    for (const auto& group : groupedPaths) {
        for (const auto& cycle : group.second) {
            uint64_t cycleId = cycle.first;
            for (const auto& path : cycle.second) {
                allPaths.push_back(path);
                pathToCycles[path].insert(cycleId);
            }
        }
    }

    // Sort by length descending
    std::sort(allPaths.begin(), allPaths.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
        return a.size() > b.size();
    });

    std::vector<std::vector<uint64_t>> globalFiltered;
    for (const auto& p : allPaths) {
        bool isSubpath = false;
        for (const auto& longer : globalFiltered) {
            if (longer.size() > p.size() && std::search(longer.begin(), longer.end(), p.begin(), p.end()) != longer.end()) {
                isSubpath = true;
                break;
            }
        }
        if (!isSubpath) {
            globalFiltered.push_back(p);
        }
    }

    // Unique the global filtered
    std::set<std::vector<uint64_t>> uniqueGlobal(globalFiltered.begin(), globalFiltered.end());
    globalFiltered.assign(uniqueGlobal.begin(), uniqueGlobal.end());

    // Assign each unique path to the cycle with the smallest ID that had it
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> newGroupedPaths;
    for (const auto& path : globalFiltered) {
        auto it = pathToCycles.find(path);
        if (it != pathToCycles.end()) {
            uint64_t minCycleId = *std::min_element(it->second.begin(), it->second.end());
            uint64_t groupId = cycleToGroup[minCycleId];
            newGroupedPaths[groupId][minCycleId].push_back(path);
        }
    }

    return newGroupedPaths;
}

void IsolateProtospacers::WritePathsToFile(const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& paths, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const auto& group_pair : paths) {
        uint64_t group_id = group_pair.first;
        const auto& cycle_map = group_pair.second;
        file << "Group " << group_id << ":\n";
        for (const auto& cycle_pair : cycle_map) {
            uint64_t cycle_id = cycle_pair.first;
            const auto& path_list = cycle_pair.second;
            if (!path_list.empty()) {
                file << "  Cycle " << cycle_id << ":\n";
                for (const auto& path : path_list) {
                    file << "    [";
                    for (size_t i = 0; i < path.size(); ++i) {
                        file << path[i];
                        if (i < path.size() - 1) file << " ";
                    }
                    file << "]\n";
                }
            }
        }
        file << "\n";
    }

    file.close();
}