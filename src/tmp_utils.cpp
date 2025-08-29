#include "tmp_utils.h"

void trim_string(string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
}

pair<string, optional<string>> get_fastq_files_from_settings(
    const Settings& settings
) {
    const bool two_fastq_files_provided = settings.input_files.find(" ") != string::npos;
    if (two_fastq_files_provided) {
        const size_t space_pos = settings.input_files.find(" ");
        string input_file1 = settings.input_files.substr(0, space_pos);
        string input_file2 = settings.input_files.substr(space_pos + 1);

        trim_string(input_file1);
        trim_string(input_file2);

        return std::make_pair(input_file1, optional<string>(input_file2));
    } else {
        return std::make_pair(settings.input_files, std::nullopt);
    }
}

size_t get_cycle_count(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    size_t count = 0;
    for (const auto& pair : cycles_map) {
        count += pair.second.size();
    }
    return count;
}

void print_sdbg_graph_to_dot_file_convert(const string& lib_file_path) {
    SDBG sdbg;
    sdbg.LoadFromFile(lib_file_path.c_str());


    std::cout << "--------------DOT-FILE------------------" << std::endl;
    std::cout << "digraph TinyGraph {\n";

    for (size_t i = 0; i < sdbg.size(); ++i) {
        std::cout << "    " << i << ";\n";
    }

    for (size_t node = 0; node < sdbg.size(); ++node) {
        int outdegree = sdbg.EdgeOutdegree(node);
        if (outdegree > 0) {
            uint64_t* outgoings = new uint64_t[outdegree];
            if (sdbg.OutgoingEdges(node, outgoings) != -1) {
                for (int i = 0; i < outdegree; ++i) {
                    const uint64_t neighbor = outgoings[i];
                    if (!sdbg.IsValidEdge(neighbor)) continue;

                    std::cout << "    " << node << " -> " << neighbor << ";\n";
                }
            }
            delete[] outgoings;
        }
    }

    std::cout << "}\n";
    std::cout << "--------------DOT-FILE-END--------------" << std::endl;
}


void rotateLeft(std::vector<uint64_t>& arr, int k) {
    int n = arr.size();
    if (n == 0) return; // Handle empty array case
    
    k = k % n; // In case k is greater than the size of the array
    if (k == 0) return; // No rotation needed if k is 0 or a multiple of n
    

    // Step 1: Reverse the first k elements
    std::reverse(arr.begin(), arr.begin() + k);
    
    // Step 2: Reverse the remaining n-k elements
    std::reverse(arr.begin() + k, arr.end());
    
    // Step 3: Reverse the entire array
    std::reverse(arr.begin(), arr.end());
}

std::vector<uint64_t> FindRepeatNodePaths(
    SDBG& sdbg,
    const vector<uint64_t>& repeat_nodes,
    const vector<vector<uint64_t>>& ordered_cycles
) {
    uint64_t start;
    uint64_t end;
    vector<uint64_t> all_the_neighbors;
    
    //LIST
    for(const auto& node : repeat_nodes) {
        uint64_t outgoings[4]; 
        int num_outgoings = sdbg.OutgoingEdges(node, outgoings);
        
        for(int i = 0; i < num_outgoings; i++) {
            all_the_neighbors.push_back(outgoings[i]);
        }
    }

    for(const auto& node : repeat_nodes) {
        if(std::find(all_the_neighbors.begin(), all_the_neighbors.end(), node) == all_the_neighbors.end()) {
            start = node;
        }
    }
        
    int maxSize = 0;
    vector<vector<uint64_t>> cycles_per_group = ordered_cycles;
    std::vector<uint64_t> arr;
    auto it = std::find(arr.begin(), arr.end(), start);
    int position_to_rotate = std::distance(arr.begin(), it);

    // Loop through the list of vectors and find the one with the max number of elements
    for (int i = 0; i < cycles_per_group.size(); i++) {
        rotateLeft(cycles_per_group[i], position_to_rotate);
        if (cycles_per_group[i].size() > maxSize) {
            maxSize = cycles_per_group[i].size();
            arr = cycles_per_group[i];  // Update the vector with the maximum size
        }
    }
    
    rotateLeft(arr, position_to_rotate);

    arr.resize(repeat_nodes.size());

    return arr;
}

pair<vector<uint64_t>, vector<vector<uint64_t>>> FindCRISPRArrayNodes(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& ordered_cycles
) {
    int threshold = static_cast<int>(ordered_cycles.size());
    
    std::unordered_map<uint64_t, int> element_count;
    for (const auto& cycle : ordered_cycles) {
        for (const auto& element : cycle) {
            element_count[element]++;
        }
    }

    std::vector<uint64_t> repeat_nodes;
    for (const auto& [element, count] : element_count) {
        if (count >= threshold) {
            repeat_nodes.push_back(element);
        }
    }

    if (repeat_nodes.size() >= 27) {
        return {{}, {}};
    }
    
    std::vector<std::vector<uint64_t>> spacer_nodes;
    
    repeat_nodes = FindRepeatNodePaths(sdbg, repeat_nodes, ordered_cycles);
    
    for (const auto& cycle : ordered_cycles) {
        // IMPORTANT: DO NOT ADD UPPER BOUNDARY
        if (cycle.size() - repeat_nodes.size() >= 23) {
            std::vector<uint64_t> spacers(cycle.begin() + repeat_nodes.size(), cycle.end());
            
            spacer_nodes.push_back(spacers);
        }
        // IMPORTANT: DO NOT ADD UPPER BOUNDARY
    }

    if(repeat_nodes.size() == 0 || spacer_nodes.size() < 3) {
        return {{}, {}};
    }

    return {repeat_nodes, spacer_nodes};
}

string FetchNodeLabel(SDBG& sdbg, const size_t& node) {
    std::string label;
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

pair<string, vector<string>> get_systems(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& ordered_cycles,
    int& number_of_spacers
) {
    auto CRISPRArrayNodes = FindCRISPRArrayNodes(sdbg, ordered_cycles);
    auto spacers_nodes = CRISPRArrayNodes.second;
    vector<uint64_t> repeat_nodes = CRISPRArrayNodes.first;

    if (CRISPRArrayNodes.first.empty() || CRISPRArrayNodes.second.empty()) {
        return std::make_pair(string(), vector<string>());
    }

    string repeat = FetchNodeLabel(sdbg, repeat_nodes[0]);
    
    for (size_t i = 1; i < repeat_nodes.size(); i++) {
        std::string node_label = FetchNodeLabel(sdbg, repeat_nodes[i]);

        // Method 1: Using back() method
        char lastChar = node_label.back();  // Get the last character
        std::string lastCharStr(1, lastChar);  // Convert char to string
        
        repeat += lastCharStr;
    }

    
    vector<string> spacers;
    for (const auto& cycle : ordered_cycles) {
        string spacer = "";

        bool passed_repeat = false;
        bool passed_spacers = false;
        int loop_counter = 0;

        // The spacers might wrap around, therefore go through cycle twice
        while (loop_counter < 2) {
            for (const auto& node : cycle) {
                bool is_repeat_node = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
                if (is_repeat_node) {
                    passed_repeat = true;
                    if (passed_spacers) {
                        break;
                    }
                } else {
                    passed_spacers = true;
                }

                if (!is_repeat_node && passed_repeat) {
                    if (spacer == "") {
                        spacer = FetchNodeLabel(sdbg, node);
                    } else {
                        char new_char = FetchNodeLabel(sdbg, node).back();
                        string new_char_str(1, new_char);
                        spacer += new_char_str;
                    }
                }
            }
            loop_counter++;
        }

        if (spacer.size() > 45) {
            spacer = spacer.substr(23, spacer.size() - 45);
        } else {
            spacer.clear();
        }
        spacers.push_back(spacer);
    }
    
    // if(spacers.size() < 3){
        //     return std::make_pair(string(), vector<string>());
        // }
        
    number_of_spacers = spacers.size();
    return std::make_pair(repeat, spacers);
}
