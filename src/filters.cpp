#include "filters.h"

char Filters::_FetchLastCharacter(size_t node) {
    vector<uint8_t> seq(sdbg.k());
    sdbg.GetLabel(node, seq.data());
    return "ACGT"[seq[sdbg.k() - 1]];
}
string Filters::_FetchNodeLabel(size_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}
Filters::Filters(SDBG& sdbg, std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles)
    : sdbg(sdbg), cycles(cycles) {}

    std::vector<uint64_t> Filters::FindRepeatNodePaths(vector<uint64_t> repeat_nodes, uint64_t start_node) {
        uint64_t start = 0;
        vector<uint64_t> all_the_neighbors;

        for (const auto& node : repeat_nodes) {
            uint64_t outgoings[4];
            int num_outgoings = this->sdbg.OutgoingEdges(node, outgoings);
            for (int i = 0; i < num_outgoings; i++) {
                all_the_neighbors.push_back(outgoings[i]);
            }
        }
        for (const auto& node : repeat_nodes)
            if (std::find(all_the_neighbors.begin(), all_the_neighbors.end(), node) == all_the_neighbors.end())
                start = node;

        int maxSize = 0;
        vector<vector<uint64_t>> cycles_per_group = this->cycles[start_node];
        std::vector<uint64_t> arr;
        auto it = std::find(arr.begin(), arr.end(), start);
        int position_to_rotate = std::distance(arr.begin(), it);

        for (int i = 0; i < cycles_per_group.size(); i++) {
            if (cycles_per_group[i].size() > maxSize) {
                maxSize = cycles_per_group[i].size();
                arr = cycles_per_group[i];
            }
        }

        arr.resize(repeat_nodes.size());
        return arr;
    }

    pair<vector<uint64_t>, vector<vector<uint64_t>>> Filters::_FindCRISPRArrayNodes(uint64_t start_node) {
        std::unordered_map<uint64_t, int> element_count;

unordered_map<string, vector<string>> Filters::ListArrays(
    vector<uint64_t> node_order,
    int& number_of_spacers
) {
    unordered_map<string, vector<string>> CRISPRArrays;
    int counter = 0;
    for (const auto& [start_node, _] : cycles) {
        auto CRISPRArrayNodes = _FindCRISPRArrayNodes(start_node);
        auto spacers_nodes = CRISPRArrayNodes.second;
        vector<uint64_t> repeat_nodes = CRISPRArrayNodes.first;
        if (!CRISPRArrayNodes.first.empty() && !CRISPRArrayNodes.second.empty()) {
            string repeat = _FetchNodeLabel(repeat_nodes[0]);
            
            for (size_t i = 1; i < repeat_nodes.size(); i++) {
                std::string node_label = _FetchNodeLabel(repeat_nodes[i]);
                    // Method 1: Using back() method
                char lastChar = node_label.back();  // Get the last character
                 std::string lastCharStr(1, lastChar);  // Convert char to string
              
                repeat += lastCharStr;
            }
            vector<vector<uint64_t>> cycles_nodes =this->cycles[start_node];
            //print the first spacer, first cycle and repeat nodes
            
            vector<string> spacers_temp;
            vector<string> spacers;
            string all_cycles_togehter = _FetchNodeLabel(node_order[0]);
            for (int i = 1; i < node_order.size(); ++i) {
                uint64_t node = node_order[i];
                std::string node_label = _FetchNodeLabel(node);
                // Method 1: Using back() method
                char lastChar = node_label.back();  // Get the last character
                std::string lastCharStr(1, lastChar);  // Convert char to string
                all_cycles_togehter += lastCharStr;
            }
            
            size_t start = 0;
            size_t end;

            // Iterate through the string and find substrings
            while ((end = all_cycles_togehter.find(repeat, start)) != std::string::npos) {
                std::string part = all_cycles_togehter.substr(start, end - start);
                if (!part.empty()) {
                    spacers_temp.push_back(part);
                }
                size_t start = 0;
                size_t end;
                while ((end = all_cycles_togehter.find(repeat, start)) != std::string::npos) {
                    std::string part = all_cycles_togehter.substr(start, end - start);
                    if (!part.empty()) {
                        spacers_temp.push_back(part);
                    }
                    start = end + repeat.size();
                }
                if (start < all_cycles_togehter.size()) {
                    spacers_temp.push_back(all_cycles_togehter.substr(start));
                }
                for (const auto& spacer : spacers_temp) {
                    if (spacer.size() < 23 || spacer.size() > 50)
                        continue;
                    spacers.push_back(spacer.substr(0, spacer.size()));
                    number_of_spacers++;
                }
                if (spacers.size() < 2) {
                    number_of_spacers -= spacers.size();
                    continue;
                }
                CRISPRArrays[repeat] = spacers;
            }
        }
        return CRISPRArrays;
    }

    int Filters::WriteToFile(const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Failed to open file: " + filename);
        }
        int number_of_spacers = 0;
        auto CRISPRArrays = ListArrays(number_of_spacers);
        for (const auto& [repeat, spacers] : CRISPRArrays) {
            file << "Repeat: " << repeat << endl;
            file << "Number of Spacers: " << spacers.size() << endl;
            file << "Spacers:" << endl;
            for (const auto& spacer : spacers) {
                file << spacer << endl;
            }
            file << "----------------------------------" << endl;
        }
        file.close();
        return number_of_spacers;
    }
    
    return CRISPRArrays;
}

int Filters::WriteToFile(vector<uint64_t> node_order, const string& filename) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open file: " + filename);
    }
    int number_of_spacers = 0;
    
    auto CRISPRArrays = ListArrays(node_order, number_of_spacers);
    
    for (const auto& [repeat, spacers] : CRISPRArrays) {
        
        file << "Repeat: " << repeat << endl;
        file << "Number of Spacers: " << spacers.size() << endl;
        file << "Spacers:" << endl;
        for (const auto& spacer : spacers) {
            file << spacer << endl;
        }
        file <<"----------------------------------" << endl; // Add a blank line between different CRISPR arrays

        
    }
    file.close();
    
    return number_of_spacers;

}
