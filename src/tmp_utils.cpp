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


string FetchNodeLabel(SDBG& sdbg, const size_t& node) {
    string label;
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

pair<vector<uint64_t>, string> Find_CRISPR_repeat_nodes_and_sequence(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& ordered_cycles
) {
    const int threshold = static_cast<int>(0.9 * static_cast<float>(ordered_cycles.size()) + 0.5);
    
    unordered_map<uint64_t, int> element_count;
    for (const auto& cycle : ordered_cycles) {
        for (const auto& element : cycle) {
            element_count[element]++;
        }
    }

    // Get all repeat_nodes
    vector<uint64_t> repeat_nodes;
    unordered_set<uint64_t> repeat_nodes_set;
    for (const auto& cycle : ordered_cycles) {
        bool passed_spacers = false;
        bool going_through_repeat_region = false;
        int loop_counter = 0;

        // The spacers might wrap around, therefore go through cycle twice:
        // Starting with S -> SSSSSRRRRRSSSSS SSSSSRRRRRSSSSS
        // Starting with R -> RRRRRSSSSSRRRRR RRRRRSSSSSRRRRR
        while (loop_counter < 2) {
            for (const auto& node : cycle) {
                bool is_spacer_node = element_count[node] < threshold;
                if (is_spacer_node) {
                    passed_spacers = true;
                    if (going_through_repeat_region) {
                        break;
                    }
                } else if (passed_spacers) {
                    going_through_repeat_region = true;
                    bool already_in_set = repeat_nodes_set.insert(node).second;
                    if (!already_in_set) {
                        repeat_nodes.push_back(node);
                    }
                }
            }

            loop_counter++;
        }
        for (const auto& [element, count] : element_count) {
            if (static_cast<float>(count) >= threshold) {
                repeat_nodes.push_back(element);
            }
        }
    }

    // Choose repeat node sequence that occur the most
    int repeat_nodes_count = 0;
    int best_idx = 0;
    for (int i = 0; i < ordered_cycles.size(); ++i) {
        const auto cycle = ordered_cycles.at(i);

        int current_repeat_nodes_count = 0;
        for (const auto& node : cycle) {
            bool is_repeat = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
            if (is_repeat) {
                current_repeat_nodes_count += element_count[node];
            }
        }

        if (current_repeat_nodes_count > repeat_nodes_count) {
            repeat_nodes_count = current_repeat_nodes_count;
            best_idx = i;
        }
    }

    // Fetch best sequence of repeat nodes
    string repeat = "";
    bool passed_spacers = false;
    bool going_through_repeat_region = false;
    int loop_counter = 0;

    // The spacers might wrap around, therefore go through cycle twice:
    // Starting with S -> SSSSSRRRRRSSSSS SSSSSRRRRRSSSSS
    // Starting with R -> RRRRRSSSSSRRRRR RRRRRSSSSSRRRRR
    while (loop_counter < 2) {
        for (const auto& node : ordered_cycles.at(best_idx)) {
            bool is_repeat_node = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
            if (!is_repeat_node) {
                passed_spacers = true;
                if (going_through_repeat_region) {
                    // Idk why, but one letter is missing
                    char new_char = FetchNodeLabel(sdbg, node).back();
                    string new_char_str(1, new_char);
                    repeat += new_char_str;
                    break;
                }
            } else if (passed_spacers) {
                going_through_repeat_region = true;
                if (repeat == "") {
                    repeat = FetchNodeLabel(sdbg, node);
                } else {
                    char new_char = FetchNodeLabel(sdbg, node).back();
                    string new_char_str(1, new_char);
                    repeat += new_char_str;
                }
            }
        }

        loop_counter++;
    }

    return {repeat_nodes, repeat};
}

tuple<string, vector<string>, string> get_systems(
    SDBG& sdbg,
    const vector<vector<uint64_t>>& ordered_cycles
) {
    auto [repeat_nodes, repeat] = Find_CRISPR_repeat_nodes_and_sequence(sdbg, ordered_cycles);
    
    vector<string> spacers;
    for (const auto& cycle : ordered_cycles) {
        string spacer = "";

        bool passed_repeat = false;
        bool going_through_spacer_region = false;
        int loop_counter = 0;

        // The spacers might wrap around, therefore go through cycle twice:
        // Starting with S -> SSSSSRRRRRSSSSS SSSSSRRRRRSSSSS
        // Starting with R -> RRRRRSSSSSRRRRR RRRRRSSSSSRRRRR
        while (loop_counter < 2) {
            for (const auto& node : cycle) {
                bool is_repeat_node = std::find(repeat_nodes.begin(), repeat_nodes.end(), node) != repeat_nodes.end();
                if (is_repeat_node) {
                    passed_repeat = true;
                    if (going_through_spacer_region) {
                        break;
                    }
                } else if (passed_repeat) {
                    going_through_spacer_region = true;
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
    
    string crispr_sequence = repeat; // roughly R-S1-R-S2-...-R-SN-R
    for (auto spacer_it = spacers.rbegin(); spacer_it != spacers.rend(); ++spacer_it) {
        crispr_sequence += *spacer_it;
        crispr_sequence += repeat;
    }

    return std::make_tuple(repeat, spacers, crispr_sequence);
}
