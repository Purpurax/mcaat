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

int get_cycle_count(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    int count = 0;
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


string fetch_node_label(SDBG& sdbg, const size_t& node) {
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
                    char new_char = fetch_node_label(sdbg, node).back();
                    string new_char_str(1, new_char);
                    repeat += new_char_str;
                    break;
                }
            } else if (passed_spacers) {
                going_through_repeat_region = true;
                if (repeat == "") {
                    repeat = fetch_node_label(sdbg, node);
                } else {
                    char new_char = fetch_node_label(sdbg, node).back();
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
    vector<vector<uint64_t>>& ordered_cycles
) {
    unordered_set<uint64_t> repeat_nodes;

    /* 1. Find all repeat nodes by thresholding */
    const int threshold = static_cast<int>(0.9 * static_cast<float>(ordered_cycles.size()));
    
    unordered_map<uint64_t, int> element_count;
    for (const auto& cycle : ordered_cycles) {
        for (const auto& element : cycle) {
            element_count[element]++;
        }
    }

    for (const auto& cycle : ordered_cycles) {
        for (const auto& element : cycle) {
            if (element_count[element] >= threshold) {
                repeat_nodes.insert(element);
            }
        }
    }

    /* 2. Align the ordered_cycles to some repeat node (most commonly appearing repeat node) */
    /* Finding cycle, probably containing consensus repeat */
    int consensus_repeat_cycle_idx = 0;
    int repeat_nodes_count = 0;
    for (int i = 0; i < ordered_cycles.size(); ++i) {
        const auto cycle = ordered_cycles.at(i);
        
        int current_repeat_nodes_count = 0;
        for (const auto& element : cycle) {
            if (repeat_nodes.find(element) != repeat_nodes.end()) {
                current_repeat_nodes_count += element_count[element];
            }
        }

        if (current_repeat_nodes_count > repeat_nodes_count) {
            repeat_nodes_count = current_repeat_nodes_count;
            consensus_repeat_cycle_idx = i;
        }
    }

    /* Naively finding first repeat node */
    const auto& consensus_cycle = ordered_cycles.at(consensus_repeat_cycle_idx);
    uint64_t start_repeat_node = std::numeric_limits<uint64_t>::max();
    for (int i = 0; i < consensus_cycle.size(); ++i) {
        const auto before_before_node = consensus_cycle.at((i - 2 + consensus_cycle.size()) % consensus_cycle.size());
        const auto before_node = consensus_cycle.at((i - 1 + consensus_cycle.size()) % consensus_cycle.size());
        const auto node = consensus_cycle.at(i);

        if (repeat_nodes.find(node) != repeat_nodes.end()) {
            if (start_repeat_node == std::numeric_limits<uint64_t>::max()) {
                // Have some repeat node in the variable
                start_repeat_node = node;
            } else if (repeat_nodes.find(before_node) == repeat_nodes.end()
                && repeat_nodes.find(before_before_node) == repeat_nodes.end()) {
                // Before and Before_before are spacers
                start_repeat_node = node;
            }
        }
    }

    /* Align cycles (that can be aligned) */
    vector<int> aligned_cycles_indices;
    for (int i = 0; i < ordered_cycles.size(); ++i) {
        const auto& cycle = ordered_cycles.at(i);
        auto start_repeat_node_position = std::find(cycle.begin(), cycle.end(), start_repeat_node);
        if (start_repeat_node_position == cycle.end()) {
            continue;
        }
        aligned_cycles_indices.push_back(i);
        
        vector<uint64_t> aligned_cycle;
        int offset = std::distance(cycle.begin(), start_repeat_node_position);
        for (int j = 0; j < cycle.size(); ++j) {
            const auto node = cycle.at((j + offset) % cycle.size());
            aligned_cycle.push_back(node);
        }
        ordered_cycles[i] = aligned_cycle;
    }

    /* 3. Turn cycles into sequence */
    vector<string> cycles_str;
    for (const auto& cycle : ordered_cycles) {
        string cycle_str = "";
        for (const auto& node : cycle) {
            const auto label = fetch_node_label(sdbg, node);
            cycle_str.push_back(label.back());
        }
        cycles_str.push_back(cycle_str);
    }

    /* 4. Barrel shift until the sequence is: RRRRRSSSSS */
    /* Get amount of repeat character matches left and right */
    const string& consensus_cycle_str = cycles_str.at(consensus_repeat_cycle_idx);
    vector<int> repeat_chr_matches_left;
    vector<int> repeat_chr_matches_right;
    for (const auto& idx : aligned_cycles_indices) {
        if (idx == consensus_repeat_cycle_idx) {
            continue;
        }

        const auto& cycle_str = cycles_str.at(idx);

        int current_match_right = 0;
        for (int i = 0; i < cycle_str.size(); ++i) {
            const char current = cycle_str.at(i);
            const char consensus = consensus_cycle_str.at(i);

            if (current != consensus) {
                break;
            }
            ++current_match_right;
        }
        repeat_chr_matches_right.push_back(current_match_right);

        int current_match_left = 0;
        for (int i = 0; i < cycle_str.size(); ++i) {
            const char current = cycle_str.at(cycle_str.size() - i - 1);
            const char consensus = consensus_cycle_str.at(consensus_cycle_str.size() - i - 1);

            if (current != consensus) {
                break;
            }
            ++current_match_left;
        }
        repeat_chr_matches_left.push_back(current_match_left);
    }

    /* Take maximum repeat match that is in 20% of all matches */
    std::sort(repeat_chr_matches_left.begin(), repeat_chr_matches_left.end());
    std::sort(repeat_chr_matches_right.begin(), repeat_chr_matches_right.end());

    int quart_idx = static_cast<int>(
        std::round(0.2 * static_cast<float>(repeat_chr_matches_left.size())));
    if (quart_idx >= repeat_chr_matches_left.size()) {
        quart_idx = repeat_chr_matches_left.size() - 1;
    }

    const int repeat_chr_match_left = repeat_chr_matches_left[quart_idx];
    const int repeat_chr_match_right = repeat_chr_matches_right[quart_idx];

    /* 5. Find consensus repeat sequence */
    string repeat = "";
    int repeat_length = repeat_chr_match_left + repeat_chr_match_right;
    for (int i = 0; i < repeat_length; ++i) {
        int offset = consensus_cycle_str.size() - repeat_chr_match_left + i;
        const char& chr = consensus_cycle_str.at(offset % consensus_cycle_str.size());

        repeat.push_back(chr);
    }

    /* 6. Reconstruct full sequence and spacer strings */
    vector<string> spacers;
    string full_sequence = "";
    for (int i = 0; i < cycles_str.size(); ++i) {
        string spacer = "";
        const auto& cycle_str = cycles_str.at(i);
        
        int offset;
        if (std::find(aligned_cycles_indices.begin(), aligned_cycles_indices.end(), i) != aligned_cycles_indices.end()) {
            // Already aligned, thereby using repeat_chr_match_left
            offset = cycle_str.size() - repeat_chr_match_left;
        } else {
            // Not aligned, thereby need to find offset
            offset = 0;
            int best_string_distance = std::numeric_limits<int>::max();

            for (int j = 0; j < cycle_str.size(); ++j) {
                string rotated_cycle_str = cycle_str.substr(j) + cycle_str.substr(0, j);

                string supposed_repeat = rotated_cycle_str.substr(0, repeat_length);
                uint16_t distance = get_levenshtein_distance(supposed_repeat, repeat);
                if (distance < best_string_distance) {
                    best_string_distance = distance;
                    offset = j;
                }
            }
        }

        for (int j = 0; j < cycle_str.size(); ++j) {
            const auto& chr = cycle_str.at((offset + j) % cycle_str.size());

            if (j >= repeat_length) {
                spacer += chr;
            }
            full_sequence += chr;
        }

        spacers.push_back(spacer);
    }

    /* crispr sequence end with a repat (usually mutated), tmp taking consensus repeat */
    full_sequence += repeat;

    return std::make_tuple(repeat, spacers, full_sequence);
}

uint16_t get_levenshtein_distance(const string& s1, const string& s2) {
    vector<vector<uint16_t>> dist;

    /* Setup matrix dist */
    vector<uint16_t> first_row;
    for (int i = 0; i < s1.size() + 1; ++i) {
        first_row.push_back(i);
    }
    dist.push_back(first_row);

    for (int i = 1; i < s2.size() + 1; ++i) {
        vector<uint16_t> row;
        row.push_back(i);
        dist.push_back(row);
    }

    /* Fill matrix with computation */
    for (int x = 1; x < s1.size() + 1; ++x) {
        for (int y = 1; y < s2.size() + 1; ++y) {
            char s1_char = s1.at(x - 1);
            char s2_char = s2.at(y - 1);

            uint16_t min;
            if (s1_char == s2_char) { // No extra distance
                min = dist.at(y - 1).at(x - 1);
            } else { // update(s1_char, s2_char)
                min = dist.at(y - 1).at(x - 1) + 1;
            }

            if (min > dist.at(y).at(x - 1) + 1) { // insert(s1_char)
                min = dist.at(y).at(x - 1) + 1;
            }

            if (min > dist.at(y - 1).at(x) + 1) { // delete()
                min = dist.at(y - 1).at(x) + 1;
            }

            dist.at(y).push_back(min);
        }
    }

    /* Return result of matrix */
    return dist.at(s2.size()).at(s1.size());
}
