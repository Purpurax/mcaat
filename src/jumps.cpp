#include "jumps.h"

std::vector<Jump> get_jumps_from_reads(SDBG& sdbg, Settings& settings) {
    std::vector<Jump> jumps;
    const size_t K = 23;
    const uint64_t INVALID_VERTEX_ID = -1;

    // Build k-mer lookup table once (for efficiency)
    std::unordered_map<std::string, uint64_t> kmer_to_node;
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        if (!sdbg.IsValidEdge(node)) continue;
        uint8_t seq_data[K];
        sdbg.GetLabel(node, seq_data);
        std::string node_label;
        for (int i = K - 1; i >= 0; --i) {
            node_label.append(1, "ACGT"[seq_data[i] - 1]);
        }
        std::reverse(node_label.begin(), node_label.end());
        kmer_to_node[node_label] = node;
    }

    try {
        std::vector<std::string> sequences;
        if (settings.input_files.find(" ") != std::string::npos) {
            size_t space_pos = settings.input_files.find(" ");
            std::string input_file1 = settings.input_files.substr(0, space_pos);
            std::string input_file2 = settings.input_files.substr(space_pos + 1);

            // Strip leading/trailing whitespace
            input_file1.erase(0, input_file1.find_first_not_of(" \t\n\r"));
            input_file1.erase(input_file1.find_last_not_of(" \t\n\r") + 1);
            input_file2.erase(0, input_file2.find_first_not_of(" \t\n\r"));
            input_file2.erase(input_file2.find_last_not_of(" \t\n\r") + 1);

            klibpp::SeqStreamIn iss1(input_file1.c_str());
            auto file_records1 = iss1.read();
            for (const auto& record : file_records1) {
                sequences.push_back(record.seq);
            }

            klibpp::SeqStreamIn iss2(input_file2.c_str());
            auto file_records2 = iss2.read();
            for (const auto& record : file_records2) {
                std::string reversed_seq = record.seq;
                std::reverse(reversed_seq.begin(), reversed_seq.end());
                // Reverse complement: reverse and complement each base
                for (char& base : reversed_seq) {
                    switch (base) {
                        case 'A': base = 'T'; break;
                        case 'T': base = 'A'; break;
                        case 'C': base = 'G'; break;
                        case 'G': base = 'C'; break;
                    }
                }
                sequences.push_back(reversed_seq);
            }
        } else {
            klibpp::SeqStreamIn iss(settings.input_files.c_str());
            auto file_records = iss.read();
            for (const auto& record : file_records) {
                sequences.push_back(record.seq);
            }
        }

        #pragma omp parallel num_threads(settings->threads_count)
        {
            #pragma omp for
            for (const auto& sequence_str : sequences) {
                size_t seq_length = sequence_str.length();
                if (seq_length < 2 * K) continue;
                
                std::string start_k_mer = sequence_str.substr(0, K);
                std::string end_k_mer = sequence_str.substr(seq_length - K, K);
                
                // Use lookup table instead of linear search
                uint64_t start_node_id = INVALID_VERTEX_ID;
                uint64_t end_node_id = INVALID_VERTEX_ID;
                
                auto start_it = kmer_to_node.find(start_k_mer);
                if (start_it != kmer_to_node.end()) {
                    start_node_id = start_it->second;
                }
                
                auto end_it = kmer_to_node.find(end_k_mer);
                if (end_it != kmer_to_node.end()) {
                    end_node_id = end_it->second;
                }
                
                if (start_node_id != INVALID_VERTEX_ID && end_node_id != INVALID_VERTEX_ID) {
                    Jump jump;
                    jump.start_k_mer_id = start_node_id;
                    jump.end_k_mer_id = end_node_id;
                    jump.nodes_in_between = (seq_length >= 2 * K - 1) ? 
                                            (seq_length - K + 1 - 2) : 0;
                    #pragma omp critical
                    {
                        jumps.push_back(jump);
                    }
                // } else {
                //     std::cout << "SEQUENCE ## " << sequence_str << " failed to be a jump" << std::endl;
                //     std::cout << "  start_k_mer: " << start_k_mer << " -> node_id: " << start_node_id << std::endl;
                //     std::cout << "  end_k_mer:   " << end_k_mer   << " -> node_id: " << end_node_id << std::endl;
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading sequences: " << e.what() << std::endl;
    }

    return jumps;
}
