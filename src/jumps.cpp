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
        // Check if the library files exist
        std::string lib_prefix = settings.graph_folder + "/outfile_prefix";
        std::string lib_info_file = lib_prefix + ".lib_info";
        std::string lib_bin_file = lib_prefix + ".bin";
        
        if (!std::filesystem::exists(lib_info_file)) {
            std::cerr << "Library info file not found: " << lib_info_file << std::endl;
            return jumps;
        }
        
        if (!std::filesystem::exists(lib_bin_file)) {
            std::cerr << "Library binary file not found: " << lib_bin_file << std::endl;
            return jumps;
        }

        // Initialize SeqPackage and SequenceLibCollection properly
        SeqPackage data_holder;
        SequenceLibCollection seq_libs;
        seq_libs.SetPath(lib_prefix);  // Set path without extension
        
        // Read the library
        seq_libs.Read(&data_holder, false);  // Don't reverse sequences
        std::cout << "Successfully loaded sequence library" << std::endl;
        
        for (size_t lib_id = 0; lib_id < seq_libs.size(); ++lib_id) {
            const SequenceLib& lib = seq_libs.GetLib(lib_id);
            
            for (size_t i = 0; i < lib.seq_count(); ++i) {
                auto seq_view = lib.GetSequenceView(i);
                
                // Build sequence string more efficiently
                std::string sequence_str;
                sequence_str.reserve(seq_view.length());
                for (unsigned j = 0; j < seq_view.length(); ++j) {
                    sequence_str += "ACGT"[seq_view.base_at(j)];
                }
                
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
                    jumps.push_back(jump);
                }
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error reading sequences: " << e.what() << std::endl;
    }

    std::cout << "Generated " << jumps.size() << " jumps from reads" << std::endl;
    return jumps;
}
