#include "jumps.h"

unordered_map<string, uint64_t> get_kmer_to_node_map(const SDBG& sdbg) {
    const uint32_t K = sdbg.k();
    unordered_map<string, uint64_t> kmer_to_node;

    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        if (!sdbg.IsValidEdge(node)) continue;

        uint8_t seq_data[K];
        sdbg.GetLabel(node, seq_data);

        string node_label;
        for (int i = K - 1; i >= 0; --i) {
            node_label.append(1, "ACGT"[seq_data[i] - 1]);
        }
        std::reverse(node_label.begin(), node_label.end());
        kmer_to_node[node_label] = node;
    }

    return kmer_to_node;
}

vector<string> extract_sequences_from_fastq_file(const string& fastq_file_path) {
    vector<string> sequences;

    try {
        klibpp::SeqStreamIn iss(fastq_file_path.c_str());
        auto file_records = iss.read();
        for (const auto& record : file_records) {
            sequences.push_back(record.seq);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading file \"" << fastq_file_path;
        std::cerr << "\" sequences because: " << e.what() << std::endl;
    }

    return sequences;
}

string reverse_pair_ends_sequence(string sequence) {
    std::reverse(sequence.begin(), sequence.end());

    for (char& base : sequence) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
        }
    }

    return sequence;
}

std::optional<Jump> create_jump_from_sequence(
    const string& sequence,
    const unordered_map<string, uint64_t>& kmer_to_node_id,
    const uint32_t K
) {
    const size_t seq_length = sequence.length();
    if (seq_length < 2 * K) {
        return std::nullopt;
    }

    const string start_k_mer = sequence.substr(0, K);
    const string end_k_mer = sequence.substr(seq_length - K, K);

    std::optional<uint64_t> start_node_id;
    std::optional<uint64_t> end_node_id;

    auto start_it = kmer_to_node_id.find(start_k_mer);
    if (start_it != kmer_to_node_id.end()) {
        start_node_id = start_it->second;
    }

    auto end_it = kmer_to_node_id.find(end_k_mer);
    if (end_it != kmer_to_node_id.end()) {
        end_node_id = end_it->second;
    }

    if (start_node_id.has_value() && end_node_id.has_value()) {
        Jump jump;
        jump.start_k_mer_id = start_node_id.value();
        jump.end_k_mer_id = end_node_id.value();
        jump.nodes_in_between = seq_length - K - 1;
        return std::make_optional(jump);
    }
    return std::nullopt;
}

vector<Jump> get_jumps_from_reads(
    const SDBG& sdbg,
    const string& fastq_file_1,
    const optional<string>& fastq_file_2,
    const size_t thread_count
) {
    vector<Jump> jumps;
    const uint32_t K = sdbg.k();

    const unordered_map<string, uint64_t> kmer_to_node_id = get_kmer_to_node_map(sdbg);
    
    vector<string> sequences = extract_sequences_from_fastq_file(fastq_file_1);
    if (fastq_file_2.has_value()) {
        const auto sequences_2 = extract_sequences_from_fastq_file(fastq_file_2.value());
        for (const auto& sequence : sequences_2) {
            const string reversed_seq = reverse_pair_ends_sequence(sequence);
            sequences.push_back(reversed_seq);
        }
    }

    #pragma omp parallel num_threads(thread_count)
    {
        #pragma omp for
        for (const auto& sequence : sequences) {
            const std::optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_id, K);
            if (jump_opt.has_value()) {
                #pragma omp critical
                jumps.push_back(jump_opt.value());
            }

            // Experimentally trying to reduce the jump size by splitting sequence
            size_t half = sequence.length() / 2;
            string seq1 = sequence.substr(0, half);
            string seq2 = sequence.substr(half);

            const std::optional<Jump> jump_opt1 = create_jump_from_sequence(seq1, kmer_to_node_id, K);
            if (jump_opt1.has_value()) {
                #pragma omp critical
                jumps.push_back(jump_opt1.value());
            }

            const std::optional<Jump> jump_opt2 = create_jump_from_sequence(seq2, kmer_to_node_id, K);
            if (jump_opt2.has_value()) {
                #pragma omp critical
                jumps.push_back(jump_opt2.value());
            }
        }
    }

    return jumps;
}
