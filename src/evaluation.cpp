#include "evaluation.h"

// float get_similarity_of_strings(const string& s1, const string& s2) {
//     uint16_t n = s1.size();
//     uint16_t m = s2.size();
//     uint16_t d = get_levenshtein_distance(s1, s2);

//     uint16_t max = n ? n > m : m;

//     return 1.0 - (static_cast<float>(d) / static_cast<float>(max));
// }

// vector<string> split_at_with_similarity_threshold(
//     const string& main_string,
//     const string& split_string,
//     const float& similarity_threshold
// ) {
//     size_t split_len = split_string.size();

//     vector<float> similarities;
//     for (int i = 0; i <= main_string.size() - split_len; ++i) {
//         uint16_t matching_characters_count = 0;
        
//         for (int j = i; j < i + split_len; ++j) {
//             if (main_string.at(j) == split_string.at(j - i)) {
//                 matching_characters_count++;
//             }
//         }

//         float similarity_to_repeat = static_cast<float>(matching_characters_count);
//         similarity_to_repeat /= static_cast<float>(split_len);

//         similarities.push_back(similarity_to_repeat);
//     }

//     vector<string> result;
//     string result_buffer = "";
//     uint16_t skip_n_characters = 0;

//     for (int i = 0; i < main_string.size(); ++i) {
//         if (skip_n_characters > 0) {
//             skip_n_characters--;
//             continue;
//         }

//         // Split, if similarity is over threshold and next similarity is not bigger
//         if (i < similarities.size() && similarities.at(i) >= similarity_threshold
//         && (i + 1 >= similarities.size() || similarities.at(i) >= similarities.at(i + 1))) {
//             skip_n_characters = split_len - 1; // First character already ignored
//             if (result_buffer != "") {
//                 result.push_back(result_buffer);
//                 result_buffer = "";
//             }
//         } else { // No Split
//             result_buffer += main_string.at(i);
//         }
//     }

//     if (result_buffer != "") {
//         result.push_back(result_buffer);
//     }

//     return result;
// }

float get_similarity(const string& s1, const string& s2) {
    uint16_t d = get_levenshtein_distance(s1, s2);
    size_t max_size = s1.size() > s2.size() ? s1.size() : s2.size();
    
    return 1.0 - (static_cast<float>(d) / static_cast<float>(max_size));
}

float compare_sequence_to_ground_of_truth(
    const string& sequence,
    vector<string>& expected_sequences
) {
    float best_similarity = -1.0;
    int best_idx = -1;
    for (int i = 0; i < expected_sequences.size(); ++i) {
        string expected_sequence = expected_sequences.at(i);
        float current_similarity = get_similarity(sequence, expected_sequence);

        if (current_similarity > best_similarity) {
            best_similarity = current_similarity;
            best_idx = i;
        }
    }

    if (best_idx != -1) {
        expected_sequences.erase(expected_sequences.begin() + best_idx);
    }
    return best_similarity;
}
