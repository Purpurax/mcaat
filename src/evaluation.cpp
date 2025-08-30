#include "evaluation.h"

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

unordered_map<string, float> compare_to_ground_of_truth(
    const vector<string>& actual_sequences,
    vector<string>& expected_sequences
) {
    unordered_map<string, float> result;

    for (const auto& actual_seq : actual_sequences) {
        if (expected_sequences.size() == 0) {
            break;
        }
        
        // finding match from actual to expected (overestimating in total)
        float best_similarity = -1.0;
        int best_idx = 0;
        for (int i = 0; i < expected_sequences.size(); ++i) {
            string expected_sequence = expected_sequences.at(i);
            float current_similarity = get_similarity(actual_seq, expected_sequence);

            if (current_similarity > best_similarity) {
                best_similarity = current_similarity;
                best_idx = i;
            }
        }

        result.insert({actual_seq, best_similarity});
        expected_sequences.erase(expected_sequences.begin() + best_idx);
    }

    return result;
}
