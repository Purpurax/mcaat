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

float get_string_similarity(const string& s1, const string& s2) {
    uint16_t d = get_levenshtein_distance(s1, s2);
    size_t max_size = s1.size() > s2.size() ? s1.size() : s2.size();
    
    return 1.0 - (static_cast<float>(d) / static_cast<float>(max_size));
}

float get_spacer_order_similarity(
    const vector<string>& spacers,
    const string& expected_sequence,
    int variant
) {
    vector<size_t> spacer_positions_in_sequence;
    for (const auto& spacer : spacers) {
        size_t pos = expected_sequence.find(spacer);
        if (pos != string::npos) {
            spacer_positions_in_sequence.push_back(pos);
        } else {
            spacer_positions_in_sequence.push_back(-1);
        }
    }

    if (variant == 0) {
        int correct = 0;
        for (int i = 0; i < spacer_positions_in_sequence.size(); ++i) {
            const size_t position = spacer_positions_in_sequence.at(i);
            
            int amount_of_smaller_elements = 0;
            for (int j = 0; j < spacer_positions_in_sequence.size(); ++j) {
                if (i != j && spacer_positions_in_sequence.at(j) < position) {
                    ++amount_of_smaller_elements;
                }
            }

            if (amount_of_smaller_elements == i) {
                ++correct;
            }
        }
        
        return static_cast<float>(correct) / static_cast<float>(spacer_positions_in_sequence.size());
    } else {
        vector<pair<size_t, size_t>> pairs;
        if (variant == 1) {
            for (int i = 0; i < spacers.size() - 1; ++i) {
                pairs.push_back(std::make_pair(i, i+1));
            }
        } else { // variant == 2
            for (int i = 0; i < spacers.size() - 1; ++i) {
                for (int j = i + 1; j < spacers.size(); ++j) {
                    pairs.push_back(std::make_pair(i, j));
                }
            }
        }

        int correct = 0;
        for (const auto& [part_a, part_b] : pairs) {
            const auto& position_a = spacer_positions_in_sequence.at(part_a);
            const auto& position_b = spacer_positions_in_sequence.at(part_b);

            if (position_a != -1 && position_b != -1 && position_a < position_b) {
                ++correct;
            }
        }

        return static_cast<float>(correct) / static_cast<float>(pairs.size());
    }
}

int get_number_of_duplicate_spacers(
    const vector<string>& spacers,
    const string& expected_sequence
) {
    int result = 0;

    for (const auto& spacer : spacers) {
        int spacer_count = 0;

        size_t pos = 0;
        while ((pos = expected_sequence.find(spacer, pos)) != string::npos) {
            ++spacer_count;
            ++pos;
        }

        if (spacer_count > 1) {
            result += spacer_count - 1;
        }
    }

    return result;
}

string pop_most_similar_sequence(
    const string& sequence,
    vector<string>& choices
) {
    if (choices.empty()) {
        return "";
    }

    float best_similarity = -1.0;
    int best_idx = -1;
    string result_sequence = "";
    for (int i = 0; i < choices.size(); ++i) {
        const string choice_sequence = choices.at(i);
        float current_similarity = get_string_similarity(sequence, choice_sequence);

        if (current_similarity > best_similarity) {
            best_similarity = current_similarity;
            best_idx = i;
            result_sequence = choice_sequence;
        }
    }

    choices.erase(choices.begin() + best_idx);
    return result_sequence;
}
