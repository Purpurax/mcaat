#include "evaluation.h"

float get_string_similarity(const string& s1, const string& s2) {
    uint16_t d = get_levenshtein_distance(s1, s2);
    size_t max_size = s1.size() > s2.size() ? s1.size() : s2.size();
    
    return 1.0 - (static_cast<float>(d) / static_cast<float>(max_size));
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

string get_most_similar_sequence(
    const string& sequence,
    const vector<string>& choices
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

    // Remove the choice from the compares:
    // choices.erase(choices.begin() + best_idx);

    return result_sequence;
}
