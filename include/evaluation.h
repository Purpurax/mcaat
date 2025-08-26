#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

/**
 * @internal
 * @brief Get the levenshtein distance of both strings
 * 
 * Using the Levenshtein distance (matrix implementation in O(|s1| |s2|)) with:
 *  - insert(c)   costs 1
 *  - update(a,b) costs 1
 *  - delete(c)   costs 1
 * 
 * @param s1 First string
 * @param s2 Second string
 * @return Distance according to the algorithm (uint16_t)
 */
uint16_t get_levenshtein_distance(const string& s1, const string& s2);

// /**
//  * @internal
//  * @brief Get the similarity of the strings
//  * 
//  * Gets the levenstheins distance and returns the ratio compared to the size.
//  * So for |s1| = n, |s2| = m, and levensthein distance d
//  *  the result = 1 - (d / max(n, m))
//  * 
//  * @todo Use Hamming distance for strings of equal length and just in general make it efficient
//  * 
//  * @param s1 First string
//  * @param s2 Second string
//  * 
//  * @return The value between 0 and 1 for how similar the strings are (float)
//  */
// float get_similarity_of_strings(const string& s1, const string& s2);

// /**
//  * @internal
//  * @brief Splitting the main_string at each substring with some similarity
//  * 
//  * @note If the similarity threshold has been reached but there is a better
//  * match just one character after, it will split at the better match.
//  * 
//  * @note It will remove exactly |split_string| characters at any split.
//  * 
//  * @param main_string The string to be splitted
//  * @param split_string The string used to perform the split
//  * @param similarity_threshold How similar a match must be for a split
//  * 
//  * @return The result with the remaining parts(vector<string>)
//  */
// vector<string> split_at_with_similarity_threshold(
//     const string& main_string,
//     const string& split_string,
//     const float& similarity_threshold
// );

/**
 * @internal
 * @brief Compares strings using levenshtein distance
 * 
 * With |s1| = n, |s2| = m, levenshtein distance d
 * The result is 1 - (d / max(n, m))
 * 
 * @param s1 First sequence
 * @param s2 Second sequence
 * 
 * @return Value between 0 and 1 for how similar the sequences are (float)
 */
float get_similarity(const string& s1, const string& s2);

/**
 * @internal
 * @brief Computes the similarities of the actual and expected sequences
 * 
 * If for one repeat multiple sequences exist on both sides, the exact computation
 * would involve checking every combination of those sequences. As this is not efficient,
 * the algorithm chooses some pair by the greedy highest similarity, and repeats until
 * all sequences are matched. => This underestimates the similarity value
 * 
 * @note The algorithm will look for matching sequences.
 * If an actual sequence has no match, the value is -1
 * If an expected sequence has no match, the value is -2
 * 
 * @param actual_repeat_seq_map A map from repeat to a vector of crispr sequences
 * @param expected_repeat_seq_map A map from repeat to a vector of crispr sequences
 * 
 * @return A map from repeat to a vector of similarity values (unordered_map<string, vector<float>>)
 */
unordered_map<string, vector<float>> compute_the_sequence_similarities(
    const unordered_map<string, vector<string>>& actual_repeat_seq_map,
    const unordered_map<string, vector<string>>& expected_repeat_seq_map
);

/**
 * @brief Gives the results of a comparison
 * 
 * @param actual_repeat_seq_map Map of repeat sequence to full crispr sequence
 * @param expected_repeat_seq_map Map of repeat sequence to full crispr sequence
 * @param amount_of_non_matching_sequences Matching is done by pairing up crispr sequences
 * that have the same repeat. This gives the number of non matches.
 * @param accuracy_for_repeat Gives the underestimating similarity for each repeat
 * @param total_average_accuracy Computes all values together to give one accuracy value
 */
void compare_to_ground_of_truth(
    const unordered_map<string, vector<string>>& actual_repeat_seq_map,
    const unordered_map<string, vector<string>>& expected_repeat_seq_map,
    size_t& amount_of_non_matching_sequences,
    unordered_map<string, float>& accuracy_for_repeat,
    float& total_average_accuracy
);

#endif
