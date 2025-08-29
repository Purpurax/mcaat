#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <stdint.h>

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
 * @brief Returns a similarity value for each sequences in actual_sequences
 * 
 * Overall it overestimates the similarity values using greedy. It cannot be exact
 * as it is difficult to find the matching sequences from actual and expected.
 * 
 * If actual_sequences.size() > expected_sequences.size() then
 * the last few entries will have -1.0
 * 
 * If actual_sequences.size() < expected_sequences.size() then
 * the rest of the expected are ignored and not used
 * 
 * @post expected_sequence will be modified and elements will be removed
 * 
 * @param actual_sequences The actual result
 * @param expected_sequences The expected true result
 * 
 * @return A map from sequence to (overestimated) similarity (unordered_map<string, float>) 
 */
unordered_map<string, float> compare_to_ground_of_truth(
    const vector<string>& actual_sequences,
    vector<string>& expected_sequences
);

#endif
