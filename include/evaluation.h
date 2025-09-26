#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <stdint.h>

#include "tmp_utils.h"

using namespace std;

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
float get_string_similarity(const string& s1, const string& s2);

/**
 * @brief Get the spacer order similarity by variant
 * 
 * Best explained with an example:
 *  - spacers: A B C D
 *  - expected_sequence: _B_C_A_D_
 * with ABCD are spacers and _ is some other sequence
 * 
 * Variant 0: Comparing position by position
 *  => N N N J
 *  => result is 25%
 * 
 * Variant 1: Compare each pair of spacers and if their order is in the expected
 *  - spacer pairs: A->B, B->C, C->D
 *  => N J J
 *  => result is 66.6%
 * 
 * Variant 2: Compare each possible pair of spacers and if their order is in there
 *  - spacer pairs: A->B, A->C, A->D, B->C, B->D, C->D
 *  => N N J J J J
 *  => result is 66.6%
 * 
 * @note if the expected sequence contains duplicate spacer entries,
 * the results become unpredictable.
 * 
 * @param spacers 
 * @param expected_sequence 
 * @param variant Number specifying the variant (0, 1, 2)
 * 
 * @return 0 <= similarity <= (float)
 */
float get_spacer_order_similarity(
    const vector<string>& spacers,
    const string& expected_sequence,
    int variant
);

/**
 * @brief Get the number of duplicate spacers
 * 
 * Simply counts the amount of occurences of every spacer.
 * If it occurs more than once, it adds it to the result
 * 
 * @param spacers 
 * @param expected_sequence 
 * 
 * @return Amount of spacer duplicates (int)
 */
int get_number_of_duplicate_spacers(
    const vector<string>& spacers,
    const string& expected_sequence
);

/**
 * @internal
 * @brief Chooses the most similar sequence and removes it from the vector
 * 
 * @post choices will be modified and one elements will be removed (if possible)
 * 
 * @param sequence Used to compare and find the most similar
 * @param choices Where to pop the sequence from
 * 
 * @return The most similar sequence (string)
 * @return "", if choices is empty
 */
string get_most_similar_sequence(
    const string& sequence,
    vector<string>& choices
);

#endif
