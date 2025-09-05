/**
 * @file tmp_utils.h
 * @brief Temporary Helper functions
 * 
 * Here are all functions that still have to be integrated into the full codebase
 * but land here because of various reasons
 */
#ifndef INCLUDE_TMP_UTILS_H_
#define INCLUDE_TMP_UTILS_H_

#include <iostream>
#include <string>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

#include "settings.h"
#include <sdbg/sdbg.h>

using namespace std;

/**
 * @brief Removes any spaces, new lines, and tabs before and after the string
 * 
 * @post The parameter s will be trimmed.
 * 
 * @param s The string to trim
 */
void trim_string(string& s);

/**
 * @brief Get the fastq files from the settings.input_files string
 * 
 * It takes the settings.input_files and looks for a space. If some exists, it expects 2 files.
 * If not only one file is expected.
 * After Splitting it into 2 files, or just taking out the 1 file, they are trimmed and returned.
 * 
 * @param settings Which contains at least one file path, or two space separated
 * 
 * @return pair<string, optional<string>> If 2 files were found
 * @return pair<settings.input_files, nullopt> If only 1 file was found
 */
pair<string, optional<string>> get_fastq_files_from_settings(
    const Settings& settings
);

/**
 * @brief Get the cycle count from the special object
 * 
 * @return The cycle count (int)
 */
int get_cycle_count(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
);

/**
 * @brief Loads in the sdbg graph and prints the nodes and edges as dot file
 * 
 * @param lib_file_path The folder in which the files are contained (data.lib, graph.sdbg.0, ..., graph.sdbg_info)
 */
void print_sdbg_graph_to_dot_file_convert(const string& lib_file_path);

/**
 * @brief Gets the repeat and the spacers in correct order as sequences
 * 
 * 1. Find all repeat nodes by thresholding
 * 2. Align the cycles to some repeat node (most commonly appearing repeat node)
 * 3. Turn cycles into sequence
 * 4. Barrel shift until the sequence is: RRRRRSSSSS
 * 5. Find consensus repeat sequence
 * 6. Reconstruct full sequence and spacer strings
 * 
 * @post ordered_cycles will be modified internally
 * 
 * @param sdbg The graph used for getting the sequences and the repeats
 * @param ordered_cycles The cycles in the correct order containing node ids
 * 
 * @return Returns the repeat sequence, spacer sequences, and the full sequence (tuple<string, vector<string>>, string)
 */
tuple<string, vector<string>, string> get_systems(
    SDBG& sdbg,
    vector<vector<uint64_t>>& ordered_cycles
);

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

#endif
