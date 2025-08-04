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

#include "settings.h"

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
 * @return std::pair<std::string, std::optional<std::string>> If 2 files were found
 * @return std::pair<settings.input_files, std::nullopt> If only 1 file was found
 */
std::pair<std::string, std::optional<std::string>> get_fastq_files_from_settings(
    const Settings& settings
);

#endif
