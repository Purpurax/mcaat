#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <string>
#include <map>
#include <thread>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>

using namespace std;
namespace fs = std::filesystem;

struct Settings {
    std::string input_files;
    double ram = 0.0; // In gigabytes
    size_t threads = 0;
    std::string output_folder = "out";
    std::string graph_folder;
    std::string cycles_folder;
    std::string output_file;

    Settings() {
        // Generate timestamp for folder structure
        string timestamp = get_timestamp();
        output_folder = "" + timestamp;
        graph_folder = output_folder + "/graph";
        cycles_folder = output_folder + "/cycles";
        output_file = output_folder + "/CRISPR_Arrays.txt";
    }

    // Generate timestamp in YYYY-MM-DD_HH-MM-SS format
    string get_timestamp() {
        auto now = chrono::system_clock::now();
        auto in_time_t = chrono::system_clock::to_time_t(now);
        stringstream ss;
        ss << put_time(localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
        return ss.str();
    }

    // Get maximum available threads
    size_t max_threads() const {
        return thread::hardware_concurrency();
    }

    // Validate settings and return status map
    map<string, pair<bool, string>> validate_settings() const {
        map<string, pair<bool, string>> settings_message_map;
        
        // Check input files
        bool input_valid = !input_files.empty();
        settings_message_map["Input Files"] = make_pair(input_valid, 
            input_valid ? input_files + " exist(s)" : "No input files specified");

        // Check RAM (using reasonable default max of 128GB if not specified)
        string ram_str = to_string(ram).substr(0, to_string(ram).find(".") + 3);
        bool ram_valid = ram > 1.0;
        settings_message_map["RAM"] = make_pair(ram_valid, 
            ram_valid ? ram_str + " GB" : 
            "Value " + ram_str + " GB is invalid (must be greater than 1 GB)");

        // Check threads
        size_t max_t = max_threads();
        bool threads_valid = threads > 0 && threads <= max_t;
        settings_message_map["Threads"] = make_pair(threads_valid, 
            threads_valid ? to_string(threads) + " thread(s)" : 
            "Value " + to_string(threads) + " is invalid (must be between 1 and " + to_string(max_t) + ")");

        // Check folder paths
        bool output_valid = !output_folder.empty();
        settings_message_map["Output Folder"] = make_pair(output_valid, 
            output_valid ? output_folder : "Invalid output folder");

        return settings_message_map;
    }

    // Print settings summary and return erroneous properties
    string print_settings() {
        string erroneous_properties;
        auto settings_map = validate_settings();

        for (const auto& [key, value] : settings_map) {
            if (value.first) {
                cout << "[✔] " << key << ": " << value.second << endl;
            } else {
                erroneous_properties += key + " ";
                cout << "[✗] " << key << ": " << value.second << endl;
            }
        }
        return erroneous_properties;
    }
};

#endif