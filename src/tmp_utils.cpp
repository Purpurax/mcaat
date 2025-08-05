#include "tmp_utils.h"

void trim_string(string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r"));
    s.erase(s.find_last_not_of(" \t\n\r") + 1);
}

pair<string, optional<string>> get_fastq_files_from_settings(
    const Settings& settings
) {
    const bool two_fastq_files_provided = settings.input_files.find(" ") != string::npos;
    if (two_fastq_files_provided) {
        const size_t space_pos = settings.input_files.find(" ");
        string input_file1 = settings.input_files.substr(0, space_pos);
        string input_file2 = settings.input_files.substr(space_pos + 1);

        trim_string(input_file1);
        trim_string(input_file2);

        return std::make_pair(input_file1, optional<string>(input_file2));
    } else {
        return std::make_pair(settings.input_files, std::nullopt);
    }
}

size_t get_cycle_count(
    const unordered_map<uint64_t, vector<vector<uint64_t>>>& cycles_map
) {
    size_t count = 0;
    for (const auto& pair : cycles_map) {
        count += pair.second.size();
    }
    return count;
}
