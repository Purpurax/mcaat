#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <rapidfuzz/fuzz.hpp>

class CRISPRAnalyzer {
private:
    std::string output_path;
    std::unordered_map<std::string, std::vector<std::string>> systems;
    int omitted_repeats = 0;
    int total_spacers = 0;
    int amount, min_sl, max_sl, min_rl, max_rl, mean_similarity;
    std::map<std::string, std::vector<std::string>> grouped_repeat_cycles;

public:
    CRISPRAnalyzer(unordered_map<string, vector<string>> systems_map,
                std::string output = "crispr_report.txt",
                int amt = 2, int minsl = 23, int maxsl = 50,
                int minrl = 23, int maxrl = 50, int mean_sim = 90)
        : systems(std::move(systems_map)), output_path(std::move(output)), amount(amt),
        min_sl(minsl), max_sl(maxsl), min_rl(minrl), max_rl(maxrl),
        mean_similarity(mean_sim) {}

    std::map<std::string, std::vector<std::string>> getSystems() const {
        return this->grouped_repeat_cycles;
    }

    void parse_input(const std::string& content) {
        std::istringstream stream(content);
        std::string line, repeat;
        while (std::getline(stream, line)) {
            if (line.empty() || line == "----------------------------------") continue;
            if (line.find("Repeat:") == 0) {
                repeat = line.substr(7);
                repeat.erase(0, repeat.find_first_not_of(" \t"));
                systems[repeat] = {};
            } else if (line.find("Number of Spacers:") == std::string::npos && line != "Spacers:") {
                systems[repeat].push_back(line);
            }
        }
    }
    std::vector<std::string> get_common_kmers(const std::vector<std::string>& kmers,
                                              const std::vector<std::string>& sequences) {
        std::unordered_map<std::string, int> count;
        for (const auto& kmer : kmers) {
            count[kmer]++;
        }

        std::vector<std::string> common;
        int threshold = sequences.size() * 0.75;
        for (const auto& pair : count) {
            if (pair.second >= threshold) {
                common.push_back(pair.first);
            }
        }
        return common;
    }

    std::vector<std::string> find_common_prefix_kmers(const std::vector<std::string>& sequences, int k) {
        std::vector<std::string> kmers;
        for (const auto& seq : sequences) {
            for (int i = 1; i <= std::min(k, (int)seq.size()); ++i) {
                kmers.push_back(seq.substr(0, i));
            }
        }
        return get_common_kmers(kmers, sequences);
    }

    std::vector<std::string> find_common_suffix_kmers(const std::vector<std::string>& sequences, int k) {
        std::vector<std::string> kmers;
        for (const auto& seq : sequences) {
            for (int i = std::max(0, (int)seq.size() - k); i < (int)seq.size(); ++i) {
                kmers.push_back(seq.substr(i));
            }
        }
        return get_common_kmers(kmers, sequences);
    }

    std::vector<std::string> trim_kmers_from_sequences(const std::vector<std::string>& sequences,
                                                       const std::vector<std::string>& prefixes,
                                                       const std::vector<std::string>& suffixes) {
        std::vector<std::string> trimmed;
        for (auto seq : sequences) {
            for (const auto& pre : prefixes) {
                if (seq.find(pre) == 0) {
                    seq = seq.substr(pre.size());
                    break;
                }
            }
            for (const auto& suf : suffixes) {
                if (seq.size() >= suf.size() &&
                    seq.compare(seq.size() - suf.size(), suf.size(), suf) == 0) {
                    seq = seq.substr(0, seq.size() - suf.size());
                    break;
                }
            }
            if ((int)seq.size() >= min_sl && (int)seq.size() <= max_sl) {
                trimmed.push_back(seq);
            }
        }
        return trimmed;
    }
    bool validate_spacer_diversity(const std::vector<std::string>& sequences) {
        std::vector<double> scores;
        for (size_t i = 0; i < sequences.size(); ++i) {
            for (size_t j = i + 1; j < sequences.size(); ++j) {
                double score = rapidfuzz::fuzz::ratio(sequences[i], sequences[j]);
                scores.push_back(score);
            }
        }
        if (scores.empty()) return false;
        double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
        double mean = sum / scores.size();
        return mean <= mean_similarity;
    }

    std::vector<std::string> filter_substring_spacers(const std::vector<std::string>& spacers) {
        std::vector<std::string> filtered;
        std::unordered_set<std::string> kept;
        // Sort by length descending to check shorter ones against longer ones
        std::vector<std::string> sorted = spacers;
        std::sort(sorted.begin(), sorted.end(), [](const std::string& a, const std::string& b) {
            return a.size() > b.size();
        });
        for (const auto& spacer : sorted) {
            bool is_substring = false;
            for (const auto& kept_spacer : kept) {
                if (rapidfuzz::fuzz::partial_ratio(spacer, kept_spacer) >= 90.0) {
                    is_substring = true;
                    break;
                }
            }
            if (!is_substring) {
                kept.insert(spacer);
                filtered.push_back(spacer);
            }
        }
        return filtered;
    }

    std::vector<std::string> filter_by_length(const std::vector<std::string>& spacers) {
        std::vector<std::string> filtered;
        for (const auto& spacer : spacers) {
            if ((int)spacer.size() >= min_sl && (int)spacer.size() <= max_sl) {
                filtered.push_back(spacer);
            }
        }
        return filtered;
    }

    std::string reconstruct_repeat(const std::string& original,
                                   const std::vector<std::string>& prefixes,
                                   const std::vector<std::string>& suffixes) {
        std::string result = original;
        if (!prefixes.empty()) result += prefixes.back();
        if (!suffixes.empty()) result = suffixes.front() + result;
        return result;
    }

    void generate_report(const std::string& repeat,
                         const std::vector<std::string>& spacers,
                         std::ofstream& out) {
        out << "--------------------------------------------------\n";
        out << repeat << "\n";
        this->grouped_repeat_cycles[repeat] = {};
        out << "--------------------------------------------------\n";
        for (const auto& spacer : spacers) {
            out << spacer << "\n";
            this->grouped_repeat_cycles[repeat].push_back(spacer);
        }
        out << "--------------------------------------------------\n";
        out << "Number of Spacers: " << spacers.size() << "\n";
        out << "--------------------------------------------------\n\n";
    }
    void run_analysis() {

        std::ofstream report(output_path);
        report << "CRISPR Analysis Report\n";
        report << "The tool was run with the following parameters:\n";
        report << "Amount of Spacers: " << amount << "\n";
        report << "[Min:Max] Length of Spacers: [" << min_sl << ":" << max_sl << "]\n";
        report << "[Min:Max] Length of Repeats: [" << min_rl << ":" << max_rl << "]\n";
        report << "Mean Similarity Between Spacers: " << mean_similarity << "\n";
        report << "Conservation Threshold: 80%\n";
        report << "--------------------------------------------------\n";
        for (const auto& pair : systems) {
            const std::string& repeat = pair.first;
            const std::vector<std::string>& spacers = pair.second;

            if (spacers.size() < 2) {
                omitted_repeats++;
                continue;
            }

            int k = this->max_rl - repeat.size();
            auto prefix_kmers = find_common_prefix_kmers(spacers, k);
            auto suffix_kmers = find_common_suffix_kmers(spacers, k);
            std::string updated_repeat = reconstruct_repeat(repeat, prefix_kmers, suffix_kmers);

            if ((int)updated_repeat.size() < min_rl || (int)updated_repeat.size() > max_rl) {
                omitted_repeats++;
                continue;
            }

            auto trimmed = trim_kmers_from_sequences(spacers, prefix_kmers, suffix_kmers);
            if ((int)trimmed.size() < amount) {
                omitted_repeats++;
                continue;
            }

            std::unordered_set<std::string> unique(trimmed.begin(), trimmed.end());
            std::vector<std::string> unique_vec(unique.begin(), unique.end());

            unique_vec = filter_substring_spacers(unique_vec);

            unique_vec = filter_by_length(unique_vec);

            if ((int)unique_vec.size() < amount) {
                omitted_repeats++;
                continue;
            }

            // Recompute kmers after filtering substrings
            auto new_prefix_kmers = find_common_prefix_kmers(unique_vec, k);
            auto new_suffix_kmers = find_common_suffix_kmers(unique_vec, k);
            updated_repeat = reconstruct_repeat(repeat, new_prefix_kmers, new_suffix_kmers);

            if ((int)updated_repeat.size() < min_rl || (int)updated_repeat.size() > max_rl) {
                omitted_repeats++;
                continue;
            }

            // Trim again with new kmers
            unique_vec = trim_kmers_from_sequences(unique_vec, new_prefix_kmers, new_suffix_kmers);

            if ((int)unique_vec.size() < amount) {
                omitted_repeats++;
                continue;
            }

            if (!validate_spacer_diversity(unique_vec)) {
                omitted_repeats++;
                continue;
            }

            generate_report(updated_repeat, unique_vec, report);
            total_spacers += unique_vec.size();
        }
        //fill the data in for the grouped_repeat_cycles


        report << "Number of Systems: " << (systems.size() - omitted_repeats) << "\n";
        report << "Number of Spacers: " << total_spacers << "\n";
        report << "Omitted Repeats: " << omitted_repeats << "\n";
    }
};