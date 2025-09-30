#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cctype>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "settings.h"
// #include "filters.h"
// #include <cstring>
// #include "sdbg_build.h"
// #include "post_processing.h"
// #include <cctype>
// #include <unordered_map>
// #include "phage_curator.h"
// #include "isolate_protospacers.h"
#ifdef __linux__
#include <sys/sysinfo.h>
#elif defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#ifdef DEBUG
#include "io_ops.h"
#endif

#include "cycle_finder.h"
#include "filters.h"
#include "reads.h"
#include "post_processing.h"
#include "sdbg/sdbg.h"
#include "sdbg_build.h"
#include "settings.h"
#include "spacer_ordering.h"
#include "tmp_utils.h"
#include "evaluation.h"

using namespace std;
namespace fs = std::filesystem;
#ifdef DEBUG
#pragma message("DEBUG is defined")
#else
#pragma message("DEBUG is NOT defined")
#endif
void print_usage(const char* program_name) {
    cout << "-------------------------------------------------------" << endl;
    cout << "\n";
    cout << "mCAAT - Metagenomic CRISPR Array Analysis Tool v. 0.1" << endl;
    cout << "\n";
    cout << "-------------------------------------------------------" << endl;
}

bool check_for_error(Settings& settings) {
    cout << "Step 1. Checking the inputs: " << endl;
    string erroneous_message = settings.print_settings();
    if (erroneous_message.empty()) {
        cout << "All inputs are correct. [✔]" << endl;
        return false;
    }
    cout << "Please check the following: " << erroneous_message << endl;
    return true;
}

// Get total system RAM in GB
double get_total_system_ram() {
    double total_bytes = 0.0;
    
#ifdef __linux__
    struct sysinfo mem_info;
    if (sysinfo(&mem_info) == 0) {
        total_bytes = static_cast<double>(mem_info.totalram) * mem_info.mem_unit;
    }
#elif defined(_WIN32)
    MEMORYSTATUSEX mem_info;
    mem_info.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&mem_info)) {
        total_bytes = static_cast<double>(mem_info.ullTotalPhys);
    }
#elif defined(__APPLE__)
    int mib[2] = {CTL_HW, HW_MEMSIZE};
    int64_t mem_size;
    size_t len = sizeof(mem_size);
    if (sysctl(mib, 2, &mem_size, &len, nullptr, 0) == 0) {
        total_bytes = static_cast<double>(mem_size);
    }
#endif
    
    // Convert bytes to GB
    return total_bytes / (1024.0 * 1024.0 * 1024.0);
}

Settings parse_arguments(int argc, char* argv[]) {
    vector<string> input_files_default;
    Settings settings;

    // Get timestamp for fallback folder naming
    string timestamp = settings.get_timestamp();

    bool output_folder_provided = false;
 bool required_files_provided = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "--help" || arg == "-h" || arg=="") {
            cout << "Usage: ./crispr_analyzer --input_files <file1> [file2] [options]\n"
                 << "\nRequired:\n"
                 << "  --input_files <file1> [file2]   One or two input FASTA/FASTQ files\n"
                 << "\nOptional:\n"
                 << "  --ram <amount>                  RAM to use (e.g., 4G, 500M). Default: 95% of system RAM\n"
                 << "  --threads <num>                 Number of threads. Default: CPU cores - 2\n"
                 << "  --output-folder <path>          Output directory. If not provided, a timestamped folder is created\n"
                 << "  --benchmark <file>              File containing expected crispr sequences line separated\n"
                 << "  --help, -h                      Show this help message\n";
            exit(0);
        }
       
        if (arg == "--input-files" || arg == "-i") {
            while (++i < argc && argv[i][0] != '-') {
                input_files_default.push_back(argv[i]);
                
            }
            --i;
            required_files_provided = true;
        } else if (arg == "--benchmark") {
            if (++i < argc) {
                settings.benchmark_file = argv[i];
            } else {
                throw runtime_error("Error: Missing value for --benchmark");
            }
            --i;
        } else if (arg == "--ram") {
            if (++i < argc) {
                string ram_input = argv[i];
                try {
                    double value = 0.0;
                    char unit = 'G';
                    size_t pos = ram_input.find_first_not_of("0123456789.");
                    if (pos != string::npos) {
                        value = stod(ram_input.substr(0, pos));
                        unit = toupper(ram_input[pos]);
                    } else {
                        value = stod(ram_input);
                    }

                    switch (unit) {
                        case 'B': settings.ram = value / (1024.0 * 1024.0 * 1024.0); break;
                        case 'K': settings.ram = value / (1024.0 * 1024.0); break;
                        case 'M': settings.ram = value / 1024.0; break;
                        case 'G': settings.ram = value; break;
                        default:
                            throw runtime_error("Error: Invalid RAM unit. Use B, K, M, or G.");
                    }

                    double total_ram = get_total_system_ram();
                    if (settings.ram < 1.0) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GB is too low (must be at least 1 GB)");
                    }
                    if (settings.ram > total_ram) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GB exceeds system total of " + to_string(total_ram) + " GB");
                    }
                } catch (const invalid_argument&) {
                    throw runtime_error("Error: Invalid RAM value provided: " + ram_input);
                }
            } else {
                throw runtime_error("Error: Missing value for --ram");
            }
        } else if (arg == "--threads") {
            if (++i < argc) {
                settings.threads = stoul(argv[i]);
                
            } else {
                throw runtime_error("Error: Missing value for --threads");
            }
        } else if (arg == "--output-folder" || arg == "--output_folder") {
            if (++i < argc) {
                settings.output_folder = string(argv[i]);
                output_folder_provided = true;
            } else {
                throw runtime_error("Error: Missing value for --output-folder");
            }
        }
    }
    if (!required_files_provided && input_files_default.empty()) {
        throw runtime_error("Error: No input files provided. Use --input-files <file1> [file2]");
    }
    // Set default output folder if not provided
    if (!output_folder_provided) {
        settings.output_folder = "mcaat_run_" + timestamp;
    }

    // Set subfolders and output file
    settings.graph_folder = settings.output_folder + "/graph";
    settings.cycles_folder = settings.output_folder + "/cycles";
    settings.output_file = settings.output_folder + "/CRISPR_Arrays.txt";

    // Debug output
    cout << "Output folder: " << settings.output_folder << endl;
    cout << "Graph folder: " << settings.graph_folder << endl;
    cout << "Cycles folder: " << settings.cycles_folder << endl;

    // Create directories
    try {
        fs::create_directories(settings.output_folder);
        fs::create_directories(settings.graph_folder);
        fs::create_directories(settings.cycles_folder);
        cout << "Created directories: " << settings.output_folder << ", "
             << settings.graph_folder << ", " << settings.cycles_folder << endl;
    } catch (const fs::filesystem_error& e) {
        throw runtime_error("Error: Could not create directories: " + string(e.what()));
    }

    // Validate input files
    if (input_files_default.size() < 1 || input_files_default.size() > 2) {
        throw runtime_error("Error: You must provide one or two input files.");
    }

    int count = 0;
    for (const auto& file : input_files_default) {
        if (!fs::exists(file)) {
            throw runtime_error("Error: Input file " + file + " does not exist.");
        }
        count++;
        if (count == 1)
            settings.input_files += file;
        if (count == 2)
            settings.input_files += " " + file + " ";
    }

    if (settings.threads == 0) {
        settings.threads = thread::hardware_concurrency() - 2;
    }

    if (settings.ram == 0.0) {
        settings.ram = get_total_system_ram() * 0.95;
    }

    return settings;
}


string fetchNodeLabel(SDBG& sdbg, uint64_t node) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

uint64_t findNodeFromKmer(const SDBG& sdbg, const std::string& kmer) {
    int k = sdbg.k();
    if (static_cast<int>(kmer.size()) != k) {
        std::cerr << "Warning: kmer size " << kmer.size() << " does not match k=" << k << std::endl;
        return UINT64_MAX; // or handle error
    }
    std::vector<uint8_t> seq(k);
    for (int i = 0; i < k; ++i) {
        seq[i] = "ACGT"s.find(kmer[i]) + 1;
    }
    int64_t result = sdbg.IndexBinarySearch(seq.data());
    return (result == -1) ? UINT64_MAX : static_cast<uint64_t>(result);
}

std::map<uint64_t, std::vector<std::vector<uint64_t>>> createRepeatToSpacerNodes(const SDBG& sdbg, const std::map<std::string, std::vector<std::string>>& systems_from_analyzer) {
    std::map<uint64_t, std::vector<std::vector<uint64_t>>> repeat_to_spacer_nodes;
    int k = sdbg.k();
    for (const auto& [repeat, spacers] : systems_from_analyzer) {
        if (static_cast<int>(repeat.size()) < k) continue;
        std::string first_kmer = repeat.substr(0, k);
        uint64_t key_node = findNodeFromKmer(sdbg, first_kmer);
        if (key_node == UINT64_MAX) continue; // invalid
        std::vector<std::vector<uint64_t>> spacer_node_vectors;
        for (const auto& spacer : spacers) {
            std::vector<uint64_t> nodes;
            int L = spacer.size();
            for (int i = 0; i <= L - k; ++i) {
                std::string kmer = spacer.substr(i, k);
                uint64_t node = findNodeFromKmer(sdbg, kmer);
                if (node != UINT64_MAX) {
                    nodes.push_back(node);
                }
            }
            if (!nodes.empty()) {
                spacer_node_vectors.push_back(std::move(nodes));
            }
        }
        if (!spacer_node_vectors.empty()) {
            repeat_to_spacer_nodes[key_node] = std::move(spacer_node_vectors);
        }
    }
    return repeat_to_spacer_nodes;
}
#ifdef DEBUG
int main(int argc, char** argv) {
    // %% PARSE ARGUMENTS %%
    Settings settings = parse_arguments(argc, argv);
    string name_of_genome = "test";
    if (check_for_error(settings)){
        //tell the user which folder we are deleting
        cout<< "Folder " << settings.output_folder << " will be deleted due to errors." << endl;
        
        cout << "Do you want that folder to be removed? (y/n): ";
        char answer;
        cin >> answer;
        if (answer != 'y' && answer != 'Y') {
            cout << "Exiting the program." << endl;
            return 1;
        }
        cout << "Removing folder: " << settings.output_folder << endl;
        fs::remove_all(settings.output_folder); 
        return 1;
    }
    // %% PARSE ARGUMENTS %%

    // %% BUILD GRAPH %%
    //SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
   
    int length_bound = 77;
    SDBG sdbg;
    vector<string> folders = {"/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-27_09-53-11/graph/graph","/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-26-39/graph/graph"};
    string graph_folder_old = settings.graph_folder;///vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-23-46
    settings.graph_folder=folders[1];
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;

    // %% LOAD GRAPH %%
    
    delete[] cstr;
    string ending_seq = "CAGAGATAGAAATTATTTTTATTATACGTTTTTTTGT";
    vector<int64_t> end_nodes;
    //using findNodeFromKmer find all the node ids in the ending_seq
    for (size_t i = 0; i <= ending_seq.size() - sdbg.k(); ++i) {
        string kmer = ending_seq.substr(i, sdbg.k());
        uint64_t node = findNodeFromKmer(sdbg, kmer);
        if (node != UINT64_MAX) {
            end_nodes.push_back(node);
        }
    }
    //print all the end nodes
    cout << "End nodes: ";
    for (const auto& node : end_nodes) {
        cout << node << " ";
    }
    cout << endl;



    // %% FBCE ALGORITHM %%
    cout << "FBCE FROM DEBUG START:" << endl;
    auto start_time = chrono::high_resolution_clock::now();
    CycleFinder cycle_finder(sdbg, length_bound, 27, settings.cycles_folder, settings.threads);
    int number_of_spacers_total = 0;
    auto cycles = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles.size() << endl;
    // %% FBCE ALGORITHM %%
    
    int number_of_spacers = 0;
    
    // %% FILTERS %%
    cout << "FILTERS START:" << endl;
    Filters filters(sdbg, cycles);
    auto  SYSTEMS = filters.ListArrays(number_of_spacers);
    cout<< "Number of spacers: " << number_of_spacers << " before cleaning"<<endl;
    // %% FILTERS %%

    //%% POST PROCESSING %%
    cout << "POST PROCESSING START:" << endl;
    CRISPRAnalyzer analyzer(SYSTEMS, settings.output_file);
    analyzer.run_analysis();
    cout << "Saved in: " << settings.output_file << endl;
    auto systems_from_analyzer = analyzer.getSystems();
    auto repeat_to_spacer_nodes = createRepeatToSpacerNodes(sdbg, systems_from_analyzer);
    cout << "Created repeat_to_spacer_nodes map with " << repeat_to_spacer_nodes.size() << " entries." << endl;
    //%% POST PROCESSING %%

    //%% PROTOSPACER ISOLATION %%
    IsolateProtospacers isolator(sdbg, repeat_to_spacer_nodes);
    pair<std::map<uint64_t,std::set<uint64_t>>,std::map<uint64_t,std::set<uint64_t>>> protospacer_nodes = isolator.getProtospacerNodes();
    auto grouped_paths_protospacers = isolator.DepthLimitedPathsFromInToOut(protospacer_nodes.first, protospacer_nodes.second, 50,1);
    isolator.WritePathsToFile(grouped_paths_protospacers, "grouped_paths_protospacers.txt");
    
    //%% PROTOSPACER ISOLATION %%
    PhageCurator phage_curator(sdbg, grouped_paths_protospacers, cycles);
    auto extended_paths = phage_curator.ExtendFromGroupedPaths(1000, 5000);
    phage_curator.ReconstructPaths(extended_paths);
    phage_curator.WriteSequencesToFasta("PhageCurator.txt");    
    //io_ops::write_nodes_gfa("output.gfa", sdbg);
    

    //%% POST PROCESSING %%
    // %% DELETE THE GRAPH FOLDER %%
    //fs::remove_all(graph_folder_old);
    // %% DELETE THE GRAPH FOLDER %%          
    
}


#else
int main(int argc, char** argv) {
    // %% PARSE ARGUMENTS %%
    Settings settings = parse_arguments(argc, argv);
    if (check_for_error(settings)){
        cout << "Folder " << settings.output_folder;
        cout << " will be deleted due to errors." << endl;
        
        cout << "Do you want that folder to be removed? (y/n): ";
        char answer;
        cin >> answer;
        if (answer != 'y' && answer != 'Y') {
            cout << "Exiting the program." << endl;
            return 1;
        }
        cout << "Removing folder: " << settings.output_folder << endl;
        fs::remove_all(settings.output_folder); 
        return 1;
    }
    // %% PARSE ARGUMENTS %%

    // %% BUILD GRAPH %%
    SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
    // %% LOAD GRAPH %%
    const int length_bound = 77;
    string graph_path = settings.graph_folder + "/graph";
    cout << "Graph folder: " << graph_path << endl;

    SDBG sdbg;
    sdbg.LoadFromFile(graph_path.c_str());
    cout << "Loaded the graph" << endl;
    // %% LOAD GRAPH %%

    // %% FBCE ALGORITHM %%
    cout << "FBCE START:" << endl;
    CycleFinder cycle_finder(sdbg, length_bound, 27, settings.cycles_folder, settings.threads);
    auto cycles_map = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles_map.size() << endl;
    // %% FBCE ALGORITHM %%

    auto cycles = cycles_map_to_cycles(cycles_map);
    
    // // %% FILTER CYCLES %%
    // cout << "FILTER CYCLES START:" << endl;
    // int amount_of_cycles_before = get_cycle_count(cycles_map);
    // // keep_relevant_cycles(cycles_map);
    // int amount_of_cycles_after = get_cycle_count(cycles_map);
    // cout << amount_of_cycles_after << " out of ";
    // cout << amount_of_cycles_before << " are kept" << endl;
    // // %% FILTER CYCLES %%

    std::chrono::_V2::system_clock::time_point start_time;
    std::chrono::_V2::system_clock::time_point end_time;

    cout << "\n══════════════════════════════════════════════" << endl;
    cout << "🔸STEP 6: Finding relevant reads" << endl;
    cout << "══════════════════════════════════════════════" << endl;
    start_time = std::chrono::high_resolution_clock::now();

    auto fastq_files = get_fastq_files_from_settings(settings);
    auto reads = get_reads(
        sdbg,
        fastq_files.first,
        fastq_files.second,
        cycles
    );
    cout << "    ▸ Found " << reads.size() << " reads" << endl;
    if (reads.size() == 0) {
        cout << "══════════════════════════════════════════════" << endl;
        return 0;
    }

    end_time = std::chrono::high_resolution_clock::now();
    cout << "\n⏳ Time elapsed: ";
    cout << std::fixed << std::setprecision(2);
    cout << std::chrono::duration<double>(end_time - start_time).count();
    cout << " seconds" << endl;

    cout << "\n══════════════════════════════════════════════" << endl;
    cout << "🔸STEP 7: Order the spacers" << endl;
    cout << "══════════════════════════════════════════════" << endl;
    start_time = std::chrono::high_resolution_clock::now();

    
    cout << "  ▸ Splitting into subproblems" << endl;
    const size_t reads_size = reads.at(0).size();
    auto subgraphs = get_crispr_regions_extended_by_k(sdbg, reads_size, cycles);

    cout << "  🔄 Filtering subproblems:" << endl;
    vector<Graph> remaining_subgraphs;
    vector<vector<vector<uint64_t>>> remaining_reads;
    vector<vector<vector<size_t>>> remaining_cycles;
    for (size_t idx = 0; idx < subgraphs.size(); ++idx) {
        const auto& subgraph = subgraphs[idx];
        const auto relevant_reads = get_relevant_reads(subgraph, reads);
        auto relevant_cycles = get_relevant_cycles(subgraph, cycles);

        get_minimum_cycles_for_full_coverage(relevant_cycles);
        
        // The assembly of megahit always assembles the graph in its reverse complement.
        // We discard the reverse complement by assuming that it won't have any relevant reads
        if (relevant_reads.size() == 0 || relevant_cycles.size() < 3) {
            continue;
        }

        remaining_subgraphs.push_back(subgraph);
        remaining_reads.push_back(relevant_reads);
        remaining_cycles.push_back(relevant_cycles);
    }
    cout << "  ✅ Filtered out " << subgraphs.size() - remaining_subgraphs.size();
    cout << "/" << subgraphs.size() << " subproblems" << endl;

    
    cout << "  🔄 Solving " << remaining_subgraphs.size();
    cout << " subproblems..." << endl;
    vector<tuple<string, string, vector<string>, float>> found_systems;
    for (size_t idx = 0; idx < remaining_subgraphs.size(); ++idx) {
        const auto& subgraph = remaining_subgraphs[idx];
        const auto& relevant_reads = remaining_reads[idx];
        const auto& relevant_cycles = remaining_cycles[idx];
        
        cout << "    Subproblem " << idx + 1 << "/";
        cout << remaining_subgraphs.size() << ":" << endl;
        
        cout << "      🛈 Graph with " << subgraph.nodes.size();
        cout << " nodes and " << subgraph.edge_count() << " edges" << endl;

        cout << "      🛈 Reads with " << relevant_reads.size() << "/";
        cout << reads.size() << " used" << endl;

        cout << "      🛈 Cycles with " << relevant_cycles.size() << "/";
        cout << get_cycle_count(cycles_map) << " used" << endl;

        float confidence = 1.0;
        auto cycle_order = order_cycles(subgraph, relevant_reads, relevant_cycles, confidence);

        cout << "      ▸ The order is ";
        for (auto node : cycle_order) {
            cout << node << " ";
        }
        cout << endl;

        cout << "      ▸ Turning the cycle order into a node order" << endl;
        auto ordered_cycles = get_ordered_cycles(cycle_order, relevant_cycles);
        
        if (ordered_cycles.size() < 2) {
            cout << "      ▸ Node order is to short and is not processed further" << endl;
            continue;
        }

        cout << "      ▸ Starting the filter process:" << endl;
        auto [repeat, spacers, full_sequence] = get_systems(sdbg, ordered_cycles);
        
        cout << "        ▸ Number of spacers: " << spacers.size() << endl;

        found_systems.push_back(std::make_tuple(full_sequence, repeat, spacers, confidence));
    }
    cout << "  ✅ Completed each subproblem" << endl;

    end_time = std::chrono::high_resolution_clock::now();
    cout << "\n⏳ Time elapsed: ";
    cout << std::fixed << std::setprecision(2);
    cout << std::chrono::duration<double>(end_time - start_time).count();
    cout << " seconds" << endl;


    if (settings.benchmark_file != "") {
        cout << "\n══════════════════════════════════════════════" << endl;
        cout << "🔸STEP 8: Compare to ground of truth using benchmark file" << endl;
        cout << "══════════════════════════════════════════════" << endl;
        start_time = std::chrono::high_resolution_clock::now();
    
        vector<string> benchmark_sequences;
        ifstream benchmark_in(settings.benchmark_file);
        if (!benchmark_in) {
            cerr << "Error: Could not open benchmark file: " << settings.benchmark_file << endl;
        } else {
            string line;
            while (getline(benchmark_in, line)) {
                if (!line.empty()) {
                    benchmark_sequences.push_back(line);
                }
            }
            benchmark_in.close();
            cout << "Loaded " << benchmark_sequences.size() << " benchmark sequences." << endl;
        }

        cout << "  ▸ " << found_systems.size();
        cout << " crispr sequences are found and benchmarked using ";
        cout << benchmark_sequences.size() << " sequences" << endl;

        size_t no_match_count = 0;
        float average_sequence_similarity = 0.0;
        for (const auto& [sequence, repeat, spacers, confidence] : found_systems) {
            // This leads to overestimation by choosing greedy
            const auto expected_sequence = get_most_similar_sequence(
                sequence, benchmark_sequences
            );
            if (expected_sequence == "") {
                cout << "    ▸ No expected match for sequence: ";
                cout << sequence << endl;
                no_match_count++;
                continue;
            }
            
            const auto sequence_similarity = get_string_similarity(
                sequence, expected_sequence
            );
            vector<float> spacer_similarity = {0.0, 0.0, 0.0};
            for (int variant = 0; variant < 3; ++variant) {
                spacer_similarity[variant] = get_spacer_order_similarity(
                    spacers, expected_sequence, variant
                );
            }

            const auto amount_of_duplicate_spacers = get_number_of_duplicate_spacers(
                spacers, expected_sequence
            );

            cout << "    ▸ ≥" << std::fixed << std::setprecision(2);
            cout << (sequence_similarity * 100) << "% sequence similarity, ";
            cout << (spacer_similarity[0] * 100) << "% spacer similarity variant 0, ";
            cout << (spacer_similarity[1] * 100) << "% spacer similarity variant 1, ";
            cout << (spacer_similarity[2] * 100) << "% spacer similarity variant 2, ";
            cout << "with ";
            cout << spacers.size() << " spacers, ";
            cout << amount_of_duplicate_spacers << " duplicate spacers, confidence of cycle resolution: ";
            cout << (confidence * 100) << "%, and the repeat: ";
            cout << repeat << ", and sequence: ";
            cout << sequence << endl;
            average_sequence_similarity += sequence_similarity;
        }
        average_sequence_similarity /= static_cast<float>(found_systems.size() - no_match_count);

        cout << "  ▸ The average sequence similarity is ";
        cout << std::fixed << std::setprecision(2);
        cout << (average_sequence_similarity * 100);
        cout << "% with " << no_match_count << "/" << found_systems.size();
        cout << " ignored" << endl;

        end_time = std::chrono::high_resolution_clock::now();
        cout << "\n⏳ Time elapsed: ";
        cout << std::fixed << std::setprecision(2);
        cout << std::chrono::duration<double>(end_time - start_time).count();
        cout << " seconds" << endl;
    }
    
    cout << "══════════════════════════════════════════════" << endl;

    // %% POST PROCESSING %%
    cout << "POST PROCESSING START:" << endl;
    unordered_map<string, vector<string>> all_systems;
    for (const auto& [_sequence, repeat, spacers, _confidence] : found_systems) {
        all_systems[repeat] = spacers;
    }
    CRISPRAnalyzer analyzer(all_systems, settings.output_file);
    analyzer.run_analysis();
    cout << "Saved in: " << settings.output_file << endl;
    // %% POST PROCESSING %%

    // %% DELETE THE GRAPH FOLDER %%
    try {
        fs::remove_all(settings.graph_folder);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Warning: Could not remove graph folder: " << e.what() << std::endl;
    }
    // %% DELETE THE GRAPH FOLDER %%            
}
#endif