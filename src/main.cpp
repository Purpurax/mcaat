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
#ifdef __linux__
#include <sys/sysinfo.h>
#elif defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#ifdef DEVELOP
#include "io_ops.h"
#endif

#include "cycle_finder.h"
#include "cycle_filter.h"
#include "filters.h"
#include "jumps.h"
#include "post_processing.h"
#include "sdbg/sdbg.h"
#include "sdbg_build.h"
#include "settings.h"
#include "spacer_ordering.h"
#include "tmp_utils.h"
#include "evaluation.h"

using namespace std;
namespace fs = std::filesystem;

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
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}
// int main(int argc, char** argv) {
//     // %% PARSE ARGUMENTS %%
//     Settings settings = parse_arguments(argc, argv);
//     string name_of_genome = "test";
//     if (check_for_error(settings)){
//         //tell the user which folder we are deleting
//         cout<< "Folder " << settings.output_folder << " will be deleted due to errors." << endl;
        
//         cout << "Do you want that folder to be removed? (y/n): ";
//         char answer;
//         cin >> answer;
//         if (answer != 'y' && answer != 'Y') {
//             cout << "Exiting the program." << endl;
//             return 1;
//         }
//         cout << "Removing folder: " << settings.output_folder << endl;
//         fs::remove_all(settings.output_folder); 
//         return 1;
//     }
//     // %% PARSE ARGUMENTS %%

//     // %% BUILD GRAPH %%
//      SDBGBuild sdbg_build(settings);
//     // %% BUILD GRAPH %%
    
    
//     int length_bound = 77;
//     SDBG sdbg;
//     string graph_folder_old = settings.graph_folder;
//     settings.graph_folder= settings.graph_folder+ "/graph"; //"/vol/d/development/git/mcaat_master/mcaat/build/mcaat_run_2025-08-27_09-53-11/graph/graph";
//     char * cstr = new char [settings.graph_folder.length()+1];
//     std::strcpy (cstr, settings.graph_folder.c_str());
//     cout << "Graph folder: " << cstr << endl;
//     sdbg.LoadFromFile(cstr);
//     cout << "Loaded the graph" << endl;
//     /*
//     string ending_kmer = "ATTTTTATTATACGTTTTTTTGT";
//     uint8_t end_seq[24];
//     for (int i = 0; i < 23; ++i) {
//         end_seq[i] = "ACGT"s.find(ending_kmer[i]) + 1;
//     }
//     int64_t end_node = sdbg.IndexBinarySearch(end_seq);
//     cout<<"EdgeMultiplicity: "<<sdbg.EdgeMultiplicity(end_node)<<endl;
//     cout<<"Indegree: "<<sdbg.EdgeIndegree(end_node)<<endl;
//     cout<<"Outdegree: "<<sdbg.EdgeOutdegree(end_node)<<endl;
//     // %% LOAD GRAPH %%
//     cout << "Loaded k-mer: " << ending_kmer << " as node: " << end_node << endl;
//     PhageCurator phage_curator(sdbg);
//     std::vector<std::vector<uint64_t>> paths = phage_curator.DepthLimitedPaths(end_node, 1000,5000);
//     phage_curator.ReconstructPaths(paths);
//     for(const auto& sequence : phage_curator.reconstructed_sequences) {
//         std::cout << sequence << std::endl;
//     }
// */
//     delete[] cstr;

    
//     // %% FBCE ALGORITHM %%
//     cout << "FBCE START:" << endl;
//     auto start_time = chrono::high_resolution_clock::now();
//     CycleFinder cycle_finder(sdbg, length_bound, 27, settings.cycles_folder, settings.threads);
    
//     int number_of_spacers_total = 0;
//     auto cycles = cycle_finder.results;
//     cout << "Number of nodes in results: " << cycles.size() << endl;
//     // %% FBCE ALGORITHM %%
    
//     // ############ DEVELOPMENT ############

//     #ifdef DEVELOP
//     //io_ops::read_cycles("cycles.json");
//     //io_ops::write_cycles("out.json", cycles);
//     //io_ops::write_nodes_gfa("out.gfa", sdbg); // CAREFUL: WILL TAKE A LOT OF SPACE IF GRAPH IS HUGE!
//     #endif

//     // ############ DEVELOPMENT ############


//     int number_of_spacers = 0;
//     // // %% FILTERS %%
//      cout << "FILTERS START:" << endl;
//      Filters filters(sdbg, cycles);
//      auto  SYSTEMS = filters.ListArrays(number_of_spacers);
//      cout<< "Number of spacers: " << number_of_spacers << " before cleaning"<<endl;
//     // // %% FILTERS %%



//     // //%% POST PROCESSING %%
//      cout << "POST PROCESSING START:" << endl;
//      CRISPRAnalyzer analyzer(SYSTEMS, settings.output_file);
//      analyzer.run_analysis();
//      cout << "Saved in: " << settings.output_file << endl;
//     //%% POST PROCESSING %%

//     // %% DELETE THE GRAPH FOLDER %%
// }

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
    settings.graph_folder=folders[0];
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;

    // %% LOAD GRAPH %%
    
    delete[] cstr;

        
    string ending_kmer = "ATTTTTATTATACGTTTTTTTGT";
    uint8_t end_seq[24];
    for (int i = 0; i < 23; ++i) {
        end_seq[i] = "ACGT"s.find(ending_kmer[i]) + 1;
    }
    int64_t end_node = sdbg.IndexBinarySearch(end_seq);
    cout<<"EdgeMultiplicity: "<<sdbg.EdgeMultiplicity(end_node)<<endl;
    cout<<"Indegree: "<<sdbg.EdgeIndegree(end_node)<<endl;
    cout<<"Outdegree: "<<sdbg.EdgeOutdegree(end_node)<<endl;
    // %% LOAD GRAPH %%
    cout << "Loaded k-mer: " << ending_kmer << " as node: " << end_node << endl;
    PhageCurator phage_curator(sdbg);
    std::vector<std::vector<uint64_t>> paths = phage_curator.DepthLimitedPaths(end_node, 1000,5000);
    phage_curator.ReconstructPaths(paths);
    for(const auto& sequence : phage_curator.reconstructed_sequences) {
        std::cout << sequence << std::endl;
    }
/*
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
    //%% POST PROCESSING %%

    // %% DELETE THE GRAPH FOLDER %%
    //fs::remove_all(graph_folder_old);
    // %% DELETE THE GRAPH FOLDER %%          
    */  
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
    
    // %% FILTER CYCLES %%
    cout << "FILTER CYCLES START:" << endl;
    int amount_of_cycles_before = get_cycle_count(cycles_map);
    // keep_relevant_cycles(cycles_map);
    int amount_of_cycles_after = get_cycle_count(cycles_map);
    cout << amount_of_cycles_after << " out of ";
    cout << amount_of_cycles_before << " are kept" << endl;
    // %% FILTER CYCLES %%

    std::chrono::_V2::system_clock::time_point start_time;
    std::chrono::_V2::system_clock::time_point end_time;

    cout << "\n══════════════════════════════════════════════" << endl;
    cout << "🔸STEP 6: Creating the jumps" << endl;
    cout << "══════════════════════════════════════════════" << endl;
    start_time = std::chrono::high_resolution_clock::now();

    auto fastq_files = get_fastq_files_from_settings(settings);
    auto jumps = get_jumps_from_reads(sdbg, fastq_files.first, fastq_files.second, settings.threads);
    cout << "    ▸ Created " << jumps.size() << " jumps" << endl;

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
    auto subgraphs = get_crispr_regions(sdbg, cycles_map);

    cout << "  🔄 Filtering subproblems:" << endl;
    vector<Graph> remaining_subgraphs;
    vector<vector<Jump>> remaining_jumps;
    vector<vector<vector<size_t>>> remaining_cycles;
    for (size_t idx = 0; idx < subgraphs.size(); ++idx) {
        const auto& subgraph = subgraphs[idx];
        const auto relevant_jumps = get_relevant_jumps(subgraph, jumps);
        auto relevant_cycles = get_relevant_cycles(subgraph, cycles_map);

        get_minimum_cycles_for_full_coverage(relevant_cycles);
        
        // The assembly of megahit always assembles the graph in its reverse complement.
        // We discard the reverse complement by assuming that it won't have any relevant jumps
        if (relevant_jumps.size() == 0 || relevant_cycles.size() < 3) {
            continue;
        }

        remaining_subgraphs.push_back(subgraph);
        remaining_jumps.push_back(relevant_jumps);
        remaining_cycles.push_back(relevant_cycles);
    }
    cout << "  ✅ Filtered out " << subgraphs.size()-remaining_subgraphs.size();
    cout << "/" << subgraphs.size() << " subproblems" << endl;

    
    cout << "  🔄 Solving " << remaining_subgraphs.size();
    cout << " subproblems..." << endl;
    vector<tuple<string, string, vector<string>>> found_systems;
    for (size_t idx = 0; idx < remaining_subgraphs.size(); ++idx) {
        const auto& subgraph = remaining_subgraphs[idx];
        const auto& relevant_jumps = remaining_jumps[idx];
        const auto& relevant_cycles = remaining_cycles[idx];
        
        cout << "    Subproblem " << idx + 1 << "/";
        cout << remaining_subgraphs.size() << ":" << endl;
        
        cout << "      🛈 Graph with " << subgraph.nodes.size();
        cout << " nodes and " << subgraph.edge_count() << " edges" << endl;

        cout << "      🛈 Jumps with " << relevant_jumps.size() << "/";
        cout << jumps.size() << " used" << endl;

        cout << "      🛈 Cycles with " << relevant_cycles.size() << "/";
        cout << get_cycle_count(cycles_map) << " used" << endl;

        auto cycle_order = order_cycles(subgraph, relevant_jumps, relevant_cycles);

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

        if (full_sequence.size() >= 23) {
            found_systems.push_back(std::make_tuple(full_sequence, repeat, spacers));
        } else {
            cout << "        ▸ The sequence is discarded as it is shorter than 23" << endl;
        }
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
        float average_similarity = 0.0;
        for (const auto& [sequence, _repeat, spacers] : found_systems) {
            const auto similarity = compare_sequence_to_ground_of_truth(
                sequence,
                benchmark_sequences
            );

            if (similarity == -1.0) {
                cout << "    ▸ No expected match for sequence: ";
                cout << sequence << endl;
                no_match_count++;
            } else {
                cout << "    ▸ ≥" << std::fixed << std::setprecision(2);
                cout << (similarity * 100) << "% similarity with ";
                cout << spacers.size() << " spacers and sequence: ";
                cout << sequence << endl;
                average_similarity += similarity;
            }
        }
        average_similarity /= static_cast<float>(found_systems.size() - no_match_count);

        cout << "  ▸ The average similarity is ";
        cout << std::fixed << std::setprecision(2);
        cout << (average_similarity * 100);
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
    for (const auto& [_sequence, repeat, spacers] : found_systems) {
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