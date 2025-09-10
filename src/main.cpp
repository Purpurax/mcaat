//test myLib.h
#include <iostream>
#include "sdbg/sdbg.h"
#include "cycle_finder.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "filters.h"
#include <cstring>
#include "sdbg_build.h"
#include "post_processing.h"
#include <cctype>
#include <unordered_map>
#include "phage_curator.h"
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
                 << "  --help, -h                      Show this help message\n";
            exit(0);
        }
       
        if (arg == "--input-files" || arg == "-i") {
            while (++i < argc && argv[i][0] != '-') {
                input_files_default.push_back(argv[i]);
                
            }
            --i;
            required_files_provided = true;
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

#ifdef RELEASE
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
    SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
   
    int length_bound = 77;
    SDBG sdbg;
    string graph_folder_old = settings.graph_folder;
    settings.graph_folder+="/graph";
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;

    // %% LOAD GRAPH %%
    
    delete[] cstr;

    
    // %% FBCE ALGORITHM %%
    cout << "FBCE START:" << endl;
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
    fs::remove_all(graph_folder_old);
    // %% DELETE THE GRAPH FOLDER %%            
}
#endif


#ifdef DEVELOP
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
    vector<string> folders = {"/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-23-46/graph/graph","/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-26-39/graph/graph"};
    string graph_folder_old = settings.graph_folder;///vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-23-46
    settings.graph_folder=folders[1];
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;

    // %% LOAD GRAPH %%
    
    delete[] cstr;

    
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
}
#endif