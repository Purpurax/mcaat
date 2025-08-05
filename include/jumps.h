/**
 * @file jumps.h
 * @brief Utilities for mapping sequencing reads to graph nodes using k-mer IDs.
 *
 * This header provides structures and functions for representing and extracting
 * simplified read mappings ("jumps") onto a graph, primarily for use in graph-based
 * sequence analysis. It includes functionality for parsing FASTQ files, manipulating
 * nucleotide sequences, and interfacing with succinct De Bruijn graphs (SDBG).
 *
 * Main features:
 * - Jump struct for abstracting read paths via node IDs.
 * - Functions for extracting sequences from FASTQ files.
 * - Utilities for mapping k-mers to graph node IDs.
 * - Methods for constructing jumps from reads and graph mappings.
 */

#ifndef INCLUDE_JUMPS_H_
#define INCLUDE_JUMPS_H_

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <optional>

#include "sdbg/sdbg.h"
#include "sequence/io/sequence_lib.h"
#include "idba/sequence.h"
#include "kseq++/seqio.hpp"

using namespace std;

/**
 * @struct Jump
 * @brief Represents a simplified mapping of a read onto a graph using node IDs.
 *
 * The Jump structure abstracts a read by storing only the start and end node IDs (k-mer IDs)
 * and the number of nodes traversed in between, omitting the actual sequence data.
 * This is useful for graph-based sequence analysis where only the path information is required.
 *
 * @note The Jump depends on a graph structure that assigns unique node IDs to k-mers.
 */
struct Jump {
    uint64_t start_k_mer_id;
    uint64_t end_k_mer_id;
    uint64_t nodes_in_between;
};

/**
 * @internal
 * @brief Get the kmer to node map using the sdbg graph
 * 
 * @param sdbg The Succinct De Bruijn Graph
 * @return A map from a kmer (String) to a node id (uint64_t)
 */
unordered_map<string, uint64_t> get_kmer_to_node_map(const SDBG& sdbg);

/**
 * @brief Reads all the sequences in the fastq file
 * 
 * Uses kseqpp to parse the fastq files and filters out only the sequences
 * 
 * @param fastq_file_path
 * @return sequences (vector<string>)
 */
vector<string> extract_sequences_from_fastq_file(const string& fastq_file_path);

/**
 * @internal
 * @brief Reverses the string and flips the nucleotides
 * 
 * The String is reversed and every 'A' - 'T' and 'C' - 'G' are flipped.
 * Any other character is not modified
 * 
 * @param sequence 
 * @return reversed sequence (string)
 */
string reverse_pair_ends_sequence(string sequence);

/**
 * @internal INTERNAL STUFF
 * @brief Try to create a jump from a sequence using the kmer to node map
 * 
 * @param sequence The sequence from which the jump shall be constructed
 * @param kmer_to_node_id The map from kmers (K-sized strings) to node_ids
 * @param K The kmer size
 * @return Returns the Jump as (std::optional<Jump>) and nothing if Jump is not valid 
 */
optional<Jump> create_jump_from_sequence(
    const string& sequence,
    const unordered_map<string, uint64_t>& kmer_to_node_id,
    const uint32_t K
);

/**
 * @brief Get the sequences from the fastq files and turns them into jumps
 * 
 * @note Reads that are shorter 2 * K (with K = sdbg.k()) are omitted
 * @note Reads where either the start sequence or the end sequence cannot be mapped to a node id, are omitted
 * 
 * @pre fastq files provided have to be valid
 * 
 * @param sdbg The graph used to map sequence k-mer to node id
 * @param fastq_file_1 Normal fastq read file
 * @param fastq_file_2 Optional second fastq file with reversed complement sequences
 * @param thread_count The amount of threads used to construct the jumps parallely
 * 
 * @return All valid Jumps that can be build using the read files
 */
vector<Jump> get_jumps_from_reads(
    const SDBG& sdbg,
    const string& fastq_file_1,
    const optional<string>& fastq_file_2,
    const size_t thread_count
);

#endif
