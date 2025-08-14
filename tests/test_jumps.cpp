#include <gtest/gtest.h>
#include "jumps.h"

#include "construct_test_data.h"
#include <string>
#include <unordered_map>

using namespace std;


TEST(JumpsTest, ReversePairEndsSequenceBasic) {
    const string input = "ACGT";
    const string expected = "ACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceBasic2) {
    const string input = "ACGTT";
    const string expected = "AACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceMixed) {
    const string input = "AAGCT";
    const string expected = "AGCTT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceNonStandardChars) {
    const string input = "AXGT";
    const string expected = "ACXT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceNonStandardChars2) {
    const string input = "XXXQUIOCPOPYM";
    const string expected = "MYPOPGOIUQXXX";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, CreateJumpFromSequenceBasic) {
    const string sequence = "GGGCAGCGTAACAAT";
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GGC", 5},
        {"GCA", 6},
        {"CAG", 9},
        {"AGC", 17},
        {"GCG", 18},
        {"GGG", 1000}, // start
        {"CGT", 40},
        {"GTA", 103},
        {"TAA", 123},
        {"AAC", 9123},
        {"AAT", 902}, // end
        {"ACA", 1123},
        {"CAA", 999}
    };
    const uint32_t K = 3;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ(1000, jump_opt.value().start_k_mer_id);
        EXPECT_EQ(902, jump_opt.value().end_k_mer_id);
        EXPECT_EQ(11, jump_opt.value().nodes_in_between);
    } else {
        EXPECT_EQ("Jump constructed", "Jump was not constructed");
    }
}

TEST(JumpsTest, CreateJumpFromSequenceBasic2) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
    // const unordered_map<string, uint64_t> kmer_to_node_map = {
    //     {"GGGCAGCGTAACAATAATTTACG", 0},
    //     {"GGCAGCGTAACAATAATTTACGC", 1},
    //     {"GCAGCGTAACAATAATTTACGCA", 2},
    //     {"CAGCGTAACAATAATTTACGCAG", 3},
    //     {"AGCGTAACAATAATTTACGCAGA", 4},
    //     {"GCGTAACAATAATTTACGCAGAA", 5},
    //     {"CGTAACAATAATTTACGCAGAAA", 6},
    //     {"GTAACAATAATTTACGCAGAAAT", 7},
    //     {"TAACAATAATTTACGCAGAAATG", 8},
    //     {"AACAATAATTTACGCAGAAATGA", 9},
    //     {"ACAATAATTTACGCAGAAATGAT", 10},
    //     {"CAATAATTTACGCAGAAATGATG", 11},
    //     {"AATAATTTACGCAGAAATGATGT", 12},
    //     {"ATAATTTACGCAGAAATGATGTG", 13},
    //     {"TAATTTACGCAGAAATGATGTGT", 14},
    //     {"AATTTACGCAGAAATGATGTGTA", 15},
    //     {"ATTTACGCAGAAATGATGTGTAT", 16},
    //     {"TTTACGCAGAAATGATGTGTATT", 17},
    //     {"TTACGCAGAAATGATGTGTATTA", 18},
    //     {"TACGCAGAAATGATGTGTATTAT", 19},
    //     {"ACGCAGAAATGATGTGTATTATT", 20},
    //     {"CGCAGAAATGATGTGTATTATTA", 21},
    //     {"GCAGAAATGATGTGTATTATTAG", 22},
    //     {"CAGAAATGATGTGTATTATTAGG", 23},
    //     {"AGAAATGATGTGTATTATTAGGG", 24},
    //     {"GAAATGATGTGTATTATTAGGGC", 25},
    //     {"AAATGATGTGTATTATTAGGGCA", 26},
    //     {"AATGATGTGTATTATTAGGGCAT", 27},
    //     {"ATGATGTGTATTATTAGGGCATA", 28},
    //     {"TGATGTGTATTATTAGGGCATAA", 29},
    //     {"GATGTGTATTATTAGGGCATAAA", 30},
    //     {"ATGTGTATTATTAGGGCATAAAT", 31},
    //     {"TGTGTATTATTAGGGCATAAATA", 32},
    //     {"GTGTATTATTAGGGCATAAATAA", 33},
    //     {"TGTATTATTAGGGCATAAATAAA", 34},
    //     {"GTATTATTAGGGCATAAATAAAT", 35},
    //     {"TATTATTAGGGCATAAATAAATA", 36},
    //     {"ATTATTAGGGCATAAATAAATAA", 37},
    //     {"TTATTAGGGCATAAATAAATAAA", 38},
    //     {"TATTAGGGCATAAATAAATAAAT", 39},
    //     {"ATTAGGGCATAAATAAATAAATG", 40},
    //     {"TTAGGGCATAAATAAATAAATGT", 41},
    //     {"TAGGGCATAAATAAATAAATGTA", 42},
    //     {"AGGGCATAAATAAATAAATGTAT", 43},
    //     {"GGGCATAAATAAATAAATGTATG", 44},
    //     {"GGCATAAATAAATAAATGTATGT", 45},
    //     {"GCATAAATAAATAAATGTATGTA", 46},
    //     {"CATAAATAAATAAATGTATGTAT", 47},
    //     {"ATAAATAAATAAATGTATGTATG", 48},
    //     {"TAAATAAATAAATGTATGTATGT", 49},
    //     {"AAATAAATAAATGTATGTATGTG", 50},
    //     {"AATAAATAAATGTATGTATGTGT", 51},
    //     {"ATAAATAAATGTATGTATGTGTG", 52},
    //     {"TAAATAAATGTATGTATGTGTGT", 53},
    //     {"AAATAAATGTATGTATGTGTGTG", 54},
    //     {"AATAAATGTATGTATGTGTGTGT", 55},
    //     {"ATAAATGTATGTATGTGTGTGTA", 56},
    //     {"TAAATGTATGTATGTGTGTGTAT", 57},
    //     {"AAATGTATGTATGTGTGTGTATG", 58},
    //     {"AATGTATGTATGTGTGTGTATGT", 59},
    //     {"ATGTATGTATGTGTGTGTATGTA", 60},
    //     {"TGTATGTATGTGTGTGTATGTAT", 61},
    //     {"GTATGTATGTGTGTGTATGTATG", 62},
    //     {"TATGTATGTGTGTGTATGTATGA", 63},
    //     {"ATGTATGTGTGTGTATGTATGAG", 64},
    //     {"TGTATGTGTGTGTATGTATGAGT", 65},
    //     {"GTATGTGTGTGTATGTATGAGTC", 66},
    //     {"TATGTGTGTGTATGTATGAGTCG", 67},
    //     {"ATGTGTGTGTATGTATGAGTCGG", 68},
    //     {"TGTGTGTGTATGTATGAGTCGGT", 69},
    //     {"GTGTGTGTATGTATGAGTCGGTG", 70},
    //     {"TGTGTGTATGTATGAGTCGGTGA", 71},
    //     {"GTGTGTATGTATGAGTCGGTGAA", 72},
    //     {"TGTGTATGTATGAGTCGGTGAAT", 73},
    //     {"GTGTATGTATGAGTCGGTGAATG", 74},
    //     {"TGTATGTATGAGTCGGTGAATGA", 75},
    //     {"GTATGTATGAGTCGGTGAATGAC", 76},
    //     {"TATGTATGAGTCGGTGAATGACT", 77},
    //     {"ATGTATGAGTCGGTGAATGACTT", 78},
    //     {"TGTATGAGTCGGTGAATGACTTT", 79},
    //     {"GTATGAGTCGGTGAATGACTTTG", 80},
    //     {"TATGAGTCGGTGAATGACTTTGG", 81},
    //     {"ATGAGTCGGTGAATGACTTTGGT", 82},
    //     {"TGAGTCGGTGAATGACTTTGGTT", 83},
    //     {"GAGTCGGTGAATGACTTTGGTTT", 84},
    //     {"AGTCGGTGAATGACTTTGGTTTT", 85},
    //     {"GTCGGTGAATGACTTTGGTTTTA", 86},
    //     {"TCGGTGAATGACTTTGGTTTTAG", 87},
    //     {"CGGTGAATGACTTTGGTTTTAGT", 88},
    //     {"GGTGAATGACTTTGGTTTTAGTG", 89},
    //     {"GTGAATGACTTTGGTTTTAGTGA", 90},
    //     {"TGAATGACTTTGGTTTTAGTGAT", 91},
    //     {"GAATGACTTTGGTTTTAGTGATA", 92},
    //     {"AATGACTTTGGTTTTAGTGATAT", 93},
    //     {"ATGACTTTGGTTTTAGTGATATG", 94},
    //     {"TGACTTTGGTTTTAGTGATATGA", 95},
    //     {"GACTTTGGTTTTAGTGATATGAG", 96},
    //     {"ACTTTGGTTTTAGTGATATGAGT", 97},
    //     {"CTTTGGTTTTAGTGATATGAGTA", 98},
    //     {"TTTGGTTTTAGTGATATGAGTAG", 99},
    //     {"TTGGTTTTAGTGATATGAGTAGC", 100},
    //     {"TGGTTTTAGTGATATGAGTAGCC", 101},
    //     {"GGTTTTAGTGATATGAGTAGCCA", 102},
    //     {"GTTTTAGTGATATGAGTAGCCAC", 103}
    // };
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GGGCAGCGTAACAATAATTTACG", 0},
        {"GCGTAACAATAATTTACGCAGAA", 5},
        {"CGTAACAATAATTTACGCAGAAA", 6},
        {"AACAATAATTTACGCAGAAATGA", 9},
        {"TTTACGCAGAAATGATGTGTATT", 17},
        {"TTACGCAGAAATGATGTGTATTA", 18},
        {"ATTAGGGCATAAATAAATAAATG", 40},
        {"GTTTTAGTGATATGAGTAGCCAC", 103}
    };
    const uint32_t K = 23;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ(0, jump_opt.value().start_k_mer_id);
        EXPECT_EQ(103, jump_opt.value().end_k_mer_id);
        EXPECT_EQ(102, jump_opt.value().nodes_in_between);
    } else {
        EXPECT_EQ("Jump constructed", "Jump was not constructed");
    }
}

TEST(JumpsTest, CreateJumpFromSequenceBasicNotMatchingK) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GGGCAGCGTAACAA", 0},
        {"GCGTAACAATAATT", 5},
        {"CGTAACAATAATTT", 6},
        {"AACAATAATTTACG", 9},
        {"TTTACGCAGAAATG", 17},
        {"TTACGCAGAAATGA", 18},
        {"ATTAGGGCATAAAT", 40},
        {"GTTTTAGTGATATG", 103}
    };
    const uint32_t K = 30;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ("Jump not constructed", "Jump was constructed tho");
    } else {
        EXPECT_EQ(true, true);
    }
}

TEST(JumpsTest, CreateJumpFromSequenceBasicNotMatchingStart) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GCGTAACAATAATTTACGCAGAA", 5},
        {"CGTAACAATAATTTACGCAGAAA", 6},
        {"AACAATAATTTACGCAGAAATGA", 9},
        {"TTTACGCAGAAATGATGTGTATT", 17},
        {"TTACGCAGAAATGATGTGTATTA", 18},
        {"ATTAGGGCATAAATAAATAAATG", 40},
        {"GTTTTAGTGATATGAGTAGCCAC", 103}
    };
    const uint32_t K = 23;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ("Jump not constructed", "Jump was constructed tho");
    } else {
        EXPECT_EQ(true, true);
    }
}

TEST(JumpsTest, CreateJumpFromSequenceBasicNotMatchingEnd) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GGGCAGCGTAACAATAATTTACG", 0},
        {"GCGTAACAATAATTTACGCAGAA", 5},
        {"CGTAACAATAATTTACGCAGAAA", 6},
        {"AACAATAATTTACGCAGAAATGA", 9},
        {"TTTACGCAGAAATGATGTGTATT", 17},
        {"TTACGCAGAAATGATGTGTATTA", 18},
        {"ATTAGGGCATAAATAAATAAATG", 40}
    };
    const uint32_t K = 23;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ("Jump not constructed", "Jump was constructed tho");
    } else {
        EXPECT_EQ(true, true);
    }
}

TEST(JumpsTest, CreateJumpFromSequenceBasicNotMatchingBoth) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
    const unordered_map<string, uint64_t> kmer_to_node_map = {
        {"GCGTAACAATAATTTACGCAGAA", 5},
        {"CGTAACAATAATTTACGCAGAAA", 6},
        {"AACAATAATTTACGCAGAAATGA", 9},
        {"TTTACGCAGAAATGATGTGTATT", 17},
        {"TTACGCAGAAATGATGTGTATTA", 18},
        {"ATTAGGGCATAAATAAATAAATG", 40}
    };
    const uint32_t K = 23;

    optional<Jump> jump_opt = create_jump_from_sequence(sequence, kmer_to_node_map, K);

    if (jump_opt.has_value()) {
        EXPECT_EQ("Jump not constructed", "Jump was constructed tho");
    } else {
        EXPECT_EQ(true, true);
    }
}
