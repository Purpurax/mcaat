#include <gtest/gtest.h>
#include "jumps.h"

#include <string>
#include <unordered_map>

using namespace std;


// Normal test
TEST(JumpsTest, ReversePairEndsSequenceBasic) {
    const string input = "ACGT";
    const string expected = "ACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Normal test
TEST(JumpsTest, ReversePairEndsSequenceBasic2) {
    const string input = "ACGTT";
    const string expected = "AACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Normal test
TEST(JumpsTest, ReversePairEndsSequenceBasic3) {
    const string input = "AAGCT";
    const string expected = "AGCTT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Special test 'X' is in string
TEST(JumpsTest, ReversePairEndsSequenceNonStandardChars) {
    const string input = "AXGT";
    const string expected = "ACXT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Special test: Many non 'ACGT' chracters are in the string
TEST(JumpsTest, ReversePairEndsSequenceNonStandardChars2) {
    const string input = "XXXQUIOCPOPYM";
    const string expected = "MYPOPGOIUQXXX";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Normal case short kmers
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

// Normal case longer kmers
TEST(JumpsTest, CreateJumpFromSequenceBasic2) {
    const string sequence = "GGGCAGCGTAACAATAATTTACGCAGAAATGATGTGTATTATTAGGGCATAAATAAATAAATGTATGTATGTGTGTGTATGTATGAGTCGGTGAATGACTTTGGTTTTAGTGATATGAGTAGCCAC";
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

// Border case: Non-matching K
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

// Border case: start-kmer is not in the kmer_to_node_map
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

// Border case: end-kmer is not in the kmer_to_node_map
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

// Border case: start-kmer and end-kmer are not in the kmer_to_node_map
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
