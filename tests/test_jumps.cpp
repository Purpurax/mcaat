#include <gtest/gtest.h>
#include "jumps.h"

TEST(JumpsTest, ReversePairEndsSequenceBasic) {
    // "ACGT" reversed is "TGCA", then complemented: T->A, G->C, C->G, A->T => "ACGT"
    std::string input = "ACGT";
    std::string expected = "ACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceMixed) {
    // "AAGCT" reversed is "TCGAA", then complemented: T->A, C->G, G->C, A->T, A->T => "AGCTT"
    std::string input = "AAGCT";
    std::string expected = "AGCTT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

TEST(JumpsTest, ReversePairEndsSequenceNonStandardChars) {
    // Non-standard chars should remain unchanged
    std::string input = "AXGT";
    std::string expected = "ACXT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}
