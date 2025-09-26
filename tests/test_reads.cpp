#include <gtest/gtest.h>
#include "reads.h"

#include <string>
#include <unordered_map>

using namespace std;


// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic) {
    const string input = "ACGT";
    const string expected = "ACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic2) {
    const string input = "ACGTT";
    const string expected = "AACGT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Normal test
TEST(ReadsTest, ReversePairEndsSequenceBasic3) {
    const string input = "AAGCT";
    const string expected = "AGCTT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Special test 'X' is in string
TEST(ReadsTest, ReversePairEndsSequenceNonStandardChars) {
    const string input = "AXGT";
    const string expected = "ACXT";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}

// Special test: Many non 'ACGT' chracters are in the string
TEST(ReadsTest, ReversePairEndsSequenceNonStandardChars2) {
    const string input = "XXXQUIOCPOPYM";
    const string expected = "MYPOPGOIUQXXX";
    EXPECT_EQ(reverse_pair_ends_sequence(input), expected);
}
