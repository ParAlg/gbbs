#include "pbbslib/sequence_ops.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "pbbslib/monoid.h"
#include "pbbslib/seq.h"

using ::testing::ElementsAre;

TEST(ScanInplace, Sequence) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  int32_t sum{pbbs::scan_inplace(numbers, pbbs::addm<int32_t>{})};
  EXPECT_EQ(sum, 60);
  EXPECT_THAT(numbers, ElementsAre(0, 10, 21, 33, 46));
}

TEST(ScanInplace, Slice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  int32_t sum{pbbs::scan_inplace(numbers.slice(1, 4), pbbs::addm<int32_t>{})};
  EXPECT_EQ(sum, 36);
  EXPECT_THAT(numbers, ElementsAre(10, 0, 11, 23, 14));
}

TEST(PackOut, OutputIsSequence) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(3, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers);
  EXPECT_THAT(even_numbers, ElementsAre(10, 12, 14));
}

TEST(PackOut, OutputIsSlice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(5, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers.slice(1, 4));
  EXPECT_THAT(even_numbers, ElementsAre(0, 10, 12, 14, 0));
}

TEST(PackOut, OutputIsReverseSlice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(5, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers.rslice(1, 4));
  EXPECT_THAT(even_numbers, ElementsAre(0, 14, 12, 10, 0));
}

TEST(FilterOut, OutputIsSequence) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  constexpr auto is_even{[](const int32_t x) { return x % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(3, 0);
  pbbs::filter_out(numbers, even_numbers, is_even);
  EXPECT_THAT(even_numbers, ElementsAre(10, 12, 14));
}

TEST(FilterOut, OutputIsSlice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  constexpr auto is_even{[](const int32_t x) { return x % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(5, 0);
  pbbs::filter_out(numbers, even_numbers.slice(1, 4), is_even);
  EXPECT_THAT(even_numbers, ElementsAre(0, 10, 12, 14, 0));
}
