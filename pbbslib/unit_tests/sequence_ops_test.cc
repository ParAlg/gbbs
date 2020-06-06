#include "pbbslib/sequence_ops.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "pbbslib/seq.h"

using ::testing::ElementsAre;

TEST(PackOut, PackIntoSequence) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(3, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers);
  EXPECT_THAT(even_numbers, ElementsAre(10, 12, 14));
}

TEST(PackOut, PackIntoSlice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(5, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers.slice(1, 4));
  EXPECT_THAT(even_numbers, ElementsAre(0, 10, 12, 14, 0));
}

TEST(PackOut, PackIntoReverseSlice) {
  pbbs::sequence<int32_t> numbers{5, [&](const size_t i) { return i + 10; }};
  pbbs::sequence<bool> even_flags{
    5, [&](const size_t i) { return i % 2 == 0; }};
  pbbs::sequence<int32_t> even_numbers(5, 0);
  pbbs::pack_out(numbers, even_flags, even_numbers.rslice(1, 4));
  EXPECT_THAT(even_numbers, ElementsAre(0, 14, 12, 10, 0));
}
