#include "pbbslib/seq.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::ElementsAre;
using ::testing::IsEmpty;

TEST(Sequence, ConstructEmptySequence) {
  const pbbs::sequence<size_t> sequence{};
  EXPECT_THAT(sequence, IsEmpty());
}

TEST(Sequence, ConstructSequenceByFunction) {
  const pbbs::sequence<size_t> sequence{
    3,
    [](const size_t i) {
      return i;
    }};
  EXPECT_THAT(sequence, ElementsAre(0, 1, 2));
}

TEST(Sequence, Slice) {
  const pbbs::sequence<size_t> sequence{
    5,
    [](const size_t i) {
      return i;
    }};

  {
    const pbbs::range<size_t*> slice{sequence.slice(2, 2)};
    EXPECT_THAT(slice, IsEmpty());
  }
  {
    const pbbs::range<size_t*> slice{sequence.slice(1, 4)};
    EXPECT_THAT(slice, ElementsAre(1, 2, 3));
  }
}
