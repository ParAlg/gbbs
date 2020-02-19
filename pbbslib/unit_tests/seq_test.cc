#include "pbbslib/seq.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::ElementsAre;
using ::testing::IsEmpty;

TEST(Sequence, ConstructEmptySequence) {
  const pbbs::sequence<size_t> sequence{};
  EXPECT_THAT(sequence, IsEmpty());
}
