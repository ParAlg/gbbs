#include "gbbs/helpers/parallel_for_with_status.h"

#include <cstddef>
#include <vector>

#include "benchmark/benchmark.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "status_macros.h"
#include "absl/base/thread_annotations.h"
#include "absl/log/absl_log.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/synchronization/mutex.h"
#include "parlay/parallel.h"

namespace gbbs {
namespace {

using ::testing::AnyOf;
using ::testing::UnorderedElementsAre;

class ThreadSafeCounter {
 public:
  void Increment() {
    absl::MutexLock lock(&mutex_);
    ++value_;
  }

  int Value() {
    absl::MutexLock lock(&mutex_);
    return value_;
  }

 private:
  absl::Mutex mutex_;
  int value_ ABSL_GUARDED_BY(mutex_) = 0;
};

class ParallelForWithStatusTest : public ::testing::Test {
// private:
//  graph_mining::in_memory::ParallelSchedulerReference scheduler_;
};

TEST_F(ParallelForWithStatusTest, EmptyRange) {
  ThreadSafeCounter counter;
  EXPECT_OK(gbbs::parallel_for_with_status(5, 5, [&](size_t i) {
    counter.Increment();
    ABSL_LOG(FATAL) << "This function should not have been called, because the "
                       "range is empty";
    return absl::OkStatus();
  }));
  EXPECT_EQ(counter.Value(), 0);
}

TEST_F(ParallelForWithStatusTest, SingleCall_Ok) {
  ThreadSafeCounter counter;
  EXPECT_OK(gbbs::parallel_for_with_status(3, 4, [&](size_t i) {
    counter.Increment();
    return absl::OkStatus();
  }));
  EXPECT_EQ(counter.Value(), 1);
}

TEST_F(ParallelForWithStatusTest, SingleCall_Error) {
  ThreadSafeCounter counter;
  EXPECT_THAT(
      gbbs::parallel_for_with_status(3, 4,
                                     [&](size_t i) {
                                       counter.Increment();
                                       return absl::NotFoundError("foo");
                                     }),
      StatusIs(absl::StatusCode::kNotFound, "foo"));
  EXPECT_EQ(counter.Value(), 1);
}

TEST_F(ParallelForWithStatusTest, MultipleCalls_Ok) {
  absl::Mutex indices_seen_mutex;
  std::vector<size_t> indices_seen;  // Guarded by `indices_seen_mutex`.
  EXPECT_OK(gbbs::parallel_for_with_status(5, 15, [&](size_t i) {
    {
      absl::MutexLock lock(&indices_seen_mutex);
      indices_seen.push_back(i);
    }
    return absl::OkStatus();
  }));
  EXPECT_THAT(indices_seen,
              UnorderedElementsAre(5, 6, 7, 8, 9, 10, 11, 12, 13, 14));
}

TEST_F(ParallelForWithStatusTest, MultipleCalls_AllErrors) {
  ThreadSafeCounter counter;
  EXPECT_THAT(gbbs::parallel_for_with_status(1, 5,
                                             [&](size_t i) {
                                               counter.Increment();
                                               return absl::InternalError(
                                                   absl::StrCat("foo_", i));
                                             }),
              StatusIs(absl::StatusCode::kInternal,
                       AnyOf("foo_1", "foo_2", "foo_3", "foo_4")));
  EXPECT_EQ(counter.Value(), 4);
}

TEST_F(ParallelForWithStatusTest, MultipleCalls_SingleError) {
  ThreadSafeCounter counter;
  EXPECT_THAT(gbbs::parallel_for_with_status(
                  2, 8,
                  [&](size_t i) {
                    counter.Increment();
                    if (i == 4) {
                      return absl::FailedPreconditionError("foo");
                    } else {
                      return absl::OkStatus();
                    }
                  }),
              StatusIs(absl::StatusCode::kFailedPrecondition, "foo"));
  EXPECT_EQ(counter.Value(), 6);
}

// Microbenchmarks comparing `parlay::parallel_for` and
// `gbbs::parallel_for_with_status`, in the case where `f` does a tiny amount
// of work (in which case the overhead of handling statuses might be
// non-negligible), and `f` always returns an OK status (which is the most
// common case).
//
// To run these microbenchmarks, use:
/*
   blaze --blazerc=/dev/null build --config=benchmark \
     //gbbs/helpers:parallel_for_with_status_test

   blaze-bin/gbbs/helpers/parallel_for_with_status_test \
   --benchmark_filter=all --benchmark_min_time=60s
*/

static void BM_ParallelForPlain(benchmark::State& state) {
//  graph_mining::in_memory::ParallelSchedulerReference scheduler;
  int n = state.range(0);
  std::vector<int> v(n);
  for (auto s : state) {
    benchmark::DoNotOptimize(v);
    parlay::parallel_for(0, n, [&](size_t i) { v[i] = i; });
  }
}

BENCHMARK(BM_ParallelForPlain)->Arg(1'000'000)->Arg(10'000'000);

static void BM_ParallelForWithStatus(benchmark::State& state) {
//  graph_mining::in_memory::ParallelSchedulerReference scheduler;
  int n = state.range(0);
  std::vector<int> v(n);
  for (auto s : state) {
    benchmark::DoNotOptimize(v);
    EXPECT_OK(gbbs::parallel_for_with_status(0, n, [&](size_t i) {
      v[i] = i;
      return absl::OkStatus();
    }));
  }
}

BENCHMARK(BM_ParallelForWithStatus)->Arg(1'000'000)->Arg(10'000'000);

}  // namespace
}  // namespace gbbs
