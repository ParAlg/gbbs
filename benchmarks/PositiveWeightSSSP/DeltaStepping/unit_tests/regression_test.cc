
#include "gbbs/bucket.h"
#include "gbbs/gbbs.h"
#include "gtest/gtest.h"
#include <condition_variable>
#include <mutex>
#include <thread>

namespace {
std::mutex gate;
std::condition_variable changed;
bool first_write = false, resume = false;
thread_local bool pause_after_write = false;
}
namespace gbbs {
template<class T> bool coordinated_write_min(T* address, T value) {
  bool result = gbbs::write_min(address, value);
  if (pause_after_write) {
    std::unique_lock<std::mutex> lock(gate);
    first_write = true;
    changed.notify_all();
    changed.wait(lock, [] { return resume; });
  }
  return result;
}
}
// Intercept only this functor's write to force the competing relaxation
// between the first decrease and its notification. Production code is unchanged.
#define write_min coordinated_write_min
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#undef write_min

TEST(DeltaStepping, NotificationCarriesRoundStartDistance) {
  gbbs::sequence<std::pair<gbbs::intE,bool>> distances(3);
  distances[0] = {0,false}; distances[1] = {0,false}; distances[2] = {100,false};
  gbbs::Visit_F<gbbs::intE,gbbs::intE> visit(distances);
  std::optional<gbbs::intE> first;
  std::thread worker([&] {
    pause_after_write = true;
    first = visit.updateAtomic(0,2,10);
  });
  {
    std::unique_lock<std::mutex> lock(gate);
    changed.wait(lock, [] { return first_write; });
  }
  auto second = visit.updateAtomic(1,2,5);
  {
    std::lock_guard<std::mutex> lock(gate);
    resume = true;
  }
  changed.notify_all();
  worker.join();
  ASSERT_TRUE(first.has_value());
  EXPECT_EQ(*first,100);
  EXPECT_FALSE(second.has_value());
  EXPECT_EQ(distances[2].first,5);
}
