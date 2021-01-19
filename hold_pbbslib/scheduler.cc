#include "scheduler.h"

namespace pbbs {
namespace {

int global_scheduler_counter;  // Zero-initialized at load time.
typename std::aligned_storage<sizeof(fork_join_scheduler),
                              alignof(fork_join_scheduler)>::type
    global_scheduler_storage;

}  // namespace

namespace internal {

SchedulerInitializer::SchedulerInitializer() {
  if (global_scheduler_counter == 0) {
    new (&global_scheduler) fork_join_scheduler{};  // placement new
  }
  global_scheduler_counter++;
}

SchedulerInitializer::~SchedulerInitializer() {
  global_scheduler_counter--;
  if (global_scheduler_counter == 0) {
    (&global_scheduler)->~fork_join_scheduler();
  }
}

}  // namespace internal

fork_join_scheduler& global_scheduler{
    reinterpret_cast<fork_join_scheduler&>(global_scheduler_storage)};

fork_join_scheduler::fork_join_scheduler() { sched = new scheduler<Job>; }

fork_join_scheduler::~fork_join_scheduler() {
  if (sched) {
    delete sched;
    sched = nullptr;
  }
}

int fork_join_scheduler::num_workers() { return sched->num_workers(); }
int fork_join_scheduler::worker_id() { return sched->worker_id(); }
void fork_join_scheduler::set_num_workers(int n) { sched->set_num_workers(n); }

}  // namespace pbbs
