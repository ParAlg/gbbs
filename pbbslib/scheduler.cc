#include "scheduler.h"

fork_join_scheduler::fork_join_scheduler() {
  sched = new scheduler<Job>;
}

fork_join_scheduler::~fork_join_scheduler() {
  if (sched) {
    delete sched;
    sched = nullptr;
  }
}

// Must be called using std::atexit(..) to free resources
void fork_join_scheduler::destroy() {
  if (sched) {
    delete sched;
    sched = nullptr;
  }
}

int fork_join_scheduler::num_workers() { return sched->num_workers(); }
int fork_join_scheduler::worker_id() { return sched->worker_id(); }
void fork_join_scheduler::set_num_workers(int n) { sched->set_num_workers(n); }
