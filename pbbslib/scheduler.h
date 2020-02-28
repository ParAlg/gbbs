#pragma once

#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <thread>

// EXAMPLE USE 1:
//
// fork_join_scheduler fj;
//
// long fib(long i) {
//   if (i <= 1) return 1;
//   long l,r;
//   fj.pardo([&] () { l = fib(i-1);},
//            [&] () { r = fib(i-2);});
//   return l + r;
// }
//
// fib(40);
//
// EXAMPLE USE 2:
//
// void init(long* x, size_t n) {
//   parfor(0, n, [&] (int i) {a[i] = i;});
// }
//
// size_t n = 1000000000;
//
// long* a = new long[n];
//
// init(a, n);

// Deque from Arora, Blumofe, and Plaxton (SPAA, 1998).
template <typename Job>
struct Deque {
  using qidx = unsigned int;
  using tag_t = unsigned int;

  // use unit for atomic access
  union age_t {
    struct {
      tag_t tag;
      qidx top;
    } pair;
    size_t unit;
  };

  // align to avoid false sharing
  struct alignas(64) padded_job {
    Job* job;
  };

  static int const q_size = 200;
  age_t age;
  qidx bot;
  padded_job deq[q_size];

  inline bool cas(size_t* ptr, size_t oldv, size_t newv) {
    return __sync_bool_compare_and_swap(ptr, oldv, newv);
  }

  inline void fence() { std::atomic_thread_fence(std::memory_order_seq_cst); }

  Deque() : bot(0) {
    age.pair.tag = 0;
    age.pair.top = 0;
  }

  void push_bottom(Job* job) {
    qidx local_bot;
    local_bot = bot;           // atomic load
    deq[local_bot].job = job;  // shared store
    local_bot += 1;
    if (local_bot == q_size) {
      std::cout << "internal error: scheduler queue overflow" << std::endl;
      abort();
    }
    bot = local_bot;  // shared store
    fence();
  }

  Job* pop_top() {
    age_t old_age, new_age;
    qidx local_bot;
    Job *job, *result;
    old_age.unit = age.unit;  // atomic load

    local_bot = bot;  // atomic load
    if (local_bot <= old_age.pair.top)
      result = NULL;
    else {
      job = deq[old_age.pair.top].job;  // atomic load
      new_age.unit = old_age.unit;
      new_age.pair.top = new_age.pair.top + 1;
      if (cas(&(age.unit), old_age.unit, new_age.unit))  // cas
        result = job;
      else
        result = NULL;
    }
    return result;
  }

  Job* pop_bottom() {
    age_t old_age, new_age;
    qidx local_bot;
    Job *job, *result;
    local_bot = bot;  // atomic load
    if (local_bot == 0)
      result = NULL;
    else {
      local_bot = local_bot - 1;
      bot = local_bot;  // shared store
      fence();
      job = deq[local_bot].job;  // atomic load
      old_age.unit = age.unit;   // atomic load
      if (local_bot > old_age.pair.top)
        result = job;
      else {
        bot = 0;  // shared store
        new_age.pair.top = 0;
        new_age.pair.tag = old_age.pair.tag + 1;
        if ((local_bot == old_age.pair.top) &&
            cas(&(age.unit), old_age.unit, new_age.unit))
          result = job;
        else {
          age.unit = new_age.unit;  // shared store
          result = NULL;
        }
        fence();
      }
    }
    return result;
  }
};

// thread_local int thread_id;

template <typename Job>
struct scheduler {
 public:
  int num_threads;

  static thread_local int thread_id;

  scheduler() {
    init_num_workers();
    num_deques = 2 * num_threads;
    deques = new Deque<Job>[num_deques];
    attempts = new attempt[num_deques];
    finished_flag = 0;

    // Spawn num_workers many threads on startup
    spawned_threads = new std::thread[num_threads - 1];
    std::function<bool()> finished = [&]() { return finished_flag == 1; };
    thread_id = 0;  // thread-local write
    for (int i = 1; i < num_threads; i++) {
      spawned_threads[i - 1] = std::thread([&, i, finished]() {
        thread_id = i;  // thread-local write
        start(finished);
      });
    }
  }

  ~scheduler() {
    finished_flag = 1;
    for (int i = 1; i < num_threads; i++) {
      spawned_threads[i - 1].join();
    }
    delete[] spawned_threads;
    delete[] deques;
    delete[] attempts;
  }

  // Push onto local stack.
  void spawn(Job* job) {
    int id = worker_id();
    deques[id].push_bottom(job);
  }

  // Wait for condition: finished().
  template <typename F>
  void wait(F finished, bool conservative = false) {
    // Conservative avoids deadlock if scheduler is used in conjunction
    // with user locks enclosing a wait.
    if (conservative)
      while (!finished()) std::this_thread::yield();
    // If not conservative, schedule within the wait.
    // Can deadlock if a stolen job uses same lock as encloses the wait.
    else
      start(finished);
  }

  // All scheduler threads quit after this is called.
  void finish() { finished_flag = 1; }

  // Pop from local stack.
  Job* try_pop() {
    int id = worker_id();
    return deques[id].pop_bottom();
  }

  void init_num_workers() {
    if (const char* env_p = std::getenv("NUM_THREADS")) {
      num_threads = std::stoi(env_p);
    } else {
      num_threads = std::thread::hardware_concurrency();
    }
  }

  int num_workers() { return num_threads; }
  int worker_id() { return thread_id; }
  void set_num_workers(int n) {
    std::cout << "Unsupported" << std::endl;
    exit(-1);
  }

 private:
  // Align to avoid false sharing.
  struct alignas(128) attempt {
    size_t val;
  };

  int num_deques;
  Deque<Job>* deques;
  attempt* attempts;
  std::thread* spawned_threads;
  int finished_flag;

  // Start an individual scheduler task.  Runs until finished().
  template <typename F>
  void start(F finished) {
    while (1) {
      Job* job = get_job(finished);
      if (!job) return;
      (*job)();
    }
  }

  Job* try_steal(size_t id) {
    // use hashing to get "random" target
    size_t target = (hash(id) + hash(attempts[id].val)) % num_deques;
    attempts[id].val++;
    return deques[target].pop_top();
  }

  // Find a job, first trying local stack, then random steals.
  template <typename F>
  Job* get_job(F finished) {
    if (finished()) return NULL;
    Job* job = try_pop();
    if (job) return job;
    size_t id = worker_id();
    while (1) {
      // By coupon collector's problem, this should touch all.
      for (int i = 0; i <= num_deques * 100; i++) {
        if (finished()) return NULL;
        job = try_steal(id);
        if (job) return job;
      }
      // If haven't found anything, take a breather.
      std::this_thread::sleep_for(std::chrono::nanoseconds(num_deques * 100));
    }
  }

  uint64_t hash(uint64_t x) {
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
  }
};

template <typename T>
thread_local int scheduler<T>::thread_id = 0;

struct fork_join_scheduler {
 public:
  // Jobs are thunks -- i.e., functions that take no arguments
  // and return nothing.   Could be a lambda, e.g. [] () {}.
  using Job = std::function<void()>;

  scheduler<Job>* sched;

  fork_join_scheduler();

  ~fork_join_scheduler();

  // Must be called using std::atexit(..) to free resources
  void destroy();

  int num_workers();
  int worker_id();
  void set_num_workers(int n);

  // Fork two thunks and wait until they both finish.
  template <typename L, typename R>
  void pardo(L left, R right, bool conservative = false) {
    bool right_done = false;
    Job right_job = [&]() {
      right();
      right_done = true;
    };
    sched->spawn(&right_job);
    left();
    if (sched->try_pop() != NULL)
      right();
    else {
      auto finished = [&]() { return right_done; };
      sched->wait(finished, conservative);
    }
  }

  template <typename F>
  int get_granularity(size_t start, size_t end, F f) {
    size_t done = 0;
    size_t size = 1;
    int ticks;
    do {
      size = std::min(size, end - (start + done));
      auto tstart = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < size; i++) f(start + done + i);
      auto tstop = std::chrono::high_resolution_clock::now();
      ticks = (tstop - tstart).count();
      done += size;
      size *= 2;
    } while (ticks < 1000 && done < (end - start));
    return done;
  }

  template <typename F>
  void parfor(size_t start, size_t end, F f, size_t granularity = 0,
              bool conservative = false) {
    if (granularity == 0) {
      size_t done = get_granularity(start, end, f);
      granularity = std::max(done, (end - start) / (128 * sched->num_threads));
      parfor_(start + done, end, f, granularity, conservative);
    } else
      parfor_(start, end, f, granularity, conservative);
  }

 private:
  template <typename F>
  void parfor_(size_t start, size_t end, F f, size_t granularity,
               bool conservative) {
    if ((end - start) <= granularity)
      for (size_t i = start; i < end; i++) f(i);
    else {
      size_t n = end - start;
      // Not in middle to avoid clashes on set-associative caches
      // on powers of 2.
      size_t mid = (start + (9 * (n + 1)) / 16);
      pardo([&]() { parfor_(start, mid, f, granularity, conservative); },
            [&]() { parfor_(mid, end, f, granularity, conservative); },
            conservative);
    }
  }
};

// Global fork-join scheduler.
extern fork_join_scheduler& global_scheduler;

namespace internal {

// Nifty Counter idiom for maintaining `global_scheduler` so that the scheduler
// is initialized before first use and destroyed after last use.
//
// This is to avoid the "static initialization/de-initialization order fiasco";
// the scheduler is liable to be used by the constructors and destructors of
// global or static variable, so we need to be very careful about the
// construction and destruction of the scheduler.
struct SchedulerInitializer {
  SchedulerInitializer();
  ~SchedulerInitializer();
};
static SchedulerInitializer global_scheduler_initializer;

}  // namespace internal
