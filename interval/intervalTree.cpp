#include <stdio.h>
#include <vector>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include <parlay/random.h>

#ifdef USE_PAM
#include <pam/pam.h>
#include <pam/get_time.h>
#include <pam/parse_command_line.h>
#else
#include <cpam/pam.h>
#include <cpam/get_time.h>
#include <cpam/parse_command_line.h>
using namespace cpam;
#endif

using namespace std;

using point = int;
using par = std::tuple<point,point>;

struct interval_map {

  using interval = std::tuple<point, point>;

  struct entry {
    using key_t = point;
    using val_t = point;
    using aug_t = interval;
    static inline bool comp(key_t a, key_t b) { return a < b;}
    static aug_t get_empty() { return interval(0,0);}
    static aug_t from_entry(key_t k, val_t v) { return interval(k,v);}
    static aug_t combine(aug_t a, aug_t b) {
      return (std::get<1>(a) > std::get<1>(b)) ? a : b;}
  };

  using amap = aug_map<entry>;
  amap m;

  static void reserve(size_t n) { amap::reserve(n); }
  static void finish() { amap::finish(); }

  interval_map(parlay::sequence<interval> const &A) {
    m = amap(A); }

  bool stab(point p) {
    return (std::get<1>(m.aug_left(p)) > p);}

//  bool check_structure() {
//    return m.check_structure();
//  }

  void insert(interval i) {m.insert(i);}
  void insertBatch(const parlay::sequence<interval>& updates) {
    m = amap::multi_insert(std::move(m), updates);
  }

  vector<interval> report_all(point p) {
    vector<interval> vec;
    amap a = m;
    interval I = a.aug_left(p);
    while (std::get<1>(I) > p) {
      vec.push_back(I);
      a = amap::remove(move(a),std::get<0>(I));
      I = a.aug_left(p); }
    return vec; }

  void remove_small(point l) {
    auto f = [&] (interval I) {
      return (std::get<1>(I) - std::get<0>(I) >= l);};
    m = amap::filter(move(m),f); }
};

void test_all(size_t n, size_t q_num, size_t rounds) {
  parlay::sequence<par> v(n);
  parlay::sequence<par> vv(n);
  size_t max_size = (((size_t) 1) << 31)-1;

  parlay::random r(0);
  parlay::parallel_for(0, n, [&] (size_t i) {
    point start = r.ith_rand(2*i)%(max_size/2);
    point end = start + r.ith_rand(2*i+1)%(max_size-start);
    v[i] = make_tuple(start, end);
  });

  parlay::sequence<bool> result(q_num);
  parlay::sequence<int> queries(q_num);
  parlay::parallel_for(0, q_num, [&] (size_t i) {
    queries[i] = r.ith_rand(6*i)%max_size;
  });

  parlay::sequence<double> build_tm(rounds);
  parlay::sequence<double> query_tm(rounds);
  const size_t threads = parlay::num_workers();
  for (size_t i=0; i < rounds+1; i++) {
    {
      parlay::parallel_for (0, n, [&] (size_t i) {vv[i] = v[i];});
      interval_map::reserve(n);
      timer t;
      t.start();
      interval_map itree(vv);
      double tm = t.stop();
      std::cout << "Finished build in " << tm << std::endl;
      if (i>0) build_tm[i-1] = tm;

      timer tq;
      tq.start();
      parlay::parallel_for(0, q_num, [&] (size_t i) {
          result[i] = itree.stab(queries[i]);});
      double tm2 = tq.stop();
      std::cout << "Finished queries in " << tm2 << std::endl;
      auto res_map = parlay::delayed_seq<size_t>(q_num, [&] (size_t i) {
        return (size_t)result[i];
      });
      std::cout << "Sum result = " << parlay::reduce(res_map) << std::endl;
      if (i>0) query_tm[i-1] = tm2;
    }
    interval_map::finish();
  }

  auto less = [] (double a, double b) {return a < b;};
  parlay::sort(build_tm,less);
  parlay::sort(query_tm,less);

  cout << "interval build"
       << ", threads = " << threads
       << ", rounds = " << rounds
       << ", n = " << n
       << ", q = " << q_num
       << ", time = " << build_tm[rounds/2] << endl;

  cout << "interval query"
       << ", threads = " << threads
       << ", rounds = " << rounds
       << ", n = " << n
       << ", q = " << q_num
       << ", time = " << query_tm[rounds/2] << endl;

}

void test_updates(size_t n, size_t m, size_t rounds) {
  parlay::sequence<par> v(n);
  parlay::sequence<par> vv(n);
  size_t max_size = (((size_t) 1) << 31)-1;

  parlay::random r(0);
  parlay::parallel_for(0, n, [&] (size_t i) {
    point start = r.ith_rand(2*i)%(max_size/2);
    point end = start + r.ith_rand(2*i+1)%(max_size-start);
    v[i] = make_tuple(start, end);
  });

  parlay::sequence<par> updates(m);
  r = r.next();
  parlay::parallel_for(0, m, [&] (size_t i) {
    point start = r.ith_rand(2*i)%(max_size/2);
    point end = start + r.ith_rand(2*i+1)%(max_size-start);
    updates[i] = make_tuple(start, end);
  });

  parlay::sequence<double> build_tm(rounds);
  parlay::sequence<double> update_tm(rounds);
  const size_t threads = parlay::num_workers();
  for (size_t i=0; i < rounds+1; i++) {
    {
      parlay::parallel_for (0, n, [&] (size_t i) {vv[i] = v[i];});
      interval_map::reserve(n);
      timer t;
      t.start();
      interval_map itree(vv);
      double tm = t.stop();
      std::cout << "Finished build in " << tm << std::endl;
      if (i>0) build_tm[i-1] = tm;

      timer tup;
      tup.start();
      itree.insertBatch(updates);
      double tm2 = tup.stop();
      std::cout << "Finished updates in " << tm2 << std::endl;
      if (i>0) update_tm[i-1] = tm2;
    }
    interval_map::finish();
  }

  auto less = [] (double a, double b) {return a < b;};
  parlay::sort(build_tm,less);
  parlay::sort(update_tm,less);

  cout << "interval build"
       << ", threads = " << threads
       << ", rounds = " << rounds
       << ", n = " << n
       << ", m = " << m
       << ", time = " << build_tm[rounds/2] << endl;

  cout << "interval update"
       << ", threads = " << threads
       << ", rounds = " << rounds
       << ", n = " << n
       << ", m = " << m
       << ", time = " << update_tm[rounds/2] << endl;
}


int main(int argc, char** argv) {
  commandLine P(argc, argv,
		"./intervalTree [-n size] [-q num_queries] [-r rounds]");
  size_t n = P.getOptionLongValue("-n", 100000000);
  size_t q_num = P.getOptionLongValue("-q", 10000000);
  size_t rounds = P.getOptionIntValue("-r", 5);
  if (P.getOptionValue("-run_updates")) {
    size_t m = P.getOptionIntValue("-m", 10000000);
    test_updates(n, m, rounds);
  } else {
    test_all(n, q_num, rounds);
  }

  return 0;
}
