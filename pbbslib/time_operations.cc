#include "time_operations.h"

typedef unsigned __int128 long_int;

#if defined(CILK)
#include "reducer.h"
double t_histogram_reducer(size_t n, bool check) {
  return 1.0;
//  pbbs::random r(0);
//  constexpr int count = 1024;
//  histogram_reducer<int, count> red;
//  using aa = std::array<size_t, 8>;
//  pbbs::sequence<aa> In(n, [&](size_t i) {
//    aa x;
//    x[0] = r.ith_rand(i) % count;
//    return x;
//  });
//  auto f = [&](size_t i) { red->add_value(In[i][0]); };
//  time(t, parallel_for(0, n, f););
//  // std::cout << red.get_value()[0] << std::endl;
//  return t;
}
#else
double t_histogram_reducer(size_t n, bool check) { return 1.0; }
#endif

double t_map_reduce_128(size_t n, bool check) {
  int stride = 16;
  pbbs::sequence<size_t> S(n * stride, (size_t)1);
  auto get = [&](size_t i) {
    // gives marginal improvement (5% or so on aware)
    __builtin_prefetch(&S[(i + 4) * stride], 0, 3);
    return S[i * stride];
  };
  auto T = pbbs::delayed_seq<size_t>(n, get);
  time(t, pbbs::reduce(T, pbbs::addm<size_t>()););
  return t;
}

double t_integer_sort_128(size_t n, bool check) {
  pbbs::random r(0);
  size_t bits = 128;
  pbbs::sequence<long_int> S(n, [&](size_t i) -> long_int {
    return r.ith_rand(2 * i) + (((long_int)r.ith_rand(2 * i + 1)) << 64);
  });
  auto identity = [](long_int a) { return a; };
  time(t, pbbs::integer_sort_inplace(S.slice(), identity, bits););
  return t;
}
