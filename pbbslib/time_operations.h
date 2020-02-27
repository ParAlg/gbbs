#pragma once

#include "bag.h"
#include "collect_reduce.h"
#include "counting_sort.h"
#include "get_time.h"
#include "hash_table.h"
#include "histogram.h"
#include "integer_sort.h"
#include "merge.h"
#include "merge_sort.h"
#include "monoid.h"
#include "random.h"
#include "random_shuffle.h"
#include "range_min.h"
#include "sample_sort.h"
#include "sequence_ops.h"
#include "sparse_mat_vec_mult.h"
#include "stlalgs.h"
#include "utilities.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <iostream>

static timer bt;
using uchar = unsigned char;

#define time(_var, _body) \
  bt.start();             \
  _body;                  \
  double _var = bt.stop();

template <typename T>
double t_tabulate(size_t n, bool check) {
  auto f = [](size_t i) { return i; };
  time(t, pbbs::sequence<T>(n, f););
  return t;
}

template <typename T>
double t_map(size_t n, bool check) {
  pbbs::sequence<T> In(n, (T)1);
  auto f = [&](size_t i) { return In[i]; };
  time(t, pbbs::sequence<T>(n, f));
  return t;
}

template <typename T>
double t_reduce_add(size_t n, bool check) {
  pbbs::sequence<T> S(n, (T)1);
  // time(t, sum(S.begin(), S.size()););
  // time(t, pbbs::reduce_add(S););
  time(t, pbbs::reduce(S, pbbs::addm<T>()));
  return t;
}

double t_map_reduce_128(size_t n, bool check);

template <typename T>
double t_scan_add(size_t n, bool check) {
  pbbs::sequence<T> In(n, (T)1);
  pbbs::sequence<T> Out;
  T sum;
  time(t, std::tie(Out, sum) = pbbs::scan(In, pbbs::addm<T>()););
  return t;
}

template <typename T>
double t_scan_add_seq(size_t n, bool check) {
  pbbs::sequence<T> In(n, (T)1);
  pbbs::sequence<T> Out(n, (T)0);
  time(t, T x = 0; for (size_t i = 0; i < n; i++) {
    Out[i] = x;
    x += In[i];
  });
  return t;
}

template <typename T>
double t_pack(size_t n, bool check) {
  pbbs::sequence<bool> flags(n, [](size_t i) -> bool { return i % 2; });
  pbbs::sequence<T> In(n, [](size_t i) -> T { return i; });
  time(t, pbbs::pack(In, flags););
  return t;
}

template <typename T>
double t_split3_old(size_t n, bool check) {
  pbbs::sequence<uchar> flags(n, [](size_t i) -> uchar { return i % 3; });
  pbbs::sequence<T> In(n, [](size_t i) -> T { return i; });
  pbbs::sequence<T> Out(n, (T)0);
  time(t, pbbs::split_three(In, Out, flags););
  return t;
}

template <typename T>
double t_split3(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> In(n, [&](size_t i) { return r.ith_rand(i); });
  pbbs::sequence<T> Out(n, (T)0);
  time(t, pbbs::p_split3(In, Out.slice(), std::less<T>()););
  return t;
}

template <typename T>
double t_partition(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> In(n,
                       [&](size_t i) -> size_t { return r.ith_rand(i) % n; });
  pbbs::sequence<T> Out;
  auto f = pbbs::delayed_seq<bool>(
      n, [&](size_t i) { return In[i] < n / ((size_t)2); });
  size_t m;
  time(t, std::tie(Out, m) = pbbs::split_two(In, f););
  return t;
}

double t_histogram_reducer(size_t n, bool check);

template <typename T>
double t_gather(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return i; });
  pbbs::sequence<T> idx(n, [&](size_t i) { return r.ith_rand(i) % n; });
  auto f = [&](size_t i) {
    // prefetching helps significantly
    __builtin_prefetch(&in[idx[i + 4]], 0, 1);
    return in[idx[i]];
  };
  // note problem with prefetching since will go over end for last 4 iterations
  time(t, pbbs::sequence<T>(n - 4, f););
  return t;
}

template <typename T>
double t_scatter(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> out(n, (T)0);
  pbbs::sequence<T> idx(n, [&](size_t i) { return r.ith_rand(i) % n; });
  auto f = [&](size_t i) {
    // prefetching makes little if any difference
    //__builtin_prefetch (&out[idx[i+4]], 1, 1);
    out[idx[i]] = i;
  };
  time(t, parallel_for(0, n - 4, f););
  return t;
}

template <typename T>
double t_write_add(size_t n, bool check) {
  pbbs::random r(0);
  // pbbs::sequence<T> out(n, (T) 0);
  pbbs::sequence<std::atomic<T>> out(n);
  parallel_for(0, n, [&](size_t i) { std::atomic_init(&out[i], (T)0); });
  pbbs::sequence<T> idx(n, [&](size_t i) { return r.ith_rand(i) % n; });
  auto f = [&](size_t i) {
    // putting write prefetch in slows it down
    //__builtin_prefetch (&out[idx[i+4]], 0, 1);
    //__sync_fetch_and_add(&out[idx[i]],1);};
    pbbs::write_add(&out[idx[i]], 1);
  };
  time(t, parallel_for(0, n - 4, f););
  return t;
}

template <typename T>
double t_write_min(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<std::atomic<T>> out(n);
  parallel_for(0, n, [&](size_t i) { std::atomic_init(&out[i], (T)n); });
  // pbbs::sequence<T> out(n, (T) n);
  pbbs::sequence<T> idx(n, [&](size_t i) { return r.ith_rand(i) % n; });
  auto f = [&](size_t i) {
    // putting write prefetch in slows it down
    //__builtin_prefetch (&out[idx[i+4]], 1, 1);
    pbbs::write_min(&out[idx[i]], (T)i, std::less<T>());
  };
  time(t, parallel_for(0, n - 4, f););
  return t;
}

template <typename T>
double t_shuffle(size_t n, bool check) {
  pbbs::sequence<T> in(n, [&](size_t i) { return i; });
  time(t, pbbs::random_shuffle(in, n););
  return t;
}

template <typename T>
bool check_histogram(pbbs::sequence<T> const &in,
                     pbbs::sequence<T> const &out) {
  size_t m = out.size();
  auto a = sort(in, std::less<T>());
  auto b = get_counts(a, [&](T x) { return x; }, m);
  size_t err_loc =
      pbbs::find_if_index(m, [&](size_t i) { return out[i] != b[i]; });
  if (err_loc != m) {
    cout << "ERROR in histogram at location " << err_loc << ", got "
         << out[err_loc] << ", expected " << b[err_loc] << endl;
    return false;
  }
  return true;
}

template <typename T>
double t_histogram(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i) % n; });
  pbbs::sequence<T> out;
  time(t, out = pbbs::histogram<T>(in, n););
  if (check) check_histogram(in, out);
  return t;
}

template <typename T>
double t_histogram_few(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i) % 256; });
  pbbs::sequence<T> out;
  time(t, out = pbbs::histogram<T>(in, 256););
  if (check) check_histogram(in, out);
  return t;
}

template <typename T>
double t_histogram_same(size_t n, bool check) {
  pbbs::sequence<T> in(n, (T)10311);
  pbbs::sequence<T> out;
  time(t, out = pbbs::histogram<T>(in, n););
  if (check) check_histogram(in, out);
  return t;
}

// this checks against a given sort
template <typename T, typename Cmp>
bool check_sort(pbbs::sequence<T> const &in, pbbs::sequence<T> const &out,
                Cmp less, std::string sort_name) {
  size_t n = in.size();
  auto a = pbbs::merge_sort(in, std::less<T>());
  size_t err_loc = pbbs::find_if_index(
      n, [&](size_t i) { return less(a[i], out[i]) || less(out[i], a[i]); });
  if (err_loc != n) {
    cout << "ERROR in " << sort_name << " at location " << err_loc << endl;
    return false;
  }
  return true;
}

template <typename T>
double t_sort(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i) % n; });
  pbbs::sequence<T> out;
  time(t, out = pbbs::sample_sort(in, std::less<T>()););
  if (check) check_sort(in, out, std::less<T>(), "sample sort");
  return t;
}

// no check since it is used for the sort for checking, and hence
// checked against the other sorts
template <typename T>
double t_merge_sort(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i) % n; });
  pbbs::sequence<T> out;
  time(t, pbbs::merge_sort_inplace(in.slice(), std::less<T>()););
  return t;
}

template <typename T>
double t_quicksort(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i) % n; });
  pbbs::sequence<T> copy;
  if (check) copy = in;
  time(t, pbbs::p_quicksort_inplace(in.slice(), std::less<T>()););
  if (check) check_sort(copy, in, std::less<T>(), "quicksort");
  return t;
}

template <typename T>
double t_count_sort_bits(size_t n, size_t bits) {
  pbbs::random r(0);
  size_t num_buckets = (1 << bits);
  size_t mask = num_buckets - 1;
  pbbs::sequence<T> in(n, [&](size_t i) { return r.ith_rand(i); });
  pbbs::sequence<T> out(n);
  auto f = [&](size_t i) { return in[i] & mask; };
  auto keys = pbbs::delayed_seq<unsigned char>(n, f);
  time(t, pbbs::count_sort(in, out.slice(), keys, num_buckets););
  for (size_t i = 1; i < n; i++) {
    if ((out[i - 1] & mask) > (out[i] & mask)) {
      cout << "ERROR in count sort at: " << i << endl;
      abort();
    }
  }
  return t;
}

template <typename T>
double t_count_sort_8(size_t n, bool check) {
  return t_count_sort_bits<T>(n, 8);
}

template <typename T>
double t_count_sort_2(size_t n, bool check) {
  return t_count_sort_bits<T>(n, 2);
}

template <typename T>
double t_collect_reduce_pair_dense(size_t n, bool check) {
  using par = std::pair<T, T>;
  pbbs::random r(0);
  pbbs::sequence<par> S(
      n, [&](size_t i) -> par { return par(r.ith_rand(i) % n, 1); });
  pbbs::sequence<T> out;
  auto get_key = [&](par a) { return a.first; };
  auto get_val = [&](par a) { return a.second; };
  time(t, out = pbbs::collect_reduce(S, get_key, get_val, pbbs::addm<T>(), n););
  if (check)
    check_histogram(pbbs::sequence<T>(n, [&](size_t i) { return S[i].first; }),
                    out);
  return t;
}

template <typename T>
double t_collect_reduce_pair_sparse(size_t n, bool check) {
  using par = std::pair<T, T>;
  pbbs::random r(0);
  struct hasheq {
    static inline size_t hash(par a) { return pbbs::hash64_2(a.first); }
    static inline bool eql(par a, par b) { return a.first == b.first; }
  };
  pbbs::sequence<par> S(
      n, [&](size_t i) -> par { return par(r.ith_rand(i) % n, 1); });
  time(t, pbbs::collect_reduce_sparse(S, hasheq(), pbbs::addm<T>()););
  return t;
}

template <typename T>
double t_collect_reduce_8(size_t n, bool check) {
  using par = std::pair<T, T>;
  pbbs::random r(0);
  size_t num_buckets = (1 << 8);
  pbbs::sequence<par> S(
      n, [&](size_t i) { return par(r.ith_rand(i) % num_buckets, 1); });
  auto get_key = [&](par a) { return a.first; };
  auto get_val = [&](par a) { return a.first; };
  time(t, pbbs::collect_reduce(S, get_key, get_val, pbbs::addm<T>(),
                               num_buckets););
  return t;
}

// template<typename T>
// double t_collect_reduce_8_tuple(size_t n, bool check) {
//   pbbs::random r(0);
//   size_t num_buckets = (1<<8);
//   size_t mask = num_buckets - 1;
//   using sums = std::tuple<float,float,float,float>;

//   auto bucket = [&] (size_t i) -> uchar { return r.ith_rand(i) & mask; };
//   auto keys = pbbs::delayed_seq<unsigned char>(n, bucket);

//   auto sum = [] (sums a, sums b) -> sums {
//     return sums(std::get<0>(a)+std::get<0>(b), std::get<1>(a)+std::get<1>(b),
// 		std::get<2>(a)+std::get<2>(b), std::get<3>(a)+std::get<3>(b));
//   };

//   pbbs::sequence<sums> in(n, [&] (size_t i) -> sums {
//       return sums(1.0,1.0,1.0,1.0);});

//   auto monoid = make_monoid(sum, sums(0.0,0.0,0.0,0.0));

//   time(t,
//        pbbs::collect_reduce<sums>(in, keys, num_buckets, monoid););
//   return t;
// }

template <typename T>
double t_integer_sort_pair(size_t n, bool check) {
  using par = std::pair<T, T>;
  pbbs::random r(0);
  size_t bits = sizeof(T) * 8;
  pbbs::sequence<par> S(n,
                        [&](size_t i) -> par { return par(r.ith_rand(i), i); });
  pbbs::sequence<par> R;
  auto first = [](par a) { return a.first; };
  time(t, R = pbbs::integer_sort(S.slice(), first, bits););

  auto less = [](par a, par b) { return a.first < b.first; };
  if (check) check_sort(S, R, less, "integer sort pair");
  return t;
}

template <typename T>
double t_integer_sort(size_t n, bool check) {
  pbbs::random r(0);
  size_t bits = sizeof(T) * 8;
  pbbs::sequence<T> S(n, [&](size_t i) -> T { return r.ith_rand(i); });
  auto identity = [](T a) { return a; };
  pbbs::sequence<T> R;
  // time(t, pbbs::integer_sort_inplace(S.slice(),identity,bits););
  time(t, R = pbbs::integer_sort(S, identity, bits););
  if (check) check_sort(S, R, std::less<T>(), "integer sort");
  return t;
}

double t_integer_sort_128(size_t n, bool check);

template <typename T>
double t_merge(size_t n, bool check) {
  pbbs::sequence<T> in1(n / 2, [&](size_t i) { return 2 * i; });
  pbbs::sequence<T> in2(n - n / 2, [&](size_t i) { return 2 * i + 1; });
  pbbs::sequence<T> out;
  time(t, out = pbbs::merge(in1, in2, std::less<T>()););
  return t;
}

template <typename T>
double t_remove_duplicates(size_t n, bool check) {
  pbbs::random r(0);
  pbbs::sequence<T> In(n, [&](size_t i) -> T { return r.ith_rand(i) % n; });
  time(t, pbbs::remove_duplicates(In););
  return t;
}

template <typename T, typename F>
static T my_reduce(pbbs::sequence<T> const &s, size_t start, size_t end, F f) {
  if (end - start == 1) return s[start];
  size_t h = (end + start) / 2;
  T r, l;
  auto left = [&]() { r = my_reduce(s, h, end, f); };
  auto right = [&]() { l = my_reduce(s, start, h, f); };
  par_do_if(h > 100, left, right);
  return f(l, r);
}

template <typename T>
double t_bag(size_t n, bool check) {
  using TB = pbbs::bag<T>;
  TB::init();
  pbbs::sequence<TB> In(n, [&](size_t i) -> TB { return TB((T)i); });
  time(t, TB x = my_reduce(In, 0, n, TB::append); x.flatten(););
  return t;
}

template <typename s_size_t, typename T>
double t_mat_vec_mult(size_t n, bool check) {
  pbbs::random r(0);
  size_t degree = 5;
  size_t m = degree * n;
  pbbs::sequence<s_size_t> starts(n + 1, [&](size_t i) { return degree * i; });
  pbbs::sequence<s_size_t> columns(m,
                                   [&](size_t i) { return r.ith_rand(i) % n; });
  pbbs::sequence<T> values(m, (T)1);
  pbbs::sequence<T> in(n, (T)1);
  pbbs::sequence<T> out(n, (T)0);
  auto add = [](T a, T b) { return a + b; };
  auto mult = [](T a, T b) { return a * b; };

  time(t, mat_vec_mult(starts, columns, values, in, out.slice(), mult, add););
  return t;
}

template <typename T>
double t_range_min(size_t n, bool check) {
  pbbs::sequence<T> In(n, [&](size_t i) { return 5; });
  In[n / 2] = 0;
  time(t, auto foo = pbbs::make_range_min(In, std::less<T>());
       parallel_for(0, n - 1, [&](size_t i) { foo.query(0, i); }););
  if (foo.query(0, n - 1) != n / 2 || foo.query(0, n / 2 - 1) != 0 ||
      foo.query(0, n / 2) != n / 2 ||
      foo.query(n / 2 + 1, n - 1) != n / 2 + 1) {
    cout << "error in range min query " << endl;
    abort();
  }
  return t;
}
