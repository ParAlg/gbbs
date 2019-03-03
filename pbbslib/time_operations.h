#include "utilities.h"
#include "get_time.h"
#include "random.h"
#include "counting_sort.h"
#include "collect_reduce.h"
#include "random_shuffle.h"
#include "histogram.h"
#include "integer_sort.h"
#include "sample_sort.h"
#include "merge.h"
#include "merge_sort.h"
#include "bag.h"
#include "hash_table.h"
#include "sparse_mat_vec_mult.h"
#include "sequence_ops.h"
#include "monoid.h"

#include <iostream>
#include <ctype.h>
#include <math.h>
#include <assert.h>

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

template<typename T>
double t_tabulate(size_t n) {
  auto f = [] (size_t i) {return i;};
  time(t, pbbslib::sequence<T>(n, f););
  return t;
}

// template<typename T>
// double t_map(size_t n) {
//   size_t l = 1 << 14;
//   size_t ll = 1 << 10;
//   sequence<T> In(l*l, (T) 1);
//   sequence<T> Out(l*l, (T) 1);
//   auto f = [&] (size_t i) {
//     size_t r = ll*(i/16);
//     size_t c = ll*(i%16);
//     for (size_t j=r+1; j < r + ll -1; j++) {
//       for (size_t k=c+1; k < c+ ll -1; k++) {
// 	Out[j + k*l] = In[j + k*l - 1 - ll] + In[j + k*l - ll] + In[j + k*l + 1 - ll] +
// 	  In[j + k*l - 1] + In[j + k*l] + In[j + k*l + 1] +
// 	  In[j + k*l - 1 + ll] + In[j + k*l + ll] + In[j + k*l + 1 + ll];
//       }
//     }
//   };
//   time(t, parallel_for(0, 256, f););
//   return t;
// }


template<typename T>
double t_map(size_t n) {
  pbbslib::sequence<T> In(n, (T) 1);
  auto f = [&] (size_t i) {return In[i];};
  time(t, pbbslib::sequence<T>(n, f));
  return t;
}

template<typename T>
double t_reduce_add(size_t n) {
  pbbslib::sequence<T> S(n, (T) 1);
  //time(t, sum(S.begin(), S.size()););
  //time(t, pbbslib::reduce_add(S););
  time(t, pbbslib::reduce(S, pbbslib::addm<T>()));
  return t;
}

double t_map_reduce_128(size_t n) {
  int stride = 16;
  pbbslib::sequence<size_t> S(n*stride, (size_t) 1);
  auto get = [&] (size_t i) {
    // gives marginal improvement (5% or so on aware)
    __builtin_prefetch (&S[(i+4)*stride], 0, 3);
    return S[i*stride];};
  auto T = pbbslib::delayed_seq<size_t>(n, get);
  time(t, pbbslib::reduce(T, pbbslib::addm<size_t>()););
  return t;
}

template<typename T>
double t_scan_add(size_t n) {
  pbbslib::sequence<T> In(n, (T) 1);
  pbbslib::sequence<T> Out;
  T sum;
  time(t, std::tie(Out,sum) = pbbslib::scan(In, pbbslib::addm<T>()););
  return t;
}

template<typename T>
double t_scan_add_seq(size_t n) {
  pbbslib::sequence<T> In(n, (T) 1);
  pbbslib::sequence<T> Out(n, (T) 0);
  time(t,
       T x = 0;
       for (size_t i=0; i < n; i++) {
	 Out[i] = x;
	 x += In[i];
       });
  return t;
}

template<typename T>
double t_pack(size_t n) {
  pbbslib::sequence<bool> flags(n, [] (size_t i) -> bool {return i%2;});
  pbbslib::sequence<T> In(n, [] (size_t i) -> T {return i;});
  time(t, pbbslib::pack(In, flags););
  return t;
}

template<typename T>
double t_split3_old(size_t n) {
  pbbslib::sequence<uchar> flags(n, [] (size_t i) -> uchar {return i%3;});
  pbbslib::sequence<T> In(n, [] (size_t i) -> T {return i;});
  pbbslib::sequence<T> Out(n, (T) 0);
  time(t, pbbslib::split_three(In, Out, flags););
  return t;
}

template<typename T>
double t_split3(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> In(n, [&] (size_t i) {return r.ith_rand(i);});
  pbbslib::sequence<T> Out(n, (T) 0);
  time(t, pbbslib::p_split3(In, Out.slice(), std::less<T>()););
  return t;
}

template<typename T>
double t_partition(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> In(n, [&] (size_t i) -> size_t {return r.ith_rand(i)%n;});
  pbbslib::sequence<T> Out;
  auto f = pbbslib::delayed_seq<bool>(n, [&] (size_t i) {return In[i] < n/((size_t) 2);});
  size_t m;
  time(t, std::tie(Out,m) = pbbslib::split_two(In, f););
  return t;
}

#if defined(CILK)
#include "reducer.h"
double t_histogram_reducer(size_t n) {
  pbbslib::random r(0);
  constexpr int count = 1024;
  histogram_reducer<int,count> red;
  using aa = std::array<size_t,8>;
  pbbslib::sequence<aa> In(n, [&] (size_t i) {aa x; x[0] = r.ith_rand(i) % count; return x;});
  auto f = [&] (size_t i) { red->add_value(In[i][0]);};
  time(t, parallel_for(0, n, f););
  //cout << red.get_value()[0] << endl;
  return t;
}
#else
double t_histogram_reducer(size_t n) {
  return 1.0;
}
#endif

template<typename T>
double t_gather(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return i;});
  pbbslib::sequence<T> idx(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  auto f = [&] (size_t i) {
    // prefetching helps significantly
    __builtin_prefetch (&in[idx[i+4]], 0, 1);
    return in[idx[i]];};
  // note problem with prefetching since will go over end for last 4 iterations
  time(t, pbbslib::sequence<T>(n-4, f););  
  return t;
}

template<typename T>
double t_scatter(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> out(n, (T) 0);
  pbbslib::sequence<T> idx(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  auto f = [&] (size_t i) {
    // prefetching makes little if any difference
    //__builtin_prefetch (&out[idx[i+4]], 1, 1);
      out[idx[i]] = i;};
  time(t, parallel_for(0, n-4, f););
  return t;
}

template<typename T>
double t_write_add(size_t n) {
  pbbslib::random r(0);
  //pbbslib::sequence<T> out(n, (T) 0);
  pbbslib::sequence<std::atomic<T>> out(n);
  parallel_for(0,n,[&] (size_t i) {std::atomic_init(&out[i], (T) 0);});
  pbbslib::sequence<T> idx(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  auto f = [&] (size_t i) {
    // putting write prefetch in slows it down
    //__builtin_prefetch (&out[idx[i+4]], 0, 1);
    //__sync_fetch_and_add(&out[idx[i]],1);};
    pbbslib::write_add(&out[idx[i]],1);};
  time(t, parallel_for(0, n-4, f););
  return t;
}

template<typename T>
double t_write_min(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<std::atomic<T>> out(n);
  parallel_for(0,n,[&] (size_t i) {std::atomic_init(&out[i], (T) n);});
  //pbbslib::sequence<T> out(n, (T) n);
  pbbslib::sequence<T> idx(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  auto f = [&] (size_t i) {
    // putting write prefetch in slows it down
    //__builtin_prefetch (&out[idx[i+4]], 1, 1);
    pbbslib::write_min(&out[idx[i]], (T) i, std::less<T>());};
  time(t, parallel_for(0, n-4, f););
  return t;
}

template<typename T>
double t_shuffle(size_t n) {
  pbbslib::sequence<T> in(n, [&] (size_t i) {return i;});
  time(t, pbbslib::random_shuffle(in,n););
  return t;
}

template<typename T>
double t_histogram(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  pbbslib::sequence<T> out;
  time(t, out = pbbslib::histogram<T>(in,n););
  return t;
}

template<typename T>
double t_histogram_few(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i)%256;});
  pbbslib::sequence<T> out;
  time(t, out = pbbslib::histogram<T>(in,256););
  return t;
}

template<typename T>
double t_histogram_same(size_t n) {
  pbbslib::sequence<T> in(n, (T) 10311);
  pbbslib::sequence<T> out;
  time(t, out = pbbslib::histogram<T>(in,n););
  return t;
}

template<typename T>
double t_sort(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  pbbslib::sequence<T> out;
  time(t, out = pbbslib::sample_sort(in, std::less<T>()););
  //for (size_t i = 1; i < n; i++)
  //  if (std::less<T>()(in[i],in[i-1])) {cout << i << endl; abort();}
  return t;
}

template<typename T>
double t_merge_sort(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  time(t, pbbslib::merge_sort_inplace(in.slice(), std::less<T>()););
  //for (size_t i = 1; i < n; i++)
  //  if (std::less<T>()(in[i],in[i-1])) {cout << i << endl; abort();}
  return t;
}

template<typename T>
double t_quicksort(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i)%n;});
  time(t, pbbslib::p_quicksort_inplace(in.slice(), std::less<T>()););
  //for (size_t i = 1; i < n; i++)
  //  if (std::less<T>()(in[i],in[i-1])) {cout << i << endl; abort();}
  return t;
}

template<typename T>
double t_count_sort_bits(size_t n, size_t bits) {
  pbbslib::random r(0);
  size_t num_buckets = (1<<bits);
  size_t mask = num_buckets - 1;
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i);});
  pbbslib::sequence<T> out(n);
  auto f = [&] (size_t i) {return in[i] & mask;};
  auto keys = pbbslib::delayed_seq<unsigned char>(n, f);
  time(t, pbbslib::count_sort(in, out.slice(), keys, num_buckets););
  for (size_t i=1; i < n; i++) {
    if ((out[i-1] & mask) > (out[i] & mask)) {
      cout << "error in count sort at: " << i << endl;
      abort();
    }
  }
  return t;
}

template<typename T>
double t_count_sort_8(size_t n) {return t_count_sort_bits<T>(n, 8);}

template<typename T>
double t_count_sort_2(size_t n) {return t_count_sort_bits<T>(n, 2);}

template<typename T>
double t_collect_reduce_pair_dense(size_t n) {
  using par = std::pair<T,T>;
  pbbslib::random r(0);
  pbbslib::sequence<par> S(n, [&] (size_t i) -> par {
      return par(r.ith_rand(i) % n, 1);});
  auto get_index = [] (par e) {return e.first;};
  auto get_val = [] (par e) {return e.second;};
  time(t, pbbslib::collect_reduce<T>(S, n, get_index, get_val, pbbslib::addm<T>()););
  return t;
}

template<typename T>
double t_collect_reduce_pair_sparse(size_t n) {
  using par = std::pair<T,T>;
  pbbslib::random r(0);
  pbbslib::sequence<par> S(n, [&] (size_t i) -> par {
      return par(r.ith_rand(i) % n, 1);});
  time(t, pbbslib::collect_reduce_pair(S, pbbslib::addm<T>()););
  return t;
}

template<typename T>
double t_collect_reduce_8(size_t n) {
  pbbslib::random r(0);
  size_t num_buckets = (1<<8);
  size_t mask = num_buckets - 1;
  pbbslib::sequence<T> in(n, [&] (size_t i) {return r.ith_rand(i);});
  auto bucket = [&] (size_t i) {return in[i] & mask;};
  auto keys = pbbslib::delayed_seq<unsigned char>(n, bucket);
  time(t, pbbslib::collect_reduce<T>(in, keys, num_buckets, pbbslib::addm<T>()););
  return t;
}

template<typename T>
double t_collect_reduce_8_tuple(size_t n) {
  pbbslib::random r(0);
  size_t num_buckets = (1<<8);
  size_t mask = num_buckets - 1;
  using sums = std::tuple<float,float,float,float>;

  auto bucket = [&] (size_t i) -> uchar { return r.ith_rand(i) & mask; };
  auto keys = pbbslib::delayed_seq<unsigned char>(n, bucket);

  auto sum = [] (sums a, sums b) -> sums {
    return sums(std::get<0>(a)+std::get<0>(b), std::get<1>(a)+std::get<1>(b),
		std::get<2>(a)+std::get<2>(b), std::get<3>(a)+std::get<3>(b));
  };

  pbbslib::sequence<sums> in(n, [&] (size_t i) -> sums {
      return sums(1.0,1.0,1.0,1.0);});

  auto monoid = make_monoid(sum, sums(0.0,0.0,0.0,0.0));

  time(t,
       pbbslib::collect_reduce<sums>(in, keys, num_buckets, monoid););
  return t;
}


template<typename T>
double t_integer_sort_pair(size_t n) {
  using par = std::pair<T,T>;
  pbbslib::random r(0);
  size_t bits = sizeof(T)*8;
  pbbslib::sequence<par> S(n, [&] (size_t i) -> par {
      return par(r.ith_rand(i),i);});
  pbbslib::sequence<par> R;
  auto first = [] (par a) {return a.first;};
  time(t, R = pbbslib::integer_sort(S.slice(),first,bits););
  return t;
}

template<typename T>
double t_integer_sort(size_t n) {
  pbbslib::random r(0);
  size_t bits = sizeof(T)*8;
  pbbslib::sequence<T> S(n, [&] (size_t i) -> T {
      return r.ith_rand(i);});
  auto identity = [] (T a) {return a;};
  time(t, pbbslib::integer_sort_inplace(S.slice(),identity,bits););
  return t;
}

typedef unsigned __int128 long_int;
double t_integer_sort_128(size_t n) {
  pbbslib::random r(0);
  size_t bits = 128;
  pbbslib::sequence<long_int> S(n, [&] (size_t i) -> long_int {
      return r.ith_rand(2*i) + (((long_int) r.ith_rand(2*i+1)) << 64) ;});
  auto identity = [] (long_int a) {return a;};
  time(t, pbbslib::integer_sort_inplace(S.slice(),identity,bits););
  return t;
}

template<typename T>
double t_merge(size_t n) {
  pbbslib::sequence<T> in1(n/2, [&] (size_t i) {return 2*i;});
  pbbslib::sequence<T> in2(n-n/2, [&] (size_t i) {return 2*i+1;});
  pbbslib::sequence<T> out;
  time(t, out = pbbslib::merge(in1, in2, std::less<T>()););
  return t;
}

template<typename T>
double t_remove_duplicates(size_t n) {
  pbbslib::random r(0);
  pbbslib::sequence<T> In(n, [&] (size_t i) -> T {return r.ith_rand(i) % n;});
  time(t, pbbslib::remove_duplicates(In););
  return t;
}

template <typename T, typename F>
static T my_reduce(pbbslib::sequence<T> const &s, size_t start, size_t end, F f) {
  if (end - start == 1) return s[start];
  size_t h = (end + start)/2;
  T r, l;
  auto left = [&] () {r = my_reduce(s, h, end, f);};
  auto right = [&] () {l = my_reduce(s, start, h, f);};
  par_do_if(h > 100, left, right);
  return f(l,r);
}

template<typename T>
double t_bag(size_t n) {
  using TB = pbbslib::bag<T>;
  TB::init();
  pbbslib::sequence<TB> In(n, [&] (size_t i) -> TB {return TB((T) i);});
  time(t, TB x = my_reduce(In, 0, n, TB::append); x.flatten(););
  return t;
}

template<typename s_size_t, typename T>
double t_mat_vec_mult(size_t n) {
  pbbslib::random r(0);
  size_t degree = 5;
  size_t m = degree*n;
  pbbslib::sequence<s_size_t> starts(n+1, [&] (size_t i) {
      return degree*i;});
  pbbslib::sequence<s_size_t> columns(m, [&] (size_t i) {
      return r.ith_rand(i)%n;});
  pbbslib::sequence<T> values(m, (T) 1);
  pbbslib::sequence<T> in(n, (T) 1);
  pbbslib::sequence<T> out(n, (T) 0);
  auto add = [] (T a, T b) { return a + b;};
  auto mult = [] (T a, T b) { return a * b;};

  time(t, mat_vec_mult(starts, columns, values, in, out.slice(), mult, add););
  return t;
}

