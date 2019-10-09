// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "WorkEfficientSDB14/Connectivity.h"
#include "UnionFind/Connectivity.h"
#include "ShiloachVishkin/Connectivity.h"

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

/* ************************* Connectivity wrappers *************************** */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);

template<typename Graph, typename F>
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<uintE>& correct, F test, commandLine& P);

template <class Graph>
double t_gbbs_cc(Graph& G, commandLine P, pbbs::sequence<uintE>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = workefficient_cc::CC(G, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_jayanti_cc(Graph& G, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto q = jayanti_rank::JayantiTBUnite<Graph>(G);
  auto CC = q.components(););
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_shiloach_vishkin_cc(Graph& G, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto CC = shiloachvishkin_cc::CC(G););
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

/* ************************* Benchmark Utils *************************** */

template<typename F>
double reduce(std::vector<double> V, F f) {
  double x = V[0];
  for (size_t i=1; i < V.size(); i++) x = f(x,V[i]);
  return x;
}

double median(std::vector<double> V) {
  std::sort(V.begin(),V.end());
  if (V.size()%2 == 1)
    return V[V.size()/2];
  else
    return (V[V.size()/2] + V[V.size()/2 - 1])/2.0;
}

double sumf(double a, double b) {return a+ b;};
double minf(double a, double b) {return (a < b) ? a : b;};
double maxf(double a, double b) {return (a > b) ? a : b;};

template<typename Graph, typename F>
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<uintE>& correct, F test, commandLine& P) {
  std::vector<double> R;
  for (size_t i=0; i < rounds; i++) {
    R.push_back(test(G, P, correct));
  }
  return R;
}

template<typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, pbbs::sequence<uintE>& correct,
		  std::string name, commandLine& P, F test) {
  std::vector<double> t = repeat(G, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  cout << name << std::setprecision(5)
       << ": r=" << rounds
       << ", med=" << med
       << " (" << mint << "," << maxt << "), "
       << endl;
  return 1;
}

template<typename Graph,
          template <class F, class U, class G> class UFSample>
bool run_multiple_uf_apply_unite(Graph& G, size_t rounds,
    pbbs::sequence<uintE>& correct, std::string sample, std::string unite,
    std::string find, commandLine& P) {
  std::vector<double> t;

  if (unite == "unite") {
    auto test = [&] (Graph& G, commandLine& P, pbbs::sequence<uintE>& correct) -> double {
      size_t sampling_rounds = P.getOptionLongValue("-sample_rounds", 2L);
      timer t; t.start();
      auto components =
          union_find::select_algorithm<
          unite_variants::Unite,
          UFSample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
      double t_out = t.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, components);
      }
      return t_out;
    };
    t = repeat(G, rounds, correct, test, P);
  } else if (unite == "unite_early") {
    auto test = [&] (Graph& G, commandLine& P, pbbs::sequence<uintE>& correct) -> double {
      size_t sampling_rounds = P.getOptionLongValue("-sample_rounds", 2L);
      timer t; t.start();
      auto components =
          union_find::select_algorithm<
          unite_variants::UniteEarly,
          UFSample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
      double t_out = t.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, components);
      }
      return t_out;
    };
    t = repeat(G, rounds, correct, test, P);
  } else if (unite == "unite_nd") {
    auto test = [&] (Graph& G, commandLine& P, pbbs::sequence<uintE>& correct) -> double {
      size_t sampling_rounds = P.getOptionLongValue("-sample_rounds", 2L);
      timer t; t.start();
      auto components =
          union_find::select_algorithm<
          unite_variants::UniteND,
          UFSample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */true);
      double t_out = t.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, components);
      }
      return t_out;
    };
    t = repeat(G, rounds, correct, test, P);
  } else if (unite == "unite_rem") {
    auto test = [&] (Graph& G, commandLine& P, pbbs::sequence<uintE>& correct) -> double {
      size_t sampling_rounds = P.getOptionLongValue("-sample_rounds", 2L);
      timer t; t.start();
      auto components =
          union_find::select_algorithm<
          unite_variants::UniteRem,
          UFSample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
      double t_out = t.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, components);
      }
      return t_out;
    };
    t = repeat(G, rounds, correct, test, P);
  } else {
    std::cout << "Unknown unite argument: " << unite << std::endl;
    exit(0);
  }

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  auto name = "union_find; sample = " + sample + "; unite = " + unite + "; find = " + find;

  cout << name << std::setprecision(5)
       << ": r=" << rounds
       << ", med=" << med
       << " (" << mint << "," << maxt << "), "
       << endl;
  return 1;
}

template<typename Graph>
bool run_multiple_uf(Graph& G, size_t rounds,
    pbbs::sequence<uintE>& correct, std::string sample, std::string unite,
    std::string find, commandLine& P) {
  std::vector<double> t;

  if (sample == "kout") {
    return run_multiple_uf_apply_unite<Graph, union_find::UnionFindSampleTemplate>(G, rounds, correct, sample, unite, find, P);
  } else if (sample == "bfs") {
    return run_multiple_uf_apply_unite<Graph, union_find::UnionFindSampledBFSTemplate>(G, rounds, correct, sample, unite, find, P);
  } else if (sample == "ldd") {
    return run_multiple_uf_apply_unite<Graph, union_find::UnionFindLDDTemplate>(G, rounds, correct, sample, unite, find, P);
  } else if (sample == "none") {
    return run_multiple_uf_apply_unite<Graph, union_find::UnionFindTemplate>(G, rounds, correct, sample, unite, find, P);
  } else {
    std::cout << "Unknown sampling argument: " << sample << std::endl;
    exit(0);
  }
}



template <class Graph>
double pick_test(Graph& G, size_t id, size_t rounds, commandLine P, pbbs::sequence<uintE>& correct) {
  switch (id) {
  case 0:
    return run_multiple(G, rounds, correct, "gbbs_cc", P, t_gbbs_cc<Graph>);

  case 1:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_compress", P);
  case 2:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_naive", P);
  case 3:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_split", P);
  case 4:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_halve", P);
  case 5:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_atomic_split", P);
  case 6:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite", /* find = */ "find_atomic_halve", P);

  case 7:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_compress", P);
  case 8:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_naive", P);
  case 9:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_split", P);
  case 10:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_halve", P);
  case 11:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_atomic_split", P);
  case 12:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_nd", /* find = */ "find_atomic_halve", P);

  case 13:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_compress", P);
  case 14:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_naive", P);
  case 15:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_split", P);
  case 16:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_halve", P);
  case 17:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_atomic_split", P);
  case 18:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_early", /* find = */ "find_atomic_halve", P);

  case 19:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_compress", P);
  case 20:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_naive", P);
  case 21:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_split", P);
  case 22:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_halve", P);
  case 23:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_atomic_split", P);
  case 24:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "kout", /* unite = */ "unite_rem", /* find = */ "find_atomic_halve", P);

  case 25:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_compress", P);
  case 26:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_naive", P);
  case 27:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_split", P);
  case 28:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_halve", P);
  case 29:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_atomic_split", P);
  case 30:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite", /* find = */ "find_atomic_halve", P);

  case 31:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_compress", P);
  case 32:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_naive", P);
  case 33:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_split", P);
  case 34:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_halve", P);
  case 35:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_atomic_split", P);
  case 36:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_nd", /* find = */ "find_atomic_halve", P);

  case 37:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_compress", P);
  case 38:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_naive", P);
  case 39:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_split", P);
  case 40:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_halve", P);
  case 41:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_atomic_split", P);
  case 42:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_early", /* find = */ "find_atomic_halve", P);

  case 43:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_compress", P);
  case 44:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_naive", P);
  case 45:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_split", P);
  case 46:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_halve", P);
  case 47:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_atomic_split", P);
  case 48:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "bfs", /* unite = */ "unite_rem", /* find = */ "find_atomic_halve", P);

  case 49:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_compress", P);
  case 50:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_naive", P);
  case 51:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_split", P);
  case 52:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_halve", P);
  case 53:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_atomic_split", P);
  case 54:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite", /* find = */ "find_atomic_halve", P);

  case 55:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_compress", P);
  case 56:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_naive", P);
  case 57:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_split", P);
  case 58:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_halve", P);
  case 59:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_atomic_split", P);
  case 60:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_nd", /* find = */ "find_atomic_halve", P);

  case 61:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_compress", P);
  case 62:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_naive", P);
  case 63:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_split", P);
  case 64:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_halve", P);
  case 65:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_atomic_split", P);
  case 66:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_early", /* find = */ "find_atomic_halve", P);

  case 67:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_compress", P);
  case 68:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_naive", P);
  case 69:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_split", P);
  case 70:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_halve", P);
  case 71:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_atomic_split", P);
  case 72:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "ldd", /* unite = */ "unite_rem", /* find = */ "find_atomic_halve", P);

  case 73:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_compress", P);
  case 74:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_naive", P);
  case 75:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_split", P);
  case 76:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_halve", P);
  case 77:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_atomic_split", P);
  case 78:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite", /* find = */ "find_atomic_halve", P);

  case 79:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_compress", P);
  case 80:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_naive", P);
  case 81:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_split", P);
  case 82:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_halve", P);
  case 83:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_atomic_split", P);
  case 84:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_nd", /* find = */ "find_atomic_halve", P);

  case 85:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_compress", P);
  case 86:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_naive", P);
  case 87:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_split", P);
  case 88:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_halve", P);
  case 89:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_atomic_split", P);
  case 90:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_early", /* find = */ "find_atomic_halve", P);

  case 91:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_compress", P);
  case 92:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_naive", P);
  case 93:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_split", P);
  case 94:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_halve", P);
  case 95:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_atomic_split", P);
  case 96:
    return run_multiple_uf(G, rounds, correct, /* sample = */ "none", /* unite = */ "unite_rem", /* find = */ "find_atomic_halve", P);

  case 97: /* Jayanti */
    return run_multiple(G, rounds, correct, "union_find: jayanti", P, t_jayanti_cc<Graph>);

  case 98: /* Shiloach-Vishkin */
    return run_multiple(G, rounds, correct, "shiloach-vishkin", P, t_shiloach_vishkin_cc<Graph>);


  default:
    assert(false);
    std::cout << "Unknown test" << std::endl;
    exit(-1);
    return 0.0 ;
  }
}

/* ************************* Utils *************************** */

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  pbbslib::scan_add_inplace(flags);
  std::cout << "n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = pbbslib::reduce_max(flags);
  std::cout << "largest_cc has size: " << sz << "\n";
  return sz;
}

template <class Seq>
inline size_t RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = pbbs::sequence<T>(n + 1, (T)0);
  T cur_comp = 1;
  for (size_t i=0; i<n; i++) {
    T comp = ids[i];
    T new_comp = cur_comp++;
    if (component_map[comp] == 0) {
      component_map[comp] = new_comp;
    }
    ids[i] = new_comp;
  }
  return cur_comp;
}

template <class S1, class S2>
inline void cc_check(S1& correct, S2& check) {
  RelabelDet(check);

  bool is_correct = true;
  uintE max_cor = 0;
  uintE max_chk = 0;
  for (size_t i=0; i<correct.size(); i++) {
    assert(correct[i] == check[i]);
    if ((correct[i] != check[i])) {
      is_correct = false;
      cout << "at i = " << i << " cor = " << correct[i] << " got: " << check[i] << endl;
    }
    if (correct[i] > max_cor) {
      pbbs::write_max(&max_cor, correct[i], std::less<uintE>());
    }
    if (check[i] > max_chk) {
      pbbs::write_max(&max_chk, check[i], std::less<uintE>());
    }
  }//);
  cout << "correctness check: " << is_correct << endl;
  cout << "max_cor = " << max_cor << " max_chk = " << max_chk << endl;
}

/* ************************* ***** *************************** */

template <class Graph>
double Benchmark_runner(Graph& G, commandLine P) {
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);
  bool symmetric = P.getOptionValue("-s");
  int num_tests = 99; // update if new algorithm is added

  cout << "rounds = " << rounds << endl;
  cout << "num threads = " << num_workers() << endl;

  auto correct = pbbs::sequence<uintE>();
  if (P.getOptionValue("-check")) {
    correct = workefficient_cc::CC(G, 0.2, false, true);
    RelabelDet(correct);
  }
  if (test_num == -1) {
    for (int i=0; i < num_tests; i++) {
      pick_test(G, i, rounds, P, correct);
    }
  } else {
    pick_test(G, test_num, rounds, P, correct);
  }
  return 1.0;
}

generate_symmetric_once_main(Benchmark_runner, false);
