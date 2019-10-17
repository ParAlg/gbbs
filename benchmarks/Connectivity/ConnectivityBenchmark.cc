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
#include "LabelPropagation/Connectivity.h"
#include "BFSCC/Connectivity.h"
#include "Framework/framework.h"

#include "common.h"

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
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<parent>& correct, F test, commandLine& P);

template <class Graph>
double t_gbbs_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = workefficient_cc::CC(G, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_bfs_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  time(t, auto CC = bfs_cc::CC(G));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_shiloach_vishkin_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  time(t, auto CC = shiloachvishkin_cc::CC(G););
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_sample_only_algorithm_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  time(t, auto CC = labelprop_cc::CC</*use_permutation = */false>(G););
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
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<parent>& correct, F test, commandLine& P) {
  std::vector<double> R;
  for (size_t i=0; i < rounds; i++) {
    auto t = test(G, P, correct);
    std::cout << "### t = " << t << std::endl;
    R.push_back(t);
  }
  return R;
}

template<typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, pbbs::sequence<parent>& correct,
		  std::string name, commandLine& P, F test) {
  std::vector<double> t = repeat(G, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  cout << "Test = {"
       << "\"name\": \"" << name << "\"" << std::setprecision(5)
       << ", \"rounds\":" << rounds
       << ", \"med\":" << med
       << ", \"min\":" << mint
       << ", \"max\":" << maxt << "}"
       << endl;
  return 1;
}



namespace connectit {

  template<
    class Graph,
    SamplingOption sampling_option,
    UniteOption    unite_option,
    FindOption     find_option>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_uf_alg<
            Graph,
            sampling_option,
            find_option,
            unite_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = uf_options_to_string<sampling_option, find_option, unite_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption sampling_option,
    UniteOption    unite_option,
    FindOption     find_option,
    SpliceOption   splice_option>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_uf_alg<
            Graph,
            sampling_option,
            find_option,
            unite_option,
            splice_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = uf_options_to_string<sampling_option, find_option, unite_option, splice_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption    sampling_option,
    JayantiFindOption find_option>
  bool run_multiple_jayanti_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_jayanti_alg<
            Graph,
            sampling_option,
            find_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = jayanti_options_to_string<sampling_option, find_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }


  template<
    class Graph,
    SamplingOption          sampling_option,
    LiuTarjanConnectOption  connect_option,
    LiuTarjanUpdateOption   update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption    alter_option>
  bool run_multiple_liu_tarjan_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_liu_tarjan_alg<
            Graph,
            sampling_option,
            connect_option,
            update_option,
            shortcut_option,
            alter_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = liu_tarjan_options_to_string<sampling_option, connect_option, update_option, shortcut_option, alter_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption    sampling_option,
    template <class G> class Algorithm>
  bool run_multiple_sample_only_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P,
      std::string name) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC = run_sample_only_alg<Graph, sampling_option, Algorithm>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto test_name = name + "; " + sampling_to_string<sampling_option>();
    return run_multiple(G, rounds, correct, test_name, P, test);
  }


  template <class Graph>
  double pick_test(Graph& G, size_t id, size_t rounds, commandLine P, pbbs::sequence<parent>& correct) {

    switch (id) {
      case 0:
        return run_multiple(G, rounds, correct, "gbbs_cc", P, t_gbbs_cc<Graph>);

      case 1:
        return run_multiple_uf_alg<Graph, no_sampling, unite, find_compress>(G, rounds, correct, P);
      case 2:
        return run_multiple_uf_alg<Graph, no_sampling, unite, find_naive>(G, 1, correct, P);
      case 3:
        return run_multiple_uf_alg<Graph, no_sampling, unite, find_atomic_split>(G, rounds, correct, P);
      case 4:
        return run_multiple_uf_alg<Graph, no_sampling, unite, find_atomic_halve>(G, rounds, correct, P);

      case 5:
        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_compress>(G, rounds, correct, P);
      case 6:
        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_naive>(G, 1, correct, P);
      case 7:
        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_atomic_split>(G, rounds, correct, P);
      case 8:
        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_atomic_halve>(G, rounds, correct, P);

      case 9:
        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_compress>(G, rounds, correct, P);
      case 10:
        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_naive>(G, 1, correct, P);
      case 11:
        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_atomic_split>(G, rounds, correct, P);
      case 12:
        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_atomic_halve>(G, rounds, correct, P);

      case 13:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_compress>(G, rounds, correct, P);
      case 14:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_naive>(G, 1, correct, P);
      case 15:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
      case 16:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);

      case 17:
        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_compress>(G, rounds, correct, P);
      case 18:
        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_naive>(G, 1, correct, P);
      case 19:
        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_split>(G, rounds, correct, P);
      case 20:
        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_halve>(G, rounds, correct, P);

      case 21:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_compress>(G, rounds, correct, P);
      case 22:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_naive>(G, 1, correct, P);
      case 23:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_atomic_split>(G, rounds, correct, P);
      case 24:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_atomic_halve>(G, rounds, correct, P);

      case 25:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_compress>(G, rounds, correct, P);
      case 26:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_naive>(G, 1, correct, P);
      case 27:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_atomic_split>(G, rounds, correct, P);
      case 28:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_atomic_halve>(G, rounds, correct, P);

      case 29:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_compress>(G, rounds, correct, P);
      case 30:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_naive>(G, 1, correct, P);
      case 31:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
      case 32:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);

      case 33:
        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_compress>(G, rounds, correct, P);
      case 34:
        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_naive>(G, 1, correct, P);
      case 35:
        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_atomic_split>(G, rounds, correct, P);
      case 36:
        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_atomic_halve>(G, rounds, correct, P);

      case 37:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_compress>(G, rounds, correct, P);
      case 38:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_naive>(G, 1, correct, P);
      case 39:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_atomic_split>(G, rounds, correct, P);
      case 40:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_atomic_halve>(G, rounds, correct, P);

      case 41:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_compress>(G, rounds, correct, P);
      case 42:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_naive>(G, 1, correct, P);
      case 43:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_atomic_split>(G, rounds, correct, P);
      case 44:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_atomic_halve>(G, rounds, correct, P);

      case 45:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_compress>(G, rounds, correct, P);
      case 46:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_naive>(G, 1, correct, P);
      case 47:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
      case 48:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);

      case 49:
        return run_multiple_uf_alg<Graph, sample_kout, unite, find_compress>(G, rounds, correct, P);
      case 50:
        return run_multiple_uf_alg<Graph, sample_kout, unite, find_naive>(G, 1, correct, P);
      case 51:
        return run_multiple_uf_alg<Graph, sample_kout, unite, find_atomic_split>(G, rounds, correct, P);
      case 52:
        return run_multiple_uf_alg<Graph, sample_kout, unite, find_atomic_halve>(G, rounds, correct, P);

      case 53:
        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_compress>(G, rounds, correct, P);
      case 54:
        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_naive>(G, 1, correct, P);
      case 55:
        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_atomic_split>(G, rounds, correct, P);
      case 56:
        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_atomic_halve>(G, rounds, correct, P);

      case 57:
        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_compress>(G, rounds, correct, P);
      case 58:
        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_naive>(G, 1, correct, P);
      case 59:
        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_atomic_split>(G, rounds, correct, P);
      case 60:
        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_atomic_halve>(G, rounds, correct, P);

      case 61:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_compress>(G, rounds, correct, P);
      case 62:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_naive>(G, 1, correct, P);
      case 63:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
      case 64:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);

      /* Jayanti strategies */
      case 65:
        return run_multiple_jayanti_alg<Graph, sample_kout, find_twotrysplit>(G, rounds, correct,  P);
      case 66:
        return run_multiple_jayanti_alg<Graph, sample_kout, find_simple>(G, rounds, correct,  P);
      case 67:
        return run_multiple_jayanti_alg<Graph, sample_bfs, find_twotrysplit>(G, rounds, correct,  P);
      case 68:
        return run_multiple_jayanti_alg<Graph, sample_bfs, find_simple>(G, rounds, correct,  P);
      case 69:
        return run_multiple_jayanti_alg<Graph, sample_ldd, find_twotrysplit>(G, rounds, correct,  P);
      case 70:
        return run_multiple_jayanti_alg<Graph, sample_ldd, find_simple>(G, rounds, correct,  P);
      case 71:
        return run_multiple_jayanti_alg<Graph, no_sampling, find_twotrysplit>(G, rounds, correct,  P);
      case 72:
        return run_multiple_jayanti_alg<Graph, no_sampling, find_simple>(G, rounds, correct,  P);

      /* Label Propagation strategies */
      case 73:
        return run_multiple_sample_only_alg<Graph, sample_kout, labelprop_cc::LPAlgorithm>(G, rounds, correct, P, "label_prop");
      case 74:
        return run_multiple_sample_only_alg<Graph, sample_bfs, labelprop_cc::LPAlgorithm>(G, rounds, correct, P, "label_prop");
      case 75:
        return run_multiple_sample_only_alg<Graph, sample_ldd, labelprop_cc::LPAlgorithm>(G, rounds, correct, P, "label_prop");
      case 76:
        return run_multiple_sample_only_alg<Graph, no_sampling, labelprop_cc::LPAlgorithm>(G, rounds, correct, P, "label_prop");

      /* UF Rem-CAS strategies */
      case 77:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_compress, split_atomic_one>(G, rounds, correct, P);
      case 78:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_compress, halve_atomic_one>(G, rounds, correct, P);
      case 79:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_compress, splice_simple>(G, rounds, correct, P);
      case 80:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_compress, splice_atomic>(G, rounds, correct, P);

      case 81:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
      case 82:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
      case 83:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
      case 84:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);

      case 85:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
      case 86:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
      case 87:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
      case 88:
        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);

      case 89:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_compress, split_atomic_one>(G, rounds, correct, P);
      case 90:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_compress, halve_atomic_one>(G, rounds, correct, P);
      case 91:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_compress, splice_simple>(G, rounds, correct, P);
      case 92:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_compress, splice_atomic>(G, rounds, correct, P);

      case 93:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
      case 94:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
      case 95:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
      case 96:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);

      case 97:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
      case 98:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
      case 99:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
      case 100:
        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);

      case 101:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_compress, split_atomic_one>(G, rounds, correct, P);
      case 102:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_compress, halve_atomic_one>(G, rounds, correct, P);
      case 103:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_compress, splice_simple>(G, rounds, correct, P);
      case 104:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_compress, splice_atomic>(G, rounds, correct, P);

      case 105:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
      case 106:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
      case 107:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
      case 108:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);

      case 109:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
      case 110:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
      case 111:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
      case 112:
        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);

      case 113:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_compress, split_atomic_one>(G, rounds, correct, P);
      case 114:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_compress, halve_atomic_one>(G, rounds, correct, P);
      case 115:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_compress, splice_simple>(G, rounds, correct, P);
      case 116:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_compress, splice_atomic>(G, rounds, correct, P);

      case 117:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
      case 118:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
      case 119:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
      case 120:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);

      case 121:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
      case 122:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
      case 123:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
      case 124:
        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);

      /* Liu-Tarjan algorithms */
      case 125:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 126:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
      case 127:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 128:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 129:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 130:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);

      case 131:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 132:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
      case 133:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 134:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 135:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 136:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);

      case 137:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 138:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
      case 139:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 140:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 141:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 142:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);

      case 143:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 144:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
      case 145:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
      case 146:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 147:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
      case 148:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);

      case 149:
        return run_multiple(G, rounds, correct, "bfs_cc", P, t_bfs_cc<Graph>);

      /* Shiloach-Vishkin strategies */
      case 150:
        return run_multiple_sample_only_alg<Graph, sample_kout, shiloachvishkin_cc::SVAlgorithm>(G, rounds, correct, P, "shiloach_vishkin");
      case 151:
        return run_multiple_sample_only_alg<Graph, sample_bfs, shiloachvishkin_cc::SVAlgorithm>(G, rounds, correct, P, "shiloach_vishkin");
      case 152:
        return run_multiple_sample_only_alg<Graph, sample_ldd, shiloachvishkin_cc::SVAlgorithm>(G, rounds, correct, P, "shiloach_vishkin");
      case 153:
        return run_multiple_sample_only_alg<Graph, no_sampling, shiloachvishkin_cc::SVAlgorithm>(G, rounds, correct, P, "shiloach_vishkin");

    case 154:
      return run_multiple_uf_alg<Graph, no_sampling, unite, find_split>(G, rounds, correct, P);

    /* UF Rem-CAS strategies */
    case 155:
      return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
    case 156:
      return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
    case 157:
      return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
    case 158:
      return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);

    case 159:
      return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
    case 160:
      return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
    case 161:
      return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
    case 162:
      return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);

    case 163:
      return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
    case 164:
      return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
    case 165:
      return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
    case 166:
      return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);

    case 167:
      return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
    case 168:
      return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
    case 169:
      return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
    case 170:
      return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);

    default:
      std::cout << "Unknown test" << std::endl;
      return 0.0 ;
    }
  }
} // namespace connectit

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
inline void RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = pbbs::sequence<T>(n + 1, (T)0);
  T cur_comp = 0;
  for (size_t i=0; i<n; i++) {
    T comp = ids[i];
    if (component_map[comp] == 0) {
      component_map[comp] = cur_comp;
      ids[i] = cur_comp;
      cur_comp++;
    } else {
      ids[i] = component_map[comp];
    }
  }
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
      std::cout << "at i = " << i << " cor = " << correct[i] << " got: " << check[i] << std::endl;
      std::cout.flush();
      abort();
    }
    if (correct[i] > max_cor) {
      pbbs::write_max(&max_cor, correct[i], std::less<uintE>());
    }
    if (check[i] > max_chk) {
      pbbs::write_max(&max_chk, check[i], std::less<uintE>());
    }
  }
  cout << "correctness check: " << is_correct << endl;
  cout << "max_cor = " << max_cor << " max_chk = " << max_chk << endl;
}

/* ************************* ***** *************************** */

void print_cpu_stats(cpu_stats& stats, commandLine& P) {
  std::cout <<
    "Stats = { \"ipc\":" + std::to_string(stats.get_ipc())
    + " ,\"total_cycles\":" + std::to_string(stats.get_total_cycles())
    + " ,\"l2_hit_ratio\":" + std::to_string(stats.get_l2_hit_ratio())
    + " ,\"l3_hit_ratio\":" + std::to_string(stats.get_l3_hit_ratio())
    + " ,\"l2_misses\":" + std::to_string(stats.get_l2_misses())
    + " ,\"l2_hits\":" + std::to_string(stats.get_l2_hits())
    + " ,\"l3_misses\":" + std::to_string(stats.get_l3_misses())
    + " ,\"l3_hits\":" + std::to_string(stats.get_l3_hits())
    + " ,\"throughput\":" + std::to_string(stats.get_throughput())
    + "}" << std::endl;
}

template <class Graph>
double Benchmark_runner(Graph& G, commandLine P) {
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);
  bool symmetric = P.getOptionValue("-s");
  int num_tests = 171; // update if new algorithm is added

  int starting_num = P.getOptionIntValue("-start", 0);

  cout << "rounds = " << rounds << endl;
  cout << "num threads = " << num_workers() << endl;

  auto correct = pbbs::sequence<parent>();
  if (P.getOptionValue("-check")) {
    correct = workefficient_cc::CC(G, 0.2, false, true);
    RelabelDet(correct);
  }
  if (test_num == -1) {
    for (int i=starting_num; i < num_tests; i++) {
      std::cout << "running test #: " << i << std::endl;
#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif
      connectit::pick_test(G, i, rounds, P, correct);
#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
  print_cpu_stats(stats, P);
#endif

    }
  } else {
#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif
    connectit::pick_test(G, test_num, rounds, P, correct);
#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
  print_cpu_stats(stats, P);
#endif
  }
  return 1.0;
}

generate_symmetric_once_main(Benchmark_runner, false);
