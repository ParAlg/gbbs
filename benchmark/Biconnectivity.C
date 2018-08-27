// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// Usage:
// numactl -i all ./Biconnectivity -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -of : the output file to write the biconnectivity labels
//     -if : the file storing bicc labels, specify if you want to compute the
//           number of biconnected components.
//
// ex computing #biccs:
// > numactl -i all ./Biconnectivity -of orkut_bcs.out -s -m com-orkut.ungraph.txt_SJ
// > numactl -i all ./Biconnectivity -if orkut_bcs.out -s -m com-orkut.ungraph.txt_SJ

#include "Biconnectivity.h"
#include "lib/sparse_additive_map.h"
#include "ligra.h"

template <template <typename W> class vertex, class W>
void BiconnectivityStats(graph<vertex<W>>& GA, char* s, uintE component_id=UINT_E_MAX) {
  size_t n = GA.n;
  _seq<char> S = readStringFromFile(s);
  auto Wo = stringToWords(S.A, S.n);
  auto labels = array_imap<tuple<uintE, uintE>>(n);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    labels[i] =
        make_tuple(atol(Wo.Strings[2 * i]), atol(Wo.Strings[2 * i + 1]));
  });

  auto bits = array_imap<uintE>(n, (uintE)0);
  auto flags = array_imap<bool>(n, false);

  auto bicc_label = [&](const uintE& u, const uintE& v) {
    auto lab_u = labels[u];
    auto lab_v = labels[v];
    uintE p_u = get<0>(lab_u);
    uintE p_v = get<0>(lab_v);
    if (v == p_u) {
      return get<1>(lab_u);
    } else if (u == p_v) {
      return get<1>(lab_v);
    } else {
      assert(get<1>(lab_u) == get<1>(lab_v));
      return get<1>(lab_u);
    }
    exit(0);
    return (uintE)1;
  };


  size_t mask = (1 << 12) - 1;
  auto empty = make_tuple(UINT_E_MAX, UINT_E_MAX);
  auto ST = sparse_additive_map<uintE, uintE>(n, empty);

  auto map_bc_label = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    auto label = bicc_label(src, ngh);
    if (!bits[label]) {
      bits[label] = 1;
    }
    size_t key = (static_cast<size_t>(src) << 32) + static_cast<size_t>(ngh);
    if ((pbbs::hash64(key) & mask) == 0) {
      ST.insert(make_tuple(label, 1));
    }
    if (component_id != UINT_E_MAX) {
      if (label == component_id && !flags[src]) {
        flags[src] = 1;
      }
    }
  };
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), { GA.V[i].mapOutNgh(i, map_bc_label); });

  if (component_id == UINT_E_MAX) {
    auto ET = ST.entries();
    auto cmp_snd = [&] (const tuple<uintE, uintE>& l, const tuple<uintE, uintE>& r) { return get<1>(l) > get<1>(r); };
    pbbs::sample_sort(ET.start(), ET.size(), cmp_snd);
    for (size_t i=0; i<std::min((size_t)10, ET.size()); i++) {
      cout << get<0>(ET[i]) << " " << get<1>(ET[i]) << endl;
    }
  } else {
    // reduce flags
    auto flags_imap = make_in_imap<size_t>(n, [&] (size_t i) { return (size_t)flags[i]; });
    cout << "Largest component size = " << pbbs::reduce_add(flags_imap) << endl;
  }

  // The size of a biconnected component as the number of vertices that
  // participate in it. #edges is another option, but #vertices matches the
  // similar metrics for SCC and CC.
  // Can use sketching to find the top components, and then do an exact
  // component on the top-k using hash-maps, space usage would be O(nk).

  // Note that this is the number of biconnected components excluding isolated
  // vertices (the definition maps edges -> components, so isolated vertices
  // don't contribute to any meaningful components).
  uintE total_biccs = pbbs::scan_add(bits, bits);
  cout << "num biconnected components = " << total_biccs << endl;
}

template <class vertex>
void Biconnectivity_runner(graph<vertex>& GA, commandLine P) {
  auto in_f = P.getOptionValue("-if");
  auto out_f = P.getOptionValue("-of");
  uintE largest_cc_id = static_cast<uintE>(P.getOptionLongValue("-largestid", UINT_E_MAX));
  assert(P.getOptionValue("-s"));
  if (in_f) {
    BiconnectivityStats(GA, in_f, largest_cc_id);
  } else {
    Biconnectivity(GA, out_f);
  }
  // Note that Biconnectivity mutates the graph, so we only run the algorithm
  // once.
  exit(0);
}

generate_main(Biconnectivity_runner, true);
