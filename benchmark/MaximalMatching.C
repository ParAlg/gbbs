// Usage:
// numactl -i all ./MaximalMatching -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "MaximalMatching.h"
#include "ligra.h"

#include "../lib/seq.h"

#include <fstream>
#include <iostream>

template <class vertex>
void MaximalMatching_runner(graph<vertex>& GA, commandLine P) {
  assert(P.getOption("-s")); // input graph must be symmetric

  auto in_f = P.getOptionValue("-if");
  if (in_f) {
    _seq<char> S = readStringFromFile(in_f);
    auto W = stringToWords(S.A, S.n);
    size_t ms = atol(W.Strings[0]);
    using edge = tuple<uintE, uintE>;
    auto matching = sequence<edge>(ms);
    parallel_for(size_t i=0; i<ms; i++) {
      matching[i] = make_tuple(atol(W.Strings[1 + 2*i]), atol(W.Strings[2*(i+1)]));
    }
    verify_matching(GA, matching);
    exit(0);
  }
  auto matching = MaximalMatching(GA);

  // Note that as we mutate the graph by packing out edges, we can't verify that
  // the returned set of edges is a valid maximal matching on the graph
  // currently in-memory. Instead, we write the matching to disk and read it
  // back in to verify correctness.
  auto of = P.getOptionValue("-of");
  if (of) {
    cout << "outfile is = " << of << endl;
    ofstream out(of, ofstream::out);
    out << matching.size() << endl;
    for (size_t i = 0; i < matching.size(); i++) {
      auto e = matching[i];
      out << (get<0>(e) & mm::VAL_MASK) << " " << get<1>(e) << endl;
    }
    out.close();
  }
  // Maximal-matching mutates the underlying graph (unless it is copied, which
  // we don't do to prevent memory issues), so we make sure the algorithm is run
  // exactly once.
  exit(0);
}

generate_main(MaximalMatching_runner, true);
