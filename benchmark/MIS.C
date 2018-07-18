// Usage:
// numactl -i all ./MIS -stats -rounds 4 -s com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc
//     -specfor : run the speculative_for based algorithm from pbbs

#include "MIS.h"
#include "ligra.h"

template <class vertex>
void MIS_runner(graph<vertex>& GA, commandLine P) {
  bool spec_for = P.getOption("-specfor");
  assert(P.getOption("-s"));
  if (spec_for) {
    auto MIS = MIS_spec_for::MIS(GA);
    // in spec_for, MIS[i] == 1 indicates that i was chosen
    auto size_imap = make_in_imap<size_t>(GA.n, [&] (size_t i) { return (MIS[i] == 1); });
    if (P.getOptionValue("-stats")) {
      cout << "MIS size: " << pbbs::reduce_add(size_imap) << endl;
    }
    if (P.getOptionValue("-verify")) {
      verify_MIS(GA, size_imap);
    }
  } else {
    auto MIS = MIS_rootset::MIS(GA);
    auto size_imap = make_in_imap<size_t>(GA.n, [&] (size_t i) { return MIS[i]; });
    if (P.getOptionValue("-stats")) {
      cout << "MIS size: " << pbbs::reduce_add(size_imap) << endl;
    }
    if (P.getOptionValue("-verify")) {
      verify_MIS(GA, size_imap);
    }
  }
}

generate_main(MIS_runner, false);
