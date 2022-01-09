#pragma once

#ifdef EMPTY_STARTING_GRAPH

namespace gbbs {
bool print_batch_time = false;

template <class Graph, bool provides_initial_graph>
void run_all_tests(Graph& G, size_t n, sequence<incremental_update>& updates,
                   size_t batch_size, size_t insert_to_query, size_t rounds,
                   commandLine P);

/* run synthetic coo */
int RunEmptyStartingGraph(int argc, char* argv[]) {
  auto P = commandLine(argc, argv, "");
  int rounds = P.getOptionIntValue("-r", 5);

  auto in_file = P.getOptionValue("-in_file", "");
  if (in_file == "") {
    std::cout << "must specify valid coo input: " << in_file << std::endl;
    abort();
  }

  sequence<char> S = parlay::chars_from_file(in_file);
  sequence<slice<char>> tokens = parlay::map_tokens(
      parlay::make_slice(S), [](auto x) { return parlay::make_slice(x); });
  // parseback to ints

  assert(tokens.size() % 2 == 0);  // m tuples, into two tokens each

  if (tokens.size() % 2 != 0) {
    std::cout << "malformed coo input" << std::endl;
    std::cout.flush();
    abort();
  }
  size_t m = tokens.size() / 2;
  auto updates = sequence<std::tuple<uintE, uintE, gbbs::empty>>(m);

  uintE n = 0;
  parallel_for(0, m, [&](size_t i) {
    uintE l = parlay::chars_to_int_t<uintE>(tokens[2 * i]);
    uintE r = parlay::chars_to_int_t<uintE>(tokens[2 * i + 1]);
    if (l > n) {
      gbbs::write_min<uintE>(&n, l, std::greater<uintE>());
    }
    if (r > n) {
      gbbs::write_min<uintE>(&n, r, std::greater<uintE>());
    }
    updates[i] = std::make_tuple(l, r, gbbs::empty());
  });
  n = n + 1; /* 0 indexed */
             //  auto sort_f = [&] (const std::tuple<uintE, uintE>& l, const
             //  std::tuple<uintE, uintE>& r) {
             //    return l < r;
             //  };
             //  parlay::sample_sort_inplace(updates.slice(), sort_f);

  size_t batch_size =
      P.getOptionLongValue("-batch_size", 1000000); /* batch size */

  /* fraction of insertions to queries: between [0, 1] */
  double insert_to_query = P.getOptionDoubleValue("-insert_to_query", 0.5);

  bool permute = P.getOptionValue("-permute");
  auto annotated_updates =
      annotate_updates(updates, insert_to_query, n, permute);

  auto FG = edge_array<gbbs::empty>();
  run_all_tests<decltype(FG), false>(FG, n, annotated_updates, batch_size,
                                     insert_to_query, rounds, P);
  return 1;
}

}  // namespace gbbs

int main(int argc, char* argv[]) {
  gbbs::RunEmptyStartingGraph(argc, argv);
  return 1;
}

#endif
