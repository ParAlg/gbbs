#pragma once

#ifdef EMPTY_STARTING_GRAPH

template <class Graph, bool provides_initial_graph>
void run_all_tests(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P);

/* run synthetic coo */
int main(int argc, char* argv[]) {
  auto P = commandLine(argc, argv, "");
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);

  auto in_file = P.getOptionValue("-in_file", "");
  if (in_file == "") {
    std::cout << "must specify valid coo input: " << in_file << std::endl;
    abort();
  }

  sequence<char> chars = pbbs::char_seq_from_file(in_file);
  auto tokens = pbbslib::tokenize(chars, [] (const char c) { return pbbs::is_space(c); });
  // parseback to ints

  assert(tokens.size() % 2 == 0); // m tuples, into two tokens each

  if (tokens.size() % 2 != 0) {
    std::cout << "malformed coo input" << std::endl;
    std::cout.flush();
    abort();
  }
  size_t m = tokens.size() / 2;
  auto updates = pbbs::sequence<std::tuple<uintE, uintE>>(m);


  uintE n = 0;
  parallel_for(0, m, [&] (size_t i) {
    uintE l = std::atoi(tokens[2*i]);
    uintE r = std::atoi(tokens[2*i + 1]);
    if (l > n) {
      pbbs::write_min<uintE>(&n, l, std::greater<uintE>());
    }
    if (r > n) {
      pbbs::write_min<uintE>(&n, r, std::greater<uintE>());
    }
    updates[i] = std::make_tuple(l, r);
  });
  n = n + 1; /* 0 indexed */

  size_t batch_size = P.getOptionLongValue("-batch_size", 1000000); /* batch size */

  /* fraction of insertions to queries: between [0, 1] */
  double insert_to_query = P.getOptionDoubleValue("-insert_to_query", 0.5);
  auto annotated_updates = annotate_updates(updates, insert_to_query);

  auto FG = edge_array<pbbs::empty>();
  run_all_tests<decltype(FG), false>(FG, n, annotated_updates, batch_size, insert_to_query, rounds, P);
}
#endif
