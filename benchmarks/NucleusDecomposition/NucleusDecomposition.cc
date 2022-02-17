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

#include "NucleusDecomposition.h"
#include "benchmarks/KTruss/KTruss.h"
#include <math.h>
#include <fstream>
namespace gbbs {

long strToTableType(std::string tt_str, std::string inverse_index_map_str) {
  if (tt_str == "ONE_LEVEL") return 1;
  else if (tt_str == "TWO_LEVEL") {
    if (inverse_index_map_str == "BINARY_SEARCH") return 2;
    else if (inverse_index_map_str == "STORED_POINTERS") return 5;
    ABORT("Unexpected inverse index map (str): " << inverse_index_map_str);
  }
  else if (tt_str == "MULTI_LEVEL") {
    if (inverse_index_map_str == "BINARY_SEARCH") return 3;
    else if (inverse_index_map_str == "STORED_POINTERS") return 4;
    ABORT("Unexpected inverse index map (str): " << inverse_index_map_str);
  }
  ABORT("Unexpected number of levels (str): " << tt_str);
}

long strToUpdateAggregation(std::string update_agg_str) {
  if (update_agg_str == "SIMPLE_ARRAY") return 0;
  else if (update_agg_str == "LIST_BUFFER") return 1;
  else if (update_agg_str == "HASH_TABLE") return 2;
  ABORT("Unexpected update aggregation (str): " << update_agg_str);
}

template <class Graph>
double AppNucleusDecomposition_runner(Graph& GA, commandLine P) {
  auto tt_str = P.getOptionValue("--numberOfLevels", "");
  long num_multi_levels = P.getOptionLongValue("--numberOfMultiLevels", 2);
  auto inverse_index_map_str = P.getOptionValue("--inverseIndexMap", "");
  auto update_agg_str = P.getOptionValue("--updateAggregation", "");
  bool relabel = P.getOptionValue("--relabel") || P.getOptionValue("-relabel"); // for true, relabel graph
  bool contiguous_space = P.getOptionValue("--contiguousSpace") || P.getOptionValue("-contig"); // for true, contiguous space
  long r_clique_size = P.getOptionLongValue("--rClique", 3); // k as in k-cliques
  long s_clique_size = P.getOptionLongValue("--sClique", 4); // k as in k-cliques
  bool compact = P.getOptionValue("--compact");

  long r = P.getOptionLongValue("-r", 0); // k as in k-cliques
  long ss = P.getOptionLongValue("-ss", 0); // k as in k-cliques
  if (r == 0) r = r_clique_size;
  if (ss == 0) ss = s_clique_size;

  long table_type = 5;
  long num_levels = 2;
  if (tt_str == "") {
    table_type = P.getOptionLongValue("-tt", 0); // 1 = 1 lvl, 2 = 2 lvls, 3 = multi; 4 = multi nosearch, 5 = 2 lvls nosearch
    num_levels = P.getOptionLongValue("-nl", 0); // only for multi, # levels
  } else {
    table_type = strToTableType(tt_str, inverse_index_map_str);
    num_levels = num_multi_levels;
  }

  long efficient = 1;
  if (update_agg_str == "") {
    efficient = P.getOptionLongValue("-efficient", 1); // for list buffer; 0 = simple array, 1 = list buffer, 2 = hash table, 3 = actually 1 but use pnd, 5 = dynamic list buffer, 4 = more confusing dynamic list buffer that's deprecated
  } else {
    efficient = strToUpdateAggregation(update_agg_str);
  }

  // Internal options
  bool verify = P.getOptionValue("-verify"); // for testing; deprecated
  bool use_compress = P.getOptionValue("-compress"); //only for 2, 3
  bool output_size = P.getOptionValue("-output_size");

  // use_compress only runs compress actually for (2, 3)
  // otherwise, it runs space efficient code, but only if using
  // twotable, twotable_nosearch, or multitable (rest is unimplemented)

  std::cout << "### Application: Nucleus Decomposition" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));

  std::cout << "Internal state: " << std::endl;
  std::cout << "tt: " << table_type << ", nl: " << num_levels << std::endl;
  std::cout << "efficient: " << efficient << ", relabel: " << relabel << ", contig: " << contiguous_space << std::endl;
  std::cout << "End internal state: " << std::endl;

  timer t; t.start();
  if (r == 2 && ss == 3 && table_type == 5 && efficient == 2) {
    KTruss_ht(GA, 16, compact);
  } else {
    NucleusDecomposition(GA, r, ss, table_type, num_levels, relabel, contiguous_space, verify, efficient, use_compress, output_size);
  }

  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  // Require solely one round
  exit(0);

  return tt;
}
}
generate_symmetric_main(AppNucleusDecomposition_runner, false);
