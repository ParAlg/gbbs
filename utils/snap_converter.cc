// TODO(tomtseng): update this file comment
// Converts a SNAP graph (http://snap.stanford.edu/data/index.html) to
// Ligra adjacency graph format. To symmetrize the graph, pass the "-s"
// flag. For undirected graphs on SNAP, the "-s" flag must be passed
// since each edge appears in only one direction
//
//
// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <exception>

#include "ligra/graph_io.h"
#include "ligra/parse_command_line.h"

namespace {

// Convert `edge_list` to a graph and write it to `output_file`.
template <typename Weight>
void WriteEdgeListAsGraph(
    const char* output_file,
    const pbbs::sequence<gbbs_io::edge<Weight>>& edge_list,
    bool is_symmetric_graph) {
  if (is_symmetric_graph) {
    auto graph{gbbs_io::edge_list_to_symmetric_graph(edge_list)};
    gbbs_io::write_graph_to_file(output_file, &graph);
  } else {
    auto graph{gbbs_io::edge_list_to_asymmetric_graph(edge_list)};
    gbbs_io::write_graph_to_file(output_file, &graph);
  }
}

}

int main(int argc, char* argv[]) {
  // TODO(tomtseng) write a more detailed help string
  const std::string kCommandLineHelpString{
    "[-s] [-w] -i <input SNAP file> -o <output AdjacencyGraph file>"};
  const std::string kInputFlag{"-i"};
  const std::string kOutputFlag{"-o"};

  const commandLine parameters{argc, argv, kCommandLineHelpString};
  const char* const input_file{parameters.getOptionValue(kInputFlag)};
  const char* const output_file{parameters.getOptionValue(kOutputFlag)};
  const bool is_symmetric_graph{parameters.getOption("-s")};
  const bool weighted{parameters.getOption("-w")};

  if (input_file == nullptr || output_file == nullptr) {
    std::cerr << "ERROR: Please specify the input SNAP file with the '" <<
      kInputFlag << "' flag and the output file with the '" << kOutputFlag <<
      "' flag.\n";
    std::terminate();
  }

  if (weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list(input_file)};
    WriteEdgeListAsGraph(output_file, edge_list, is_symmetric_graph);
  } else {
    const auto edge_list{gbbs_io::read_unweighted_edge_list(input_file)};
    WriteEdgeListAsGraph(output_file, edge_list, is_symmetric_graph);
  }

  return 0;
}
