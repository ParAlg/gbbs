// TODO(tomtseng): update this file comment
// Converts a SNAP graph (http://snap.stanford.edu/data/index.html) to
// Ligra adjacency graph format. To symmetrize the graph, pass the "-s"
// flag. For undirected graphs on SNAP, the "-s" flag must be passed
// since each edge appears in only one direction
#include <exception>
#include <vector>

#include "ligra/graph_io.h"
#include "ligra/parse_command_line.h"

namespace {

// Convert `edge_list` to a graph and write it to `output_file`.
template <typename Weight>
void WriteEdgeListAsGraph(
    const char* output_file,
    const std::vector<gbbs_io::Edge<Weight>>& edge_list,
    bool is_symmetric_graph) {
  if (is_symmetric_graph) {
    auto graph{gbbs_io::edge_list_to_symmetric_graph(edge_list)};
    gbbs_io::write_graph_to_file(output_file, graph);
  } else {
    auto graph{gbbs_io::edge_list_to_asymmetric_graph(edge_list)};
    gbbs_io::write_graph_to_file(output_file, graph);
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
