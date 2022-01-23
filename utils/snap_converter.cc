// Converts a graph in edge-list format to the adjacency graph format that
// GBBS uses. In particular, graph files downloaded from the SNAP dataset
// collection (http://snap.stanford.edu/data/index.html) are often in edge-list
// format and may be converted using this tool.
//
// Pass the `--help` flag to see usage directions.
//
// The contents of an unweighted edge-list file should be as follows:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint>
//     <edge 2 first endpoint> <edge 2 second endpoint>
//     <edge 3 first endpoint> <edge 3 second endpoint>
//     ...
//     <edge m first endpoint> <edge m second endpoint>
//
// The contents of a weighted edge-list file should be as follows:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint> <edge 1 weight>
//     <edge 2 first endpoint> <edge 2 second endpoint> <edge 2 weight>
//     <edge 3 first endpoint> <edge 3 second endpoint> <edge 3 weight>
//     ...
//     <edge m first endpoint> <edge m second endpoint> <edge m weight>
#include <exception>
#include <vector>

#include "gbbs/graph_io.h"
#include "gbbs/helpers/parse_command_line.h"

namespace gbbs {
namespace {

// Convert `edge_list` to a graph and write it to `output_file`.
template <typename Weight>
void WriteEdgeListAsGraph(const char* output_file,
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

int RunSnapConverter(int argc, char* argv[]) {
  const std::string kCommandLineHelpString{
      "Usage: ./snap_converter [-s] [-w] -i <input file> -o <output file>\n"
      "\n"
      "Converts a graph in edge-list format to the adjacency graph format "
      "that\n"
      "GBBS uses. In particular, graph files downloaded from the SNAP dataset\n"
      "collection (http://snap.stanford.edu/data/index.html) are often in\n"
      "edge-list format and may be converted using this tool.\n"
      "\n"
      "Arguments:\n"
      "  -i <filename>: Path to the input edge-list file.\n"
      "  -o <filename>: Path to the output adjacency graph file.\n"
      "Optional arguments:\n"
      "  -s: Treat the edges as a list of undirected edges and create a\n"
      "      symmetric graph. (Without this flag, the edges are treated as a\n"
      "      list of directed edges.)\n"
      "  -w: Use this flag if the edge list is weighted with 32-bit integers.\n"
      "  -wf: Use this flag if the edge list is weighted with 32-bit "
      "floats.\n"};
  const std::string kInputFlag{"-i"};
  const std::string kOutputFlag{"-o"};

  const commandLine parameters{argc, argv, kCommandLineHelpString};
  const char* const input_file{parameters.getOptionValue(kInputFlag)};
  const char* const output_file{parameters.getOptionValue(kOutputFlag)};
  const bool is_symmetric_graph{parameters.getOption("-s")};
  const bool integer_weighted{parameters.getOption("-w")};
  const bool float_weighted{parameters.getOption("-wf")};

  if (argc < 2 || std::string(argv[1]) == "-h" ||
      std::string(argv[1]) == "--help") {
    std::cout << kCommandLineHelpString << '\n';
    return 0;
  }

  if (input_file == nullptr || output_file == nullptr) {
    std::cerr << "ERROR: Please specify the input SNAP file with the '"
              << kInputFlag << "' flag and the output file with the '"
              << kOutputFlag << "' flag.\n";
    std::terminate();
  }
  if (integer_weighted && float_weighted) {
    std::cerr << "ERROR: Please only specify one weight type.\n";
    std::terminate();
  }

  if (integer_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<int32_t>(input_file)};
    WriteEdgeListAsGraph(output_file, edge_list, is_symmetric_graph);
  } else if (float_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<float>(input_file)};
    WriteEdgeListAsGraph(output_file, edge_list, is_symmetric_graph);
  } else {
    const auto edge_list{gbbs_io::read_unweighted_edge_list(input_file)};
    WriteEdgeListAsGraph(output_file, edge_list, is_symmetric_graph);
  }
  return 0;
}

}  // namespace
}  // namespace gbbs

int main(int argc, char* argv[]) { return gbbs::RunSnapConverter(argc, argv); }
