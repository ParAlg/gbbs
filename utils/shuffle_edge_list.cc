// Shuffle the order of edge-list
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
#include "gbbs/parse_command_line.h"
#include "to_char_arr.h"

namespace gbbs {
namespace {

// Convert `edge_list` to a graph and write it to `output_file`.
template <typename Weight>
void WriteShuffledEdgeList(const char* output_file, const std::vector<gbbs_io::Edge<Weight>>& edge_list) {
    size_t n = edge_list.size();
    auto perm = pbbs::random_permutation<uintE>(n);
    pbbs::sequence<gbbs_io::Edge<Weight>> edge_array = pbbs::sequence<gbbs_io::Edge<Weight>>::no_init(n);
    size_t i = 0;
    for (auto el = edge_list.begin() + i; el != edge_list.begin() + n + 1; ++el) {
        edge_array[i] =  gbbs_io::Edge<Weight>(el->from, el->to, el->weight);
    } 
    // std::vector<gbbs_io::Edge<Weight>>& edge_out(n);
    // i = 0;
    // for (i = 0; i < n; ++i) {
    //     edge_out.emplace_back(edge_array[perm[i]].from, edge_array[perm[i]].to, edge_array[perm[i]].weight);
    // } 
    
    std::ofstream out{output_file};
    // writeArrayToStream(out, edge_array); // m edge triples
    for (i = 0; i < n; ++i) {
       out << edge_array[perm[i]].from << " "<< edge_array[perm[i]].to << " " << edge_array[perm[i]].weight << std::endl;
    } 
    out.close();
    std::cout << "# Wrote file." << std::endl;
}

void WriteShuffledEdgeListUnweighted(const char* output_file, const std::vector<gbbs_io::Edge<pbbs::empty>>& edge_list) {
    size_t n = edge_list.size();
    auto perm = pbbs::random_permutation<uintE>(n);
    pbbs::sequence<gbbs_io::Edge<pbbs::empty>> edge_array = pbbs::sequence<gbbs_io::Edge<pbbs::empty>>::no_init(n);
    size_t i = 0;
    for (auto el = edge_list.begin(); el != edge_list.begin() + n + 1; ++el) {
        edge_array[i] =  gbbs_io::Edge<pbbs::empty>(el->from, el->to);
        i++;
    } 
    // std::vector<gbbs_io::Edge<Weight>>& edge_out(n);
    // i = 0;
    // for (i = 0; i < n; ++i) {
    //     edge_out.emplace_back(edge_array[perm[i]].from, edge_array[perm[i]].to, edge_array[perm[i]].weight);
    // } 
    
    std::ofstream out{output_file};
    // writeArrayToStream(out, edge_array); // m edge triples
    for (i = 0; i < n; ++i) {
       out << edge_array[perm[i]].from << " "<< edge_array[perm[i]].to << std::endl;
    } 
    out.close();
    std::cout << "# Wrote file." << std::endl;
}

int RunShuffleEdgeList(int argc, char* argv[]) {
  const std::string kCommandLineHelpString{
    "Usage: ./shuffle_edge_list [-w] -i <input file> -o <output file>\n"
    "\n"
    "Shuffle the order of edge-list"
    "\n"
    "Arguments:\n"
    "  -i <filename>: Path to the input edge-list file.\n"
    "  -o <filename>: Path to the output adjacency graph file.\n"
    "Optional arguments:\n"
    "  -w: Use this flag if the edge list is weighted with 32-bit integers.\n"
    "  -wf: Use this flag if the edge list is weighted with 32-bit floats.\n"
  };
  const std::string kInputFlag{"-i"};
  const std::string kOutputFlag{"-o"};

  const commandLine parameters{argc, argv, kCommandLineHelpString};
  const char* const input_file{parameters.getOptionValue(kInputFlag)};
  const char* const output_file{parameters.getOptionValue(kOutputFlag)};
  const bool integer_weighted{parameters.getOption("-w")};
  const bool float_weighted{parameters.getOption("-wf")};

  if (argc < 2 ||
      std::string(argv[1]) == "-h" ||
      std::string(argv[1]) == "--help") {
    std::cout << kCommandLineHelpString << '\n';
    return 0;
  }

  if (input_file == nullptr || output_file == nullptr) {
    std::cerr << "ERROR: Please specify the input edge list file with the '" <<
      kInputFlag << "' flag and the output file with the '" << kOutputFlag <<
      "' flag.\n";
    std::terminate();
  }
  if (integer_weighted && float_weighted) {
    std::cerr << "ERROR: Please only specify one weight type.\n";
    std::terminate();
  }

  if (integer_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<int32_t>(input_file)};
    WriteShuffledEdgeList(output_file, edge_list);
  } else if (float_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<float>(input_file)};
    WriteShuffledEdgeList(output_file, edge_list);
  } else {
    const auto edge_list{gbbs_io::read_unweighted_edge_list(input_file)};
    WriteShuffledEdgeListUnweighted(output_file, edge_list);
  }
  return 0;
}

}  // namespace
}  // namespace gbbs

int main(int argc, char* argv[]) {
  return gbbs::RunShuffleEdgeList(argc, argv);
}

