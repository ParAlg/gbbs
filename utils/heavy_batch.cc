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
#include "gbbs/macros.h"

namespace gbbs {
namespace {


void SampleHighDegEdges(const char* output_file, const std::vector<gbbs_io::Edge<pbbs::empty>>& edge_list, int frac, size_t numEdges) {
 using weight_type = pbbs::empty;

  pbbs::sequence<gbbs_io::Edge<weight_type>> edges_both_directions(2 * edge_list.size());
  par_for(0, edge_list.size(), [&](const size_t i) {
      const gbbs_io::Edge<weight_type>& edge = edge_list[i];
      edges_both_directions[2 * i] = edge;
      edges_both_directions[2 * i + 1] =
        gbbs_io::Edge<weight_type>{edge.to, edge.from, edge.weight};
  });
  constexpr auto compare_endpoints = [](
      const gbbs_io::Edge<weight_type>& left,
      const gbbs_io::Edge<weight_type>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  pbbs::sequence<gbbs_io::Edge<weight_type>> edges =
    pbbs::remove_duplicates_ordered(edges_both_directions, compare_endpoints);
  // const size_t num_edges = edges.size();
  const size_t num_vertices = gbbs_io::internal::get_num_vertices_from_edges(edges);
  pbbs::sequence<vertex_data> vertex_data =
    gbbs_io::internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

    pbbs::sequence<size_t> inds = pbbs::sequence<size_t>(vertex_data.size(), [&](size_t i){return i;});


    pbbs::sequence<size_t> sorted_inds = pbbs::sample_sort(inds, 
    [&](const size_t &i, const  size_t &j){
        return vertex_data[i].degree > vertex_data[j].degree;
    });

    std::cout <<  "sorted" << std::endl;

    size_t top = sorted_inds.size() / frac;

    std::ofstream out{output_file};
    // writeArrayToStream(out, edge_array); // m edge triples
    for (size_t i = 0; i < numEdges; ++i) {
        uintE u,v = 0;
      do{u = rand()%top;v = rand()%top;}while (u == v);

       out << inds[u] << " "<< inds[v] << std::endl;
    } 
    out.close();
    std::cout << "# Wrote file." << std::endl;
}

int RunSampleHighDegEdges(int argc, char* argv[]) {
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
  // const bool integer_weighted{parameters.getOption("-w")};
  // const bool float_weighted{parameters.getOption("-wf")};
  const int frac{parameters.getOptionIntValue("-frac", 100)};
  const size_t numEdges{parameters.getOptionLongValue("-ne", 10000)};


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
  // if (integer_weighted && float_weighted) {
  //   std::cerr << "ERROR: Please only specify one weight type.\n";
  //   std::terminate();
  // }

//   if (integer_weighted) {
//     const auto edge_list{gbbs_io::read_weighted_edge_list<int32_t>(input_file)};
//     SampleHighDegEdges(output_file, edge_list, frac, numEdges);
//   } else if (float_weighted) {
//     const auto edge_list{gbbs_io::read_weighted_edge_list<float>(input_file)};
//     SampleHighDegEdges(output_file, edge_list, frac, numEdges);
//   } else {
    const auto edge_list{gbbs_io::read_unweighted_edge_list(input_file)};
    SampleHighDegEdges(output_file, edge_list, frac, numEdges);
//   }
  return 0;
}

}  // namespace
}  // namespace gbbs

int main(int argc, char* argv[]) {
  return gbbs::RunSampleHighDegEdges(argc, argv);
}

