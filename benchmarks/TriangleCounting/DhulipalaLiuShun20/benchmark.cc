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

#include "shared.h"
#include "dynamic_graph.h"
#include "preprocess.h"
#include "rebalancing.h"

namespace gbbs {
namespace DBTInternal {


  template <class EdgeT>
  pbbs::sequence<pair<EdgeT, bool>> generateEdgeUpdates(size_t nv, size_t numEdges, const commandLine& P){
    return pbbs::sequence<pair<EdgeT, bool>>(numEdges, [&](size_t i) { 
      uintE u,v = 0;
      do{u = rand()%nv;v = rand()%nv;}while (u == v);

      return make_pair(EdgeT(u,v), true);
     });

  }

  ///// PRINTING UTILITIES //////
  inline void PrintCaption(std::string t_in) {
    std::cout << "========= " << t_in << " =========" << std::endl;
  }

  inline void PrintSubcaption(std::string t_in) {
    std::cout << "========= " << t_in << std::endl;
  }

  template <class T>
  inline void PrintFunctionTitle(std::string t_func, T t_suffix) {
    std::cout << "[" << t_func << "] " << t_suffix << std::endl;
  }

  template <class T>
  inline void PrintFunctionItem(std::string t_func, std::string t_item, T t_suffix) {
    std::cout << "[" << t_func << "] " << t_item << " = " << t_suffix << std::endl;
  }

  template <class T>
  inline void PrintSubtimer(std::string t_item, T t_suffix) {
    std::cout << "::" << t_item << ": " << t_suffix << std::endl;
  }

  inline void PrintBreak() {
    std::cout << std::endl;
  }

  template <typename A>
  inline void PrintVec(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	      std::cout << t_vec[i] << ' ';
      }
      std::cout << std::endl;
    }

  template <typename A>
  inline void PrintVec(A t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	      std::cout << t_vec[i] << ' ';
      }
      std::cout << std::endl;
    }

  template <typename A>
  inline void PrintPairSeq(A t_vec) {
      for (size_t i = 0; i < t_vec.size(); ++ i) {
	      std::cout << "(" << t_vec[i].first << "," << t_vec[i].second << ") ";
      }
      std::cout << std::endl;
  }

  template <typename A>
  inline void PrintVec2(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	     t_vec[i].print();
      }
      std::cout << std::endl;
    }


//num edge is multiple of 10, already randomly shuffled
template <class UT>
inline void staticCount(const std::vector<UT>& edges, int num_batch, commandLine& P) {
//   size_t block_size = 0; 
  size_t batch_size = edges.size()/10;
  PrintFunctionItem("Static", "batchsize", batch_size);
  DBTGraph::DyGraph<DBTGraph::SymGraph> DGnew;
  for(int i = 0; i< num_batch; ++i){
    size_t batch_end = min(batch_size  * (i+1), edges.size());
    timer t; t.start();
    DBTGraph::majorRebalancing(edges, 0,  batch_end, 0, DGnew, P, false);
    PrintFunctionItem("Static", "batch", i);
    t.stop();t.reportTotal("");
    PrintBreak();

  }
}


// example: ./benchmark -i ../../../inputs/rMatEdgesShuffle
int RunBenchmark(int argc, char* argv[]) {
  const std::string kCommandLineHelpString{
    "Usage: ./benchmark [-static] [-w] [-n 10] -i <input file>\n"
    "\n"
    "Benchmark Triangle Counting"
    "\n"
    "Arguments:\n"
    "  -i <filename>: Path to the input edge-list file.\n"
    "Optional arguments:\n"
    "  -static: run static algorithm (ShunTangwongsan15).\n"
    "  -n: number of batches. default to 10.\n"
    "  -w: Use this flag if the edge list is weighted with 32-bit integers.\n"
    "  -wf: Use this flag if the edge list is weighted with 32-bit floats.\n"
  };

    const std::string kInputFlag{"-i"};

  commandLine parameters{argc, argv, kCommandLineHelpString};
  const char* const input_file{parameters.getOptionValue(kInputFlag)};
  const bool integer_weighted{parameters.getOption("-w")};
  const bool float_weighted{parameters.getOption("-wf")};
  const int num_batch{parameters.getOptionIntValue("-n", 10)};

//   if (argc < 2 ||
//       std::string(argv[1]) == "-h" ||
//       std::string(argv[1]) == "--help") {
//     std::cout << kCommandLineHelpString << '\n';
//     return 0;
//   }

//   if (input_file == nullptr || output_file == nullptr) {
//     std::cerr << "ERROR: Please specify the input edge list file with the '" <<
//       kInputFlag << "' flag and the output file with the '" << kOutputFlag <<
//       "' flag.\n";
//     std::terminate();
//   }
  if (integer_weighted && float_weighted) {
    std::cerr << "ERROR: Please only specify one weight type.\n";
    std::terminate();
  }

  if (integer_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<int32_t>(input_file)};
    staticCount(edge_list, num_batch, parameters);
  } else if (float_weighted) {
    const auto edge_list{gbbs_io::read_weighted_edge_list<float>(input_file)};
    staticCount(edge_list, num_batch, parameters);
  } else {
    const auto edge_list{gbbs_io::read_unweighted_edge_list(input_file)};
    staticCount(edge_list, num_batch, parameters);
  }
  return 0;
} 

}  // namespace DBTInternal
}  // namespace gbbs

int main(int argc, char* argv[]) {
  return gbbs::DBTInternal::RunBenchmark(argc, argv);
}

