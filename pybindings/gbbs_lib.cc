#include "gbbs/gbbs.h"
#include "gbbs/compressed_vertex.h"
#include "gbbs/graph.h"
#include "gbbs/graph_io.h"
#include "gbbs/vertex.h"
#include "gbbs/vertex_subset.h"

#include "ApproximateSetCover_lib.h"
#include "BFS_lib.h"
#include "BellmanFord_lib.h"
#include "CC_lib.h"
#include "DeltaStepping_lib.h"
#include "HAC_lib.h"
#include "KCore_lib.h"
#include "MinimumSpanningForest_lib.h"
#include "PageRank_lib.h"

//#include "benchmarks/Biconnectivity/TarjanVishkin/Biconnectivity.h"
//#include "benchmarks/Clustering/SeqHAC/HAC_api.h"
//#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
//#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
//#include "benchmarks/KCore/JulienneDBS17/KCore.h"
//#include "benchmarks/CoSimRank/CoSimRank.h"
//#include "benchmarks/GeneralWeightSSSP/BellmanFord/BellmanFord.h"
//#include "benchmarks/PageRank/PageRank.h"
//#include
//"benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16/StronglyConnectedComponents.h"

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

#include <sys/stat.h>

namespace gbbs {
namespace gbbs_lib {
namespace py = ::pybind11;

// Basic example of handling np arrays
py::array_t<uint32_t> test_array(py::array_t<uint32_t> input) {
  py::buffer_info buf1 = input.request();
  py::array_t<uint32_t> result = py::array_t<uint32_t>(buf1.size);
  py::buffer_info buf2 = result.request();
  uint32_t *ptr1 = (uint32_t *)buf1.ptr, *ptr2 = (uint32_t *)buf2.ptr;
  size_t X = buf1.shape[0];
  size_t Y = buf1.shape[1];
  parallel_for(0, X, [&](size_t i) {
    for (size_t j = 0; j < Y; j++) {
      ptr2[i * Y + j] = ptr1[i * Y + j] + 1;
    }
  });
  result.resize({X, Y});
  return result;
}

template <class W>
auto edgeListToSymmetricWeightedGraph(py::array_t<W> input) {
  py::buffer_info buf = input.request();
  W* ptr = (W*)buf.ptr;

  size_t m = buf.shape[0], cols = buf.shape[1];
  if (cols != 3) {
    std::cerr
        << "Expected [m x 3]-dim array of edges (in CSC format). Quitting."
        << std::endl;
    exit(0);
  }

  using edge = gbbs::gbbs_io::Edge<W>;
  std::vector<edge> edges;
  edges.resize(m);

  parallel_for(0, m, [&](size_t i) {
    //  for (size_t i=0; i<m; i++) {
    uint32_t u = uint32_t(ptr[i * 3]);
    uint32_t v = uint32_t(ptr[i * 3 + 1]);
    edges[i].from = u;
    edges[i].to = v;
    edges[i].weight = ptr[i * 3 + 2];
    //    std::cout << u << " " << v << " " << edges[i].weight << std::endl;
  });

  auto graph{gbbs_io::edge_list_to_symmetric_graph(edges)};

  return graph;
}

/* Defines symmetric vertex functions */
template <template <class W> class vertex_type, class W>
void SymVertexRegister(py::module& m, std::string vertex_name) {
  using vertex = vertex_type<W>;
  /* register vertex */
  py::class_<vertex>(m, vertex_name.c_str()).def("getDegree", [&](vertex& v) {
    return v.out_degree();
  });
}

template <class W>
struct parent_ptr {
  uint32_t parent;
  W wgh;
};

using uint_parent_ptr = parent_ptr<uint32_t>;
using float_parent_ptr = parent_ptr<float>;

template <class W>
auto build_dendrogram(const sequence<std::pair<uintE, W>>& S) {
  size_t n = S.size();
  using par_type = parent_ptr<W>;
  auto arr = (par_type*)(malloc(n * sizeof(par_type)));

  parallel_for(0, n, [&](size_t i) {
    arr[i].parent = S[i].first;
    arr[i].wgh = S[i].second;
  });

  py::capsule free_when_done(arr, [](void* f) { free(f); });

  return py::array_t<par_type>(
      {n},                 // shape
      {sizeof(par_type)},  // C-style contiguous strides
      arr,                 // the data pointer
      free_when_done);     // numpy array references this parent
}

template <class W>
struct edge {
  uint32_t u;
  uint32_t v;
  W wgh;
};

using uint_edge = edge<uint32_t>;
using float_edge = edge<float>;
using double_edge = edge<double>;

template <class W, class EdgeType, class Seq>
auto build_wgh_edgelist(const Seq& S) {
  size_t n = S.size();
  auto arr = (EdgeType*)(malloc(n * sizeof(EdgeType)));

  parallel_for(0, n, [&](size_t i) {
    arr[i].u = std::get<0>(S[i]);
    arr[i].v = std::get<1>(S[i]);
    arr[i].wgh = std::get<2>(S[i]);
  });

  py::capsule free_when_done(arr, [](void* f) { free(f); });

  return py::array_t<EdgeType>(
      {n},                 // shape
      {sizeof(EdgeType)},  // C-style contiguous strides
      arr,                 // the data pointer
      free_when_done);     // numpy array references this parent
}

template <class W, class Seq>
auto build_edgelist(const Seq& S) {
  // Create a Python object that will free the allocated
  // memory when destroyed:
  size_t n = S.size();

  if
    constexpr(std::is_same<W, gbbs::empty>()) {
      py::array_t<uint32_t> result = py::array_t<uint32_t>(2 * n);
      auto buf = (uint32_t*)result.request().ptr;
      parallel_for(0, n, [&](size_t i) {
        buf[2 * i] = S[i].first;
        buf[2 * i + 1] = S[i].second;
      });
      return result;
    }
  else if
    constexpr(std::is_same<W, uint32_t>()) {
      return build_wgh_edgelist<W, uint_edge>(S);
    }
  else if
    constexpr(std::is_same<W, float>()) {
      return build_wgh_edgelist<W, float_edge>(S);
    }
  else {
    static_assert(std::is_same<W, double>());  // otherwise unknown type
    return build_wgh_edgelist<W, double_edge>(S);
  }
}

template <class Seq>
auto wrap_array(const Seq& S) {
  using E = typename Seq::value_type;
  // Create a Python object that will free the allocated
  // memory when destroyed:
  size_t n = S.size();
  auto arr = (E*)(malloc(n * sizeof(E)));
  parallel_for(0, n,
               [&](size_t i) { parlay::assign_uninitialized(arr[i], S[i]); });

  py::capsule free_when_done(arr, [](void* f) { free(f); });

  return py::array_t<E>({n},          // shape
                        {sizeof(E)},  // C-style contiguous strides for double
                        arr,          // the data pointer
                        free_when_done);  // numpy array references this parent
}

/* Defines symmetric graph functions */
template <template <class W> class vertex_type, class W>
void SymGraphRegister(py::module& m, std::string graph_name) {
  /* register graph */
  using graph = symmetric_graph<vertex_type, W>;

  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE,
                                W>::type;

  py::class_<graph>(m, graph_name.c_str())
      .def("numVertices", [](const graph& G) -> size_t { return G.n; })
      .def("numEdges", [](const graph& G) -> size_t { return G.m; })
      .def("writeGraph",
           [](graph& G, std::string& filename) -> void {
             gbbs_io::write_graph_to_file(filename.c_str(), G);
           })
      .def("BFS",
           [&](graph& G, const size_t src) {
             auto parents = compiled::BFS(G, src);
             return wrap_array(parents);
           },
           py::arg("src"))
      .def("Connectivity",
           [&](graph& G) {
             auto ccs = compiled::Connectivity(G);
             return wrap_array(ccs);
           })
      .def("ApproximateSetCover",
           [&](graph& G, size_t num_buckets) {
             auto ccs = compiled::ApproximateSetCover(G, num_buckets);
             return wrap_array(ccs);
           })
      .def("KCore",
           [&](graph& G) {
             auto cores = compiled::KCore(G);
             return wrap_array(cores);
           })
      .def("PageRank",
           [&](graph& G) {
             auto ranks = compiled::PageRank(G);
             return wrap_array(ranks);
           })
      .def("MinimumSpanningForest",
           [&](graph& G) {
             auto G_copy = G;
             auto edges = compiled::MinimumSpanningForest(G_copy);
             return build_edgelist<W>(edges);
           })
      .def("BellmanFord",
           [&](graph& G, uintE source) {
             auto distances = compiled::BellmanFord(G, source);
             return wrap_array(distances);
           })
      .def("DeltaStepping",
           [&](graph& G, uintE source, Distance delta) {
             auto distances = compiled::DeltaStepping(G, source, delta);
             return wrap_array(distances);
           })
      .def("HierarchicalAgglomerativeClustering",
           [&](graph& G, std::string& linkage, bool similarity = true) {
             if
               constexpr(!std::is_same<W, gbbs::empty>()) {
                 auto dendrogram = compiled::HAC(G, linkage, similarity);
                 return build_dendrogram<W>(dendrogram);
               }
             else {
               std::cerr << "Only supported for weighted graphs. (HAC invoked "
                            "with parameters "
                         << linkage << " similarity = " << similarity << ")"
                         << std::endl;
               exit(0);
             }
           });
  //    .def("CoSimRank", [&] (graph& G, const size_t src, const size_t dest) {
  //      CoSimRank(G, src, dest);
  //      return 1.0;
  //    }, py::arg("src"), py::arg("dest"));
}

/* Defines asymmetric vertex functions */
template <template <class W> class vertex_type, class W>
void AsymVertexRegister(py::module& m, std::string vertex_name) {
  using vertex = vertex_type<W>;
  using graph = asymmetric_graph<vertex_type, W>;
  /* register vertex */
  py::class_<vertex>(m, vertex_name.c_str())
      .def("numVertices", [](const graph& G) -> size_t { return G.n; })
      .def("numEdges", [](const graph& G) -> size_t { return G.m; });
}

/* Defines asymmetric graph functions */
template <template <class W> class vertex_type, class W>
void AsymGraphRegister(py::module& m, std::string graph_name) {
  /* register graph */
  using graph = asymmetric_graph<vertex_type, W>;
  py::class_<graph>(m, graph_name.c_str())
      .def("numVertices", [](const graph& G) -> size_t { return G.n; })
      .def("numEdges", [](const graph& G) -> size_t { return G.m; })
      .def("BFS",
           [&](graph& G, const size_t src) {
             auto parents = compiled::BFS(G, src);
             return 1.0;
           },
           py::arg("src"));
}

PYBIND11_MODULE(gbbs_lib, m) {
  m.doc() = "Python module exporting core gbbs types and core data structures.";

  PYBIND11_NUMPY_DTYPE(uint_edge, u, v, wgh);
  PYBIND11_NUMPY_DTYPE(float_edge, u, v, wgh);
  PYBIND11_NUMPY_DTYPE(double_edge, u, v, wgh);

  PYBIND11_NUMPY_DTYPE(uint_parent_ptr, parent, wgh);
  PYBIND11_NUMPY_DTYPE(float_parent_ptr, parent, wgh);

  py::class_<vertexSubset>(m, "VertexSubset")
      .def(py::init<int>(), py::arg("n"))
      .def("size", [](const vertexSubset& vs) -> size_t { return vs.size(); })
      .def("isEmpty",
           [](const vertexSubset& vs) -> bool { return vs.isEmpty(); })
      .def("isDense",
           [](const vertexSubset& vs) -> bool { return vs.dense(); });

  SymVertexRegister<symmetric_vertex, gbbs::empty>(m, "SymmetricVertexEmpty");
  SymVertexRegister<csv_bytepd_amortized, gbbs::empty>(
      m, "CompressedSymmetricVertexEmpty");
  SymGraphRegister<symmetric_vertex, gbbs::empty>(m, "SymmetricGraphEmpty");
  SymGraphRegister<csv_bytepd_amortized, gbbs::empty>(
      m, "CompressedSymmetricGraphEmpty");

  SymVertexRegister<symmetric_vertex, uint32_t>(m, "SymmetricVertexInt");
  SymGraphRegister<symmetric_vertex, uint32_t>(m, "SymmetricGraphInt");

  SymVertexRegister<symmetric_vertex, float>(m, "SymmetricVertexFloat");
  SymGraphRegister<symmetric_vertex, float>(m, "SymmetricGraphFloat");

  AsymVertexRegister<asymmetric_vertex, gbbs::empty>(m,
                                                     "AsymmetricVertexEmpty");
  AsymGraphRegister<asymmetric_vertex, gbbs::empty>(m, "AsymmetricGraphEmpty");

  AsymVertexRegister<asymmetric_vertex, uint32_t>(m, "AsymmetricVertexInt");
  AsymGraphRegister<asymmetric_vertex, uint32_t>(m, "AsymmetricGraphInt");

  AsymVertexRegister<asymmetric_vertex, float>(m, "AsymmetricVertexFloat");
  AsymGraphRegister<asymmetric_vertex, float>(m, "AsymmetricGraphFloat");

  AsymVertexRegister<cav_bytepd_amortized, gbbs::empty>(
      m, "CompressedAsymmetricVertexEmpty");
  AsymGraphRegister<cav_bytepd_amortized, gbbs::empty>(
      m, "CompressedAsymmetricGraphEmpty");

  /* ============================== Graph IO ============================= */
  m.def("readSymmetricUnweightedGraph",
        [&](std::string& path, bool binary = false) {
          auto G = gbbs_io::read_unweighted_symmetric_graph(path.c_str(),
                                                            /* mmap = */ true,
                                                            binary);
          return G;
        });

  m.def("readSymmetricFloatWeightedGraph",
        [&](std::string& path, bool binary = false) {
          auto G = gbbs_io::read_weighted_symmetric_graph<float>(
              path.c_str(),
              /* mmap = */ true, binary);
          return G;
        });

  m.def("readAsymmetricUnweightedGraph",
        [&](std::string& path, bool binary = false) {
          auto G = gbbs_io::read_unweighted_asymmetric_graph(path.c_str(),
                                                             /* mmap = */ true,
                                                             binary);
          return G;
        });

  m.def("readAsymmetricFloatWeightedGraph",
        [&](std::string& path, bool binary = false) {
          auto G = gbbs_io::read_weighted_asymmetric_graph<float>(
              path.c_str(),
              /* mmap = */ true, binary);
          return G;
        });

  m.def("readCompressedSymmetricUnweightedGraph", [&](std::string& path) {
    auto G = gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(
        path.c_str(),
        /* mmap = */ true);
    return G;
  });

  m.def("readCompressedAsymmetricUnweightedGraph", [&](std::string& path) {
    auto G = gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(
        path.c_str(),
        /* mmap = */ true);
    return G;
  });

  m.def("testNumpyArray", &test_array, "testing numpy array");

  m.def("numpyEdgeListToSymmetricUnweightedGraph", [&](py::array_t<uint32_t>
                                                           input) {
    py::buffer_info buf = input.request();
    uint32_t* ptr = (uint32_t*)buf.ptr;

    size_t m = buf.shape[0], cols = buf.shape[1];
    if (cols != 2) {
      std::cerr
          << "Expected [m x 2]-dim array of edges (in CSC format). Quitting."
          << std::endl;
      exit(0);
    }

    using edge = gbbs::gbbs_io::Edge<gbbs::empty>;
    std::vector<edge> edges;
    edges.resize(m);

    for (size_t i = 0; i < m; i++) {
      uint32_t u = ptr[i * 2];
      uint32_t v = ptr[i * 2 + 1];
      edges[i].from = u;
      edges[i].to = v;
    }

    auto graph{gbbs_io::edge_list_to_symmetric_graph(edges)};

    return graph;
  });

  // Uint weighted graph.
  m.def("numpyUintEdgeListToSymmetricWeightedGraph",
        [&](py::array_t<uint32_t> input) {
          std::cout << "Constructing uint weighted graph." << std::endl;
          auto G = edgeListToSymmetricWeightedGraph<uint32_t>(input);
          return G;
        });

  // Uint weighted graph.
  m.def("numpyFloatEdgeListToSymmetricWeightedGraph",
        [&](py::array_t<float> input) {
          std::cout << "Constructing float weighted graph." << std::endl;
          auto G = edgeListToSymmetricWeightedGraph<float>(input);
          return G;
        });

  m.def("loadSymmetricEdgeListAsGraph", [&](std::string& inpath,
                                            std::string& outpath) {
    struct stat buffer;
    bool exists = (stat(outpath.c_str(), &buffer) == 0);
    if (!exists) {
      const auto edge_list{gbbs_io::read_unweighted_edge_list(inpath.c_str())};
      auto graph{gbbs_io::edge_list_to_symmetric_graph(edge_list)};
      gbbs_io::write_graph_to_file(outpath.c_str(), graph);
      return graph;
    }
    auto G = gbbs_io::read_unweighted_symmetric_graph(
        outpath.c_str(),
        /* mmap = */ true,
        /* binary = */ false); /* TODO: use binary */
    return G;
  });

  m.def("loadAsymmetricEdgeListAsGraph", [&](std::string& inpath,
                                             std::string& outpath) {
    struct stat buffer;
    bool exists = (stat(outpath.c_str(), &buffer) == 0);
    if (!exists) {
      const auto edge_list{gbbs_io::read_unweighted_edge_list(inpath.c_str())};
      auto graph{gbbs_io::edge_list_to_asymmetric_graph(edge_list)};
      gbbs_io::write_graph_to_file(outpath.c_str(), graph);
      return graph;
    }
    auto G = gbbs_io::read_unweighted_asymmetric_graph(
        outpath.c_str(),
        /* mmap = */ true,
        /* binary = */ false); /* TODO: use binary */
    return G;
  });

  m.def("loadSymmetricFloatEdgeListAsGraph", [&](std::string& inpath) {
    const auto edge_list{
        gbbs_io::read_weighted_edge_list<float>(inpath.c_str())};
    auto graph{gbbs_io::edge_list_to_symmetric_graph(edge_list)};
    return graph;
  });

  m.def("loadAsymmetricFloatEdgeListAsGraph", [&](std::string& inpath) {
    const auto edge_list{
        gbbs_io::read_weighted_edge_list<float>(inpath.c_str())};
    auto graph{gbbs_io::edge_list_to_asymmetric_graph(edge_list)};
    return graph;
  });
}

}  // namespace gbbs_lib
}  // namespace gbbs
