#include "gbbs/gbbs.h"
#include "gbbs/vertex_subset.h"
#include "gbbs/vertex.h"
#include "gbbs/compressed_vertex.h"
#include "gbbs/graph.h"
#include "gbbs/graph_io.h"

#include "benchmarks/BFS/NonDeterministicBFS/BFS.h"
#include "benchmarks/Biconnectivity/TarjanVishkin/Biconnectivity.h"
#include "benchmarks/CliqueCounting/Clique.h"
#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/CoSimRank/CoSimRank.h"
#include "benchmarks/PageRank/PageRank.h"
#include "benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16/StronglyConnectedComponents.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include <sys/stat.h>

namespace gbbs {
namespace gbbs_lib {
namespace py = ::pybind11;

// Basic example of handling np arrays
py::array_t<uint32_t> test_array(py::array_t<uint32_t> input) {
  py::buffer_info buf1 = input.request();
  py::array_t<uint32_t> result = py::array_t<uint32_t>(buf1.size);
  py::buffer_info buf2 = result.request();
  uint32_t *ptr1 = (uint32_t *) buf1.ptr,
           *ptr2 = (uint32_t *) buf2.ptr;
  int X = buf1.shape[0];
  int Y = buf1.shape[1];
  parallel_for(0, X, [&] (size_t i) {
    for (size_t j=0; j<Y; j++) {
      ptr2[i*Y + j] = ptr1[i*Y + j] + 1;
    }
  });
  result.resize({X,Y});
  return result;
}

template <class W>
auto edgeListToSymmetricWeightedGraph(py::array_t<W> input) {
  py::buffer_info buf = input.request();
  W* ptr = (W*)buf.ptr;

  size_t m = buf.shape[0], cols = buf.shape[1];
  if (cols != 3) {
    std::cout << "Expected [m x 3]-dim array of edges (in CSC format). Quitting." << std::endl;
    exit(0);
  }

  using edge = gbbs::gbbs_io::Edge<W>;
  std::vector<edge> edges;
  edges.resize(m);

  for (size_t i=0; i<m; i++) {
    uint32_t u = uint32_t(ptr[i*3]);
    uint32_t v = uint32_t(ptr[i*3 + 1]);
    edges[i].from = u;
    edges[i].to = v;
    edges[i].weight = ptr[i*3 + 2];
  }

  auto graph{gbbs_io::edge_list_to_symmetric_graph(edges)};

  return graph;
}


/* Defines symmetric vertex functions */
template <template <class W> class vertex_type, class W>
void SymVertexRegister(py::module& m, std::string vertex_name) {
  using vertex = vertex_type<W>;
  /* register vertex */
  py::class_<vertex>(m, vertex_name.c_str())
    .def("getDegree", [&] (vertex& v) {
      return v.out_degree();
    });
}

template <class Seq>
auto wrap_array(const Seq& S) {
  using E = typename Seq::value_type;
  // Create a Python object that will free the allocated
  // memory when destroyed:
  size_t n = S.size();
  auto arr = (E*)(malloc(n*sizeof(E)));
  parallel_for(0, n, [&] (size_t i) { pbbslib::assign_uninitialized(arr[i], S[i]); });

  py::capsule free_when_done(arr, [](void *f) {
      free(f);
  });

  return py::array_t<E>(
      {n}, // shape
      {sizeof(E)}, // C-style contiguous strides for double
      arr, // the data pointer
      free_when_done); // numpy array references this parent
}

/* Defines symmetric graph functions */
template <template <class W> class vertex_type, class W>
void SymGraphRegister(py::module& m, std::string graph_name) {
  /* register graph */
  using graph = symmetric_graph<vertex_type, W>;
  py::class_<graph>(m, graph_name.c_str())
    .def("numVertices", [](const graph& G) -> size_t {
      return G.n;
    })
    .def("numEdges", [](const graph& G) -> size_t {
      return G.m;
    })
    .def("BFS", [&] (graph& G, const size_t src) {
      auto parents = BFS(G, src);
      return 1.0;
    }, py::arg("src"))
    .def("Connectivity", [&] (graph& G) {
      auto ccs = workefficient_cc::CC(G);
      return wrap_array(ccs);
    })
    .def("KCore", [&] (graph& G) {
      auto cores = KCore(G);
      return wrap_array(cores);
    })
    .def("PageRank", [&] (graph& G) {
      auto ranks = PageRank(G);
      return wrap_array(ranks);
    })
    .def("CoSimRank", [&] (graph& G, const size_t src, const size_t dest) {
      CoSimRank(G, src, dest);
      return 1.0;
    }, py::arg("src"), py::arg("dest"));
}

/* Defines asymmetric vertex functions */
template <template <class W> class vertex_type, class W>
void AsymVertexRegister(py::module& m, std::string vertex_name) {
  using vertex = vertex_type<W>;
  using graph = asymmetric_graph<vertex_type, W>;
  /* register vertex */
  py::class_<vertex>(m, vertex_name.c_str())
    .def("numVertices", [](const graph& G) -> size_t {
      return G.n;
    })
    .def("numEdges", [](const graph& G) -> size_t {
      return G.m;
    })
    .def("BFS", [&] (graph& G, const size_t src) {
      auto parents = BFS(G, src);
      return 1.0;
    });
}

/* Defines asymmetric graph functions */
template <template <class W> class vertex_type, class W>
void AsymGraphRegister(py::module& m, std::string graph_name) {
  /* register graph */
  using graph = asymmetric_graph<vertex_type, W>;
  py::class_<graph>(m, graph_name.c_str())
    .def("numVertices", [](const graph& G) -> size_t {
      return G.n;
    })
    .def("numEdges", [](const graph& G) -> size_t {
      return G.m;
    })
    .def("BFS", [&] (graph& G, const size_t src) {
      auto parents = BFS(G, src);
      return 1.0;
    }, py::arg("src"));
}

PYBIND11_MODULE(gbbs_lib, m) {
  m.doc() = "Python module exporting core gbbs types and core data structures.";

  py::class_<vertexSubset>(m, "VertexSubset")
    .def(py::init<int>(), py::arg("n"))
    .def("size", [](const vertexSubset& vs) -> size_t {
      return vs.size();
    })
    .def("isEmpty", [](const vertexSubset& vs) -> bool {
      return vs.isEmpty();
    })
    .def("isDense", [](const vertexSubset& vs) -> bool {
      return vs.dense();
    });

  SymVertexRegister<symmetric_vertex, gbbs::empty>(m, "SymmetricVertexEmpty");
  SymVertexRegister<csv_bytepd_amortized, gbbs::empty>(m, "CompressedSymmetricVertexEmpty");
  SymGraphRegister<symmetric_vertex, gbbs::empty>(m, "SymmetricGraphEmpty");
  SymGraphRegister<csv_bytepd_amortized, gbbs::empty>(m, "CompressedSymmetricGraphEmpty");

  SymVertexRegister<symmetric_vertex, uint32_t>(m, "SymmetricVertexInt");
  SymGraphRegister<symmetric_vertex, uint32_t>(m, "SymmetricGraphInt");

  AsymVertexRegister<asymmetric_vertex, gbbs::empty>(m, "AsymmetricVertexEmpty");
  AsymVertexRegister<cav_bytepd_amortized, gbbs::empty>(m, "CompressedAsymmetricVertexEmpty");
  AsymGraphRegister<asymmetric_vertex, gbbs::empty>(m, "AsymmetricGraphEmpty");
  AsymGraphRegister<cav_bytepd_amortized, gbbs::empty>(m, "CompressedAsymmetricGraphEmpty");

  /* ============================== Graph IO ============================= */
  m.def("readSymmetricUnweightedGraph", [&] (std::string& path, bool binary=false) {
    auto G = gbbs_io::read_unweighted_symmetric_graph(
        path.c_str(),
        /* mmap = */true,
        binary);
    alloc_init(G);
    return G;
  });

  m.def("readAsymmetricUnweightedGraph", [&] (std::string& path, bool binary=false) {
    auto G = gbbs_io::read_unweighted_asymmetric_graph(
        path.c_str(),
        /* mmap = */true,
        binary);
    alloc_init(G);
    return G;
  });

  m.def("readCompressedSymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(
        path.c_str(),
        /* mmap = */true,
        /* mmap_copy = */false);
    alloc_init(G);
    return G;
  });

  m.def("readCompressedAsymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(
        path.c_str(),
        /* mmap = */true,
        /* mmap_copy = */false);
    alloc_init(G);
    return G;
  });

  m.def("testNumpyArray", &test_array, "testing numpy array");

  m.def("numpyEdgeListToSymmetricUnweightedGraph", [&] (py::array_t<uint32_t> input) {
    py::buffer_info buf = input.request();
    uint32_t* ptr = (uint32_t*)buf.ptr;

    size_t m = buf.shape[0], cols = buf.shape[1];
    if (cols != 2) {
      std::cout << "Expected [m x 2]-dim array of edges (in CSC format). Quitting." << std::endl;
      exit(0);
    }

    using edge = gbbs::gbbs_io::Edge<gbbs::empty>;
    std::vector<edge> edges;
    edges.resize(m);

    for (size_t i=0; i<m; i++) {
      uint32_t u = ptr[i*2];
      uint32_t v = ptr[i*2 + 1];
      edges[i].from = u;
      edges[i].to = v;
    }

    auto graph{gbbs_io::edge_list_to_symmetric_graph(edges)};

    return graph;
  });

  // Integer weighted graphs.
  m.def("numpyEdgeListToSymmetricWeightedGraph", [&] (py::array_t<uint32_t> input) {
      return edgeListToSymmetricWeightedGraph(input);
  });


  m.def("loadSymmetricEdgeListAsGraph", [&] (std::string& inpath, std::string& outpath) {
    struct stat buffer;
    bool exists = (stat (outpath.c_str(), &buffer) == 0);
    if (!exists) {
      const auto edge_list{gbbs_io::read_unweighted_edge_list(inpath.c_str())};
      auto graph{gbbs_io::edge_list_to_symmetric_graph(edge_list)};
      gbbs_io::write_graph_to_file(outpath.c_str(), graph);
      return graph;
    }
    auto G = gbbs_io::read_unweighted_symmetric_graph(
        outpath.c_str(),
        /* mmap = */true,
        /* binary = */false);  /* TODO: use binary */
    alloc_init(G);
    return G;
  });

  m.def("loadAsymmetricEdgeListAsGraph", [&] (std::string& inpath, std::string& outpath) {
    struct stat buffer;
    bool exists = (stat (outpath.c_str(), &buffer) == 0);
    if (!exists) {
      const auto edge_list{gbbs_io::read_unweighted_edge_list(inpath.c_str())};
      auto graph{gbbs_io::edge_list_to_asymmetric_graph(edge_list)};
      gbbs_io::write_graph_to_file(outpath.c_str(), graph);
      return graph;
    }
    auto G = gbbs_io::read_unweighted_asymmetric_graph(
        outpath.c_str(),
        /* mmap = */true,
        /* binary = */ false);  /* TODO: use binary */
    alloc_init(G);
    return G;
  });

}

}  // namespace gbbs_lib
}  // namespace gbbs
