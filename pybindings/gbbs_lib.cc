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
  size_t size = sizeof(size_t) + (n * sizeof(E));
  auto byte_arr = pbbslib::new_array_no_init<uint8_t>(size);
  *(reinterpret_cast<size_t*>(byte_arr)) = n;
  auto arr = reinterpret_cast<E*>(byte_arr + sizeof(size_t));
  parallel_for(0, n, [&] (size_t i) { pbbslib::assign_uninitialized(arr[i], S[i]); });

  py::capsule free_when_done(byte_arr, [](void *f) {
      auto byte_arr = reinterpret_cast<uint8_t*>(f) - sizeof(size_t);
      size_t n = *(byte_arr);
      pbbslib::free_array(byte_arr, n);
      // TODO: should destruct non-trivial objects
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
    }, py::arg("src"))
    .def("StronglyConnectedComponents", [&] (graph& G) {
      auto sccs = StronglyConnectedComponents(G);
      return wrap_array(sccs);
    });
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
  SymGraphRegister<symmetric_vertex, gbbs::empty>(m, "SymmetricGraph");
  SymGraphRegister<csv_bytepd_amortized, gbbs::empty>(m, "CompressedSymmetricGraph");

  AsymVertexRegister<asymmetric_vertex, gbbs::empty>(m, "AsymmetricVertexEmpty");
  AsymVertexRegister<cav_bytepd_amortized, gbbs::empty>(m, "CompressedAsymmetricVertexEmpty");
  AsymGraphRegister<asymmetric_vertex, gbbs::empty>(m, "AsymmetricGraph");
  AsymGraphRegister<cav_bytepd_amortized, gbbs::empty>(m, "CompressedAsymmetricGraph");

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
