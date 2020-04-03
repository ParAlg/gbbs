#include "ligra/vertex_subset.h"
#include "ligra/vertex.h"
#include "ligra/compressed_vertex.h"
#include "ligra/graph.h"
#include "ligra/graph_io.h"
#include "ligra/ligra.h"

#include "pybind11/pybind11.h"

namespace gbbs_lib {
namespace py = ::pybind11;

/* Defines symmetric vertex functions */
template <template <class W> class vertex_type, class W>
void SymVertexRegister(py::module& m, std::string vertex_name) {
  using vertex = vertex_type<W>;
  /* register vertex */
  py::class_<vertex>(m, vertex_name.c_str());
}

/* Defines symmetric graph functions */
template <template <class W> class vertex_type, class W>
void SymGraphRegister(py::module& m, std::string graph_name) {
  /* register graph */
  using graph = symmetric_graph<vertex_type, W>;
  py::class_<graph>(m, graph_name.c_str());
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

  SymVertexRegister<symmetric_vertex, pbbs::empty>(m, "SymmetricVertexEmpty");
  SymVertexRegister<csv_bytepd_amortized, pbbs::empty>(m, "CompressedSymmetricVertexEmpty");
  SymGraphRegister<symmetric_vertex, pbbs::empty>(m, "SymmetricGraph");
//  SymGraphRegister<csv_bytepd_amortized, pbbs::empty>(m, "SymmetricVertexEmpty", "SymmetricGraph");

  m.def("GetSum", [&] () { return 1; });

  m.def("readSymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_unweighted_symmetric_graph(
        path.c_str(),
        /* mmap = */true);
    alloc_init(G);
    return G;
  });

  m.def("readAsymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_unweighted_asymmetric_graph(
        path.c_str(),
        /* mmap = */true);
    alloc_init(G);
    return G;
  });

  m.def("readCompressedSymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(
        path.c_str(),
        /* mmap = */true,
        /* mmap_copy = */false);
    alloc_init(G);
    return G;
  });

  m.def("readCompressedAsymmetricUnweightedGraph", [&] (std::string& path) {
    auto G = gbbs_io::read_compressed_asymmetric_graph<pbbslib::empty>(
        path.c_str(),
        /* mmap = */true,
        /* mmap_copy = */false);
    alloc_init(G);
    return G;
  });

}

}  // namespace gbbs_lib
