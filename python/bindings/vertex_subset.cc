#include "ligra/vertex_subset.h"

#include "pybind11/pybind11.h"

namespace vertex_subset {
namespace py = ::pybind11;

PYBIND11_MODULE(vertex_subset, m) {
  m.doc() = "My extension module";
  m.def("GetSum", [&] () { return 1; });
}

}  // namespace vertex_subset
