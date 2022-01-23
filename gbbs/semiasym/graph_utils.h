#pragma once

namespace gbbs {

struct vertex_data {
  size_t offset;  // offset into edges
  uintE degree;   // vertex degree
};

}  // namespace gbbs
