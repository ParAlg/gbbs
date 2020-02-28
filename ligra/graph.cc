#include "graph.h"

std::function<void()> get_deletion_fn(void* a, void* b) {
  auto df = [&](void* aa, void* bb) {
    pbbslib::free_array(aa);
    pbbslib::free_array(bb);
  };
  return std::bind(df, a, b);
}

std::function<void()> get_deletion_fn(void* a, void* b, void* c) {
  auto df = [&](void* aa, void* bb, void* cc) {
    pbbslib::free_array(aa);
    pbbslib::free_array(bb);
    pbbslib::free_array(cc);
  };
  return std::bind(df, a, b, c);
}

std::function<void()> get_deletion_fn(void* a, void* b, void* c, void* d) {
  auto df = [&](void* aa, void* bb, void* cc, void* dd) {
    pbbslib::free_array(aa);
    pbbslib::free_array(bb);
    pbbslib::free_array(cc);
    pbbslib::free_array(dd);
  };
  return std::bind(df, a, b, c, d);
}
