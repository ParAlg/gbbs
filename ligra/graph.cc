#include "graph.h"

std::function<void()> get_deletion_fn(void* a, void* b) {
  auto df = [&](void* a, void* b) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
  };
  return std::bind(df, a, b);
}

std::function<void()> get_deletion_fn(void* a, void* b, void* c) {
  auto df = [&](void* a, void* b, void* c) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
    pbbslib::free_array(c);
  };
  return std::bind(df, a, b, c);
}

std::function<void()> get_deletion_fn(void* a, void* b, void* c, void* d) {
  auto df = [&](void* a, void* b, void* c, void* d) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
    pbbslib::free_array(c);
    pbbslib::free_array(d);
  };
  return std::bind(df, a, b, c, d);
}
