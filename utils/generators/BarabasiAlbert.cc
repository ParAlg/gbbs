#include "barabasi_albert.h"

#include "gbbs/gbbs.h"

#include <fstream>
#include <iostream>

namespace gbbs {

int BuildBarabasiAlbert(int argc, char* argv[]) {
  commandLine P(argc, argv, "");

  size_t n = P.getOptionLongValue("-n", 1UL << 20);
  size_t edges_per_vertex = P.getOptionLongValue("-edges_per_vertex", 10);

  auto out_f = P.getOptionValue("-outfile", "");

  if (out_f == "") {
    std::cout << "specify a valid outfile using -outfile" << std::endl;
    abort();
  }

  auto updates = barabasi_albert::generate_updates(n, edges_per_vertex);

  auto C = parlay::sequence_to_string(updates);

  size_t nn = C.size();
  std::ofstream file(out_f.c_str(), std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file for writing: " << out_f << std::endl;
    return -1;
  }

  file.write(C.begin(), nn);
  file.close();

  std::cout << "done" << std::endl;
  return 0;
}

}  // namespace gbbs

int main(int argc, char* argv[]) {
  return gbbs::BuildBarabasiAlbert(argc, argv);
}
