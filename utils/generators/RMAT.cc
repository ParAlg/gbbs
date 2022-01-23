#include "rmat.h"

#include "gbbs/gbbs.h"

#include <fstream>
#include <iostream>

namespace gbbs {

int BuildRMAT(int argc, char* argv[]) {
  commandLine P(argc, argv, "");

  size_t n = P.getOptionLongValue("-n", 1UL << 20);
  size_t m = P.getOptionLongValue("-m", 1UL << 23);

  double a = P.getOptionDoubleValue("-a", 0.5);
  double b = P.getOptionDoubleValue("-b", 0.1);
  double c = P.getOptionDoubleValue("-c", 0.1);

  auto out_f = P.getOptionValue("-outfile", "");

  if (out_f == "") {
    std::cout << "specify a valid outfile using -outfile" << std::endl;
    abort();
  }

  uintE seed = 4;
  std::cout << "Generating updates" << std::endl;
  auto updates = rmat::generate_updates(n, m, seed, a, b, c);
  std::cout << "Generated updates" << std::endl;

  if (n != (1UL << (parlay::log2_up(n)))) {
    std::cout << "n must be a power of two" << std::endl;
    abort();
  }

  auto C = parlay::sequence_to_string(updates);
  for (size_t i = 0; i < 100; i++) {
    std::cout << std::get<0>(updates[i]) << " " << std::get<1>(updates[i])
              << std::endl;
  }

  size_t nn = C.size();
  std::ofstream file(out_f.c_str(), std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file for writing: " << out_f << std::endl;
    return -1;
  }
  //  file << "# COO Format" << std::endl;
  //  file << "# n = " << n << std::endl;
  //  file << "# m = " << m << std::endl;

  file.write(C.begin(), nn);
  file.close();

  std::cout << "done" << std::endl;
  return 0;
}

}  // namespace gbbs

int main(int argc, char* argv[]) { return gbbs::BuildRMAT(argc, argv); }
