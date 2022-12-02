
#include <math.h>

#include "gbbs/gbbs.h"

//#include "KCore.h"
//#include "gbbs/dynamic_graph_io.h"
//#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

void split(const std::string &s, char delim, std::vector<gbbs::uintE> &elems) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(std::stoi(item));
  }
}

void read_cores(std::string& filename, sequence<std::pair<std::vector<gbbs::uintE>, gbbs::uintE>>& cores) {
  std::ifstream infile(filename);
  std::string line;
  size_t index = 0;
  while (std::getline(infile, line)) {
    std::pair<std::vector<gbbs::uintE>, gbbs::uintE> row_values;
    split(line, '\t', row_values.first);
    row_values.second = row_values.first.back();
    row_values.first.pop_back();
    std::sort(row_values.first.begin(), row_values.first.end());
    cores[index] = row_values;
    index++;
  }
  std::cout << "Size: " << index << std::endl;
  infile.close();
}

// Given approximate cores, output comparisons to exact k-core
void print_stats(std::string& exact_filename, std::string& approx_filename){
  size_t number_of_lines = 0;
  std::string line;
  std::ifstream myfile(exact_filename);
  while (std::getline(myfile, line)) ++number_of_lines;
  myfile.close();

  using PairType = std::pair<std::vector<gbbs::uintE>, gbbs::uintE>;

  sequence<PairType> approx_cores(number_of_lines);
  sequence<PairType> exact_cores(number_of_lines);
  read_cores(exact_filename, exact_cores);
  read_cores(approx_filename, approx_cores);
  //double mult_appx = (2 + 2*eps);

  auto get_core = [&](PairType& p, PairType& q){
    // We want to compare the two vectors in p.first and q.first
    for (size_t i = 0; i < p.first.size(); i++) {
      if (p.first[i] != q.first[i]) return p.first[i] < q.first[i];
    }
    return p.first[0] < q.first[0];
  };
  parlay::sample_sort_inplace(make_slice(exact_cores), get_core);
  parlay::sample_sort_inplace(make_slice(approx_cores), get_core);

    double total_error = 0.0;
    double max_error = 0.0;
    double min_error = std::numeric_limits<double>::max();
    double denominator = 0;
    uintE max_approx_core = 0;
    uintE max_exact_core = 0;
    for (size_t i=0; i<exact_cores.size(); i++) {
      double true_core = exact_cores[i].second;
      double appx_core = approx_cores[i].second;
      if (max_exact_core < true_core) max_exact_core = true_core;
      if (max_approx_core < appx_core) max_approx_core = appx_core;
      if (true_core != 0 && appx_core != 0) {
        auto this_error = std::max(true_core, appx_core) /
          std::min(true_core, appx_core);
        total_error += this_error;
        max_error = std::max(max_error, this_error);
        min_error = std::min(min_error, this_error);
        denominator++;
      }
    }
    min_error = min_error == std::numeric_limits<double>::max() ? 0 : min_error;
    auto avg_error = (denominator == 0) ? 0 : (total_error / denominator);

    std::cout << "### Coreness Estimate (Approx): " << max_approx_core << std::endl;
    std::cout << "### Coreness Exact: " << max_exact_core << std::endl;

    std::cout << "### Per Vertex Average Coreness Error: " << avg_error
      << std::endl;
    std::cout << "### Per Vertex Min Coreness Error: " << min_error << std::endl;
    std::cout << "### Per Vertex Max Coreness Error: " << max_error << std::endl;
}

template <class Graph>
double Compare_runner(Graph& G, commandLine P) {
  std::cout << "### Application: KCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));

  auto exact_str = P.getOptionValue("-exact", "");
  auto approx_str = P.getOptionValue("-approx", "");
  print_stats(exact_str, approx_str);
  return 0;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::Compare_runner, false);