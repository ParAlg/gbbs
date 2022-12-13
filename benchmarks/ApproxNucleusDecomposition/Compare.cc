
#include <math.h>
#include <string>
#include <iomanip>
#include <sstream>

#include "gbbs/gbbs.h"

//#include "KCore.h"
//#include "gbbs/dynamic_graph_io.h"
//#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

void split(const std::string& inputString, std::vector<gbbs::uintE>& elems) {
  std::stringstream stringStream(inputString);
  std::string line;
  while(std::getline(stringStream, line)) {
    std::size_t prev = 0, pos;
    while ((pos = line.find_first_of(" ':", prev)) != std::string::npos) {
      if (pos > prev) elems.push_back(stoi(line.substr(prev, pos-prev)));
      prev = pos+1;
    }
    if (prev < line.length()) elems.push_back(stoi(line.substr(prev, std::string::npos)));
  }
}

sequence<std::pair<std::vector<gbbs::uintE>, gbbs::uintE>> read_cores_ktruss(std::string& filename, size_t number_of_lines) {
  sequence<std::pair<std::vector<gbbs::uintE>, gbbs::uintE>> cores(number_of_lines);
  std::ifstream infile(filename);
  std::string line;
  size_t index = 0;
  while (std::getline(infile, line)) {
    std::pair<std::vector<gbbs::uintE>, gbbs::uintE> row_values;
    split(line, row_values.first);
    row_values.second = row_values.first.back();
    row_values.first.pop_back();
    if (row_values[1] > row_values[0]) {
      std::sort(row_values.first.begin(), row_values.first.end());
      cores[index] = row_values;
      index++;
    }
  }
  std::cout << "Size: " << index << std::endl;
  infile.close();
  return cores;
}

sequence<std::pair<std::vector<gbbs::uintE>, gbbs::uintE>> read_cores(std::string& filename, size_t number_of_lines) {
  sequence<std::pair<std::vector<gbbs::uintE>, gbbs::uintE>> cores(number_of_lines);
  std::ifstream infile(filename);
  std::string line;
  size_t index = 0;
  while (std::getline(infile, line)) {
    std::pair<std::vector<gbbs::uintE>, gbbs::uintE> row_values;
    split(line, row_values.first);
    row_values.second = row_values.first.back();
    row_values.first.pop_back();
    std::sort(row_values.first.begin(), row_values.first.end());
    cores[index] = row_values;
    index++;
  }
  std::cout << "Size: " << index << std::endl;
  infile.close();
  return cores;
}

bool check_equal_vectors(std::vector<gbbs::uintE>& a, std::vector<gbbs::uintE>& b) {
  for (size_t i =0; i < a.size(); i++) {
    if (a[i] != b[i]) {
      std::cout << "Mismatch: " << a[i] << ", " << b[i] << std::endl;
      fflush(stdout);
      return false;
    }
  }
  return true;
}

// Given approximate cores, output comparisons to exact k-core
void print_stats(std::string& exact_filename, std::string& approx_filename, bool ktruss=false){
  size_t number_of_lines = 0;
  std::string line;
  std::ifstream myfile(exact_filename);
  while (std::getline(myfile, line)) ++number_of_lines;
  myfile.close();

  std::cout << "Number of lines: " << number_of_lines << std::endl; fflush(stdout);

  using PairType = std::pair<std::vector<gbbs::uintE>, gbbs::uintE>;

  sequence<PairType> exact_cores = read_cores(exact_filename, number_of_lines);
  sequence<PairType> approx_cores;
  if (!ktruss) approx_cores = read_cores(approx_filename, number_of_lines);
  else approx_cores = read_cores_ktruss(approx_filename, number_of_lines);
  //double mult_appx = (2 + 2*eps);

  auto get_core = [&](const PairType& p, const PairType& q) -> bool{
    // We want to compare the two vectors in p.first and q.first
    //std::cout << "P size: " << p.first.size() << ", Q size: " << q.first.size() << std::endl; fflush(stdout);
    for (size_t i = 0; i < p.first.size(); i++) {
      if (p.first[i] != q.first[i]) return p.first[i] < q.first[i];
    }
    return p.first[0] < q.first[0];
  };
  parlay::sample_sort_inplace(make_slice(exact_cores), get_core);
  parlay::sample_sort_inplace(make_slice(approx_cores), get_core);

  std::cout << "Finish sorting" << std::endl; fflush(stdout);

    double total_error = 0.0;
    double max_error = 0.0;
    double min_error = std::numeric_limits<double>::max();
    double denominator = 0;
    uintE max_approx_core = 0;
    uintE max_exact_core = 0;
    for (size_t i=0; i<exact_cores.size(); i++) {
      double true_core = exact_cores[i].second;
      double appx_core = approx_cores[i].second;
      if (!check_equal_vectors(exact_cores[i].first, approx_cores[i].first)) {
        std::cout << "Issue at: " << i << std::endl;
        fflush(stdout);
        exit(0);
      }
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
  bool ktruss = P.getOptionValue("-ktruss");
  print_stats(exact_str, approx_str, ktruss);
  return 0;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::Compare_runner, false);