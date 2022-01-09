#include <algorithm>
#include <fstream>
#include <ostream>
#include <string>
#include <unordered_map>

#include "gbbs/helpers/parse_command_line.h"

int main(int argc, char* argv[]) {
  const std::string kCommandLineHelpString{
      "Usage: ./remove_singleton_vertices [-wf] -i <input file> -o <output "
      "file>"
      "\n"
      "\n"
      "Removes zero-degree vertices from an edge list by renumbering the\n"
      "vertices.\n"
      "\n"
      "Arguments:\n"
      "  -i <filename>: Path to the input edge-list file.\n"
      "  -o <filename>: Path to the output file.\n"
      "Optional arguments:\n"
      "  --skip <int>: Number of lines in beginning of file to skip\n"
      "  -wf: Use this flag if the edge list is weighted with 32-bit "
      "floats.\n"};
  if (argc < 2 || std::string(argv[1]) == "-h" ||
      std::string(argv[1]) == "--help") {
    std::cout << kCommandLineHelpString << '\n';
    return 0;
  }
  const std::string kInputFlag{"-i"};
  const std::string kOutputFlag{"-o"};
  const gbbs::commandLine parameters{argc, argv, kCommandLineHelpString};
  const char* const input_file{parameters.getOptionValue(kInputFlag)};
  const char* const output_file{parameters.getOptionValue(kOutputFlag)};
  const bool is_float_weighted{parameters.getOption("-wf")};
  const size_t num_skipped_lines{parameters.getOptionLongValue("--skip", 0)};

  std::ifstream input{input_file};
  std::ofstream output{output_file};
  for (size_t i{0}; i < num_skipped_lines; i++) {
    std::string unused;
    std::getline(input, unused);
  }
  uint64_t num_vertices{0};
  uint64_t original_num_vertices{0};
  uint64_t u, v;
  std::unordered_map<uint64_t, uint64_t> vertex_map;
  const auto update_map{[&](const uint64_t vertex_id) {
    if (vertex_map.count(vertex_id) == 0) {
      vertex_map[vertex_id] = num_vertices++;
      if (vertex_id > original_num_vertices) {
        original_num_vertices = vertex_id;
      }
    }
  }};
  uint64_t num_edges{0};

  if (is_float_weighted) {
    float weight;
    while (input >> u >> v >> weight) {
      update_map(u);
      update_map(v);
      output << vertex_map[u] << ' ' << vertex_map[v] << ' ' << weight << '\n';
      num_edges++;
    }
  } else {
    while (input >> u >> v) {
      update_map(u);
      update_map(v);
      output << vertex_map[u] << ' ' << vertex_map[v] << '\n';
      num_edges++;
    }
  }
  original_num_vertices++;

  std::cout << "original number of vertices: " << original_num_vertices << '\n';
  std::cout << "number of vertices         : " << num_vertices << '\n';
  std::cout << "number of edges            : " << num_edges << "\n";
  return 0;
}
