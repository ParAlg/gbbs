#include "string_basics.h"

namespace pbbs {

sequence<char> char_seq_from_file(
    const std::string& filename, size_t start, size_t end) {
  std::ifstream file (filename, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << filename << std::endl;
    exit(1);
  }
  size_t length = file.tellg();
  start = std::min(start,length);
  if (end == 0) end = length;
  else end = std::min(end,length);
  size_t n = end - start;
  file.seekg (start, std::ios::beg);
  char* bytes = new_array<char>(n+1);
  file.read (bytes,n);
  file.close();
  return sequence<char>(bytes,n);
}

bool is_space(char c) {
  switch (c)  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

}  // namespace pbbs
