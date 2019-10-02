#define NOTMAIN 1
#include "../suffix_array.h"
#include "../lcp.h"

using indexT = unsigned int;
using uchar = unsigned char;

indexT* suffixArray(unsigned char *s, size_t n) {
  return pbbs::suffix_array<indexT>(pbbs::make_range(s, s+n)).to_array();
}

indexT* LCP(unsigned char *s, indexT* SA, size_t n) {
  return pbbs::lcp(pbbs::make_range(s, s+n), pbbs::make_range(SA, SA+n)).to_array();
}


