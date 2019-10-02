#include "../../sequence.h"
#include "../string_basics.h"
#include <string>

using indexT = unsigned int;
using uchar = unsigned char;
indexT* suffixArray(unsigned char *s, size_t n);
indexT* LCP(unsigned char *s, indexT* SA, size_t n);

int main (int argc, char *argv[]) {
  if (argv[1] == NULL) {
    cout << "requires a filename" << endl;
    return 1;
  }

  std::string filename(argv[1]);

  // read file
  pbbs::sequence<char> str = pbbs::char_seq_from_file(filename);
  size_t n = str.size();

  // convert to unsigned char string
  uchar *s = (uchar*) str.begin();
  indexT *SA, *lcp;

  cout << "Testing length = " << n << endl;
  timer t("SA", true);
  for (int i=0; i < 3; i++) {
    t.start();
    SA = suffixArray(s, n);
    t.next("suffix array");
    lcp = LCP(s, SA, n);
    t.next("LCP");
    pbbs::free_array(SA);
    pbbs::free_array(lcp);
  }
}
