#include "sequence.h"
#include "suffix_tree.h"
#include "string_basics.h"

using namespace pbbs;

using Uchar = unsigned char;
using Uint = uint;
//using Node = node<Uint>;

int main (int argc, char *argv[]) {
  string fname = "/home/guyb/tmp/pbbsbench/testData/sequenceData/data/chr22.dna";
  sequence<char> teststr = char_seq_from_file(fname);
  //string teststr = "aabcabd";
  //string teststr = "bacacdaba";
  size_t n = teststr.size();
  sequence<Uchar> a(n, [&] (size_t i) -> uchar {return teststr[i];});
  //cout << teststr << endl;
  //sequence<Node> R =
  timer t("total Suffix Tree", true);
  for (int i=0; i < 7; i++) {
    t.start();
    suffix_tree<Uint> T(a);
    t.next("");
  }
  //for (size_t i = 0; i < n; i++) cout << R[i].location << ", " << R[i].value << ", " << R[i].parent << endl;
}
