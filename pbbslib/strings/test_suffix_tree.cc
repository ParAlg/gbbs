#include "sequence.h"
#include "suffix_tree.h"
//#include "suffix_array.h"
//#include "lcp.h"
#include "string_basics.h"
#include <string>
#include <ctype.h>

using namespace pbbs;
using namespace std;

using Uchar = unsigned char;
using Uint = uint;
//using Node = node<Uint>;

int main (int argc, char *argv[]) {

  string fname = "/home/guyb/tmp/pbbsbench/testData/sequenceData/data/wikisamp.xml";
  //string fname = "/ssd0/text/HG18.fasta";
  //string fname = "foo";
  //string fname = "bar";
  sequence<char> teststr = char_seq_from_file(fname);
  //auto r = tokens(teststr, [&] (char c) {return c != 'N';});
  //uint max = reduce(sequence<uint>(r.size(), [&] (size_t i) {return r[i].size();}), maxm<uint>());
  //cout << "max: " << max << endl;
  //sequence<char> out = filter(teststr, [&] (char c) {return (c == 'A') || c == 'C' || c == 'T' || c == 'G' || c == 'M' || c == 'N';});
  //char_seq_to_file(out, "bar");
  //string teststr = "aabcabd";
  //string teststr = "bacacdaba";
  size_t n = teststr.size();
  sequence<Uchar> a(n, [&] (size_t i) -> uchar {return teststr[i];});
  //cout << teststr << endl;
  //sequence<Node> R =
  timer t("total Suffix Tree", true);
  for (int i=0; i < 3; i++) {
    t.start();
    //sequence<Uint> sa = suffix_array<Uint>(a);
    //t.next("sa");
    //sequence<Uint> lcp = LCP(a, sa);
    //t.next("lcp");
    suffix_tree<Uint> T(a);
    t.next("");
    //cout << "wonderful: " << endl;
    //maybe<Uint> r = T.find("wonder");
    //if (r) cout << *r << endl;
  }
  //for (size_t i = 0; i < n; i++) cout << R[i].location << ", " << R[i].value << ", " << R[i].parent << endl;
}
