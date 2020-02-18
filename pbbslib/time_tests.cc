//#include <jemalloc/jemalloc.h>
#include "alloc.h"
#include "get_time.h"
#include "ligra/parse_command_line.h"
#include "parallel.h"
#include "sequence_ops.h"
#include "time_operations.h"
#include "utilities.h"

#include <ctype.h>
#include <math.h>
#include <algorithm>
#include <atomic>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

template <typename F>
std::vector<double> repeat(size_t n, size_t rounds, bool check, F test) {
  if (check) test(n, true);
  std::vector<double> R;
  for (size_t i = 0; i < rounds; i++) R.push_back(test(n, false));
  return R;
}

template <typename F>
double reduce(std::vector<double> V, F f) {
  double x = V[0];
  for (size_t i = 1; i < V.size(); i++) x = f(x, V[i]);
  return x;
}

double median(std::vector<double> V) {
  std::sort(V.begin(), V.end());
  if (V.size() % 2 == 1)
    return V[V.size() / 2];
  else
    return (V[V.size() / 2] + V[V.size() / 2 - 1]) / 2.0;
}

double minf(double a, double b) { return (a < b) ? a : b; };
double maxf(double a, double b) { return (a > b) ? a : b; };

bool global_check = false;

template <typename F>
bool run_multiple(size_t n, size_t rounds, float bytes_per_elt,
                  std::string name, F test, bool half_length = 1,
                  std::string x = "bw") {
  std::vector<double> t = repeat(n, rounds, global_check, test);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);
  double rate = n / mint;
  double l = n;
  double tt;
  if (half_length) do {
      l = round(l * .8);
      tt = reduce(repeat(l, rounds, global_check, test), minf);
    } while (tt != 0.0 && l / tt > rate / 2 && l > 1);

  double bandwidth = rate * bytes_per_elt / 1e9;

  cout << name << std::setprecision(3) << ": r=" << rounds << ", med=" << med
       << " (" << mint << "," << maxt << "), "
       << "hlen=" << round(l) << ", " << x << " = " << bandwidth << endl;
  return 1;
}

float bytes_per_read = 1.0;
float bytes_per_write_back = .70;

// the effective number of bytes assuming assymetry
float ebytes(int reads, int write_backs) {
  return reads * bytes_per_read + write_backs * bytes_per_write_back;
}

double pick_test(size_t id, size_t n, size_t rounds, bool half_length) {
  my_mem_pool.clear();
  switch (id) {
    case 0:
      return run_multiple(n, rounds, ebytes(16, 8), "map long", t_map<long>,
                          half_length);
    case 1:
      return run_multiple(n, rounds, ebytes(8, 8), "tabulate long",
                          t_tabulate<long>, half_length);
    case 2:
      return run_multiple(n, rounds, ebytes(8, 0), "reduce add long",
                          t_reduce_add<long>, half_length);
    case 3:
      return run_multiple(n, rounds, ebytes(24, 8), "scan add long",
                          t_scan_add<long>, half_length);
    case 4:
      return run_multiple(n, rounds, ebytes(14, 4), "pack long", t_pack<long>,
                          half_length);
    case 5:
      return run_multiple(n, rounds, ebytes(80, 8), "gather long",
                          t_gather<long>, half_length);
    case 6:
      return run_multiple(n, rounds, ebytes(72, 64), "scatter long",
                          t_scatter<long>, half_length);
    case 7:
      return run_multiple(n, rounds, ebytes(72, 64), "write add long",
                          t_write_add<long>, half_length);
    case 8:
      return run_multiple(n, rounds, ebytes(72, 64), "write min long",
                          t_write_min<long>, half_length);
    case 9:
      return run_multiple(n, rounds, 1, "count sort 8bit long",
                          t_count_sort_8<long>, half_length, "Gelts/sec");
    case 10:
      return run_multiple(n, rounds, 1, "random shuffle long", t_shuffle<long>,
                          half_length, "Gelts/sec");
    case 11:
      return run_multiple(n, rounds, 1, "histogram uint", t_histogram<uint>,
                          half_length, "Gelts/sec");
    case 12:
      return run_multiple(n, rounds, 1, "histogram same uint",
                          t_histogram_same<uint>, half_length, "Gelts/sec");
    case 13:
      return run_multiple(n, rounds, 1, "histogram few uint",
                          t_histogram_few<uint>, half_length, "Gelts/sec");
    case 14:
      return run_multiple(n, rounds, 1, "integer sort<uint,uint>",
                          t_integer_sort_pair<uint>, half_length, "Gelts/sec");
    case 15:
      return run_multiple(n, rounds, 1, "integer sort uint",
                          t_integer_sort<uint>, half_length, "Gelts/sec");
    case 16:
      return run_multiple(n, rounds, 1, "integer sort 128 bits",
                          t_integer_sort_128, half_length, "Gelts/sec");
    case 17:
      return run_multiple(n, rounds, 1, "sort long", t_sort<long>, half_length,
                          "Gelts/sec");
    case 18:
      return run_multiple(n, rounds, 1, "sort uint", t_sort<uint>, half_length,
                          "Gelts/sec");
    case 19:
      return run_multiple(n, rounds, 1, "sort 128 bits", t_sort<__int128>,
                          half_length, "Gelts/sec");
    case 20:
      return run_multiple(n, rounds, ebytes(16, 8), "merge long", t_merge<long>,
                          half_length);
    case 21:
      return run_multiple(n, rounds, ebytes(16 + 5 * 80, 8), "mat vect mult",
                          t_mat_vec_mult<size_t, double>, half_length);
    case 22:
      return run_multiple(n, rounds, ebytes(68, 64), "scatter int",
                          t_scatter<uint>, half_length);
    case 23:
      return run_multiple(n, rounds, 1, "merge sort long", t_merge_sort<long>,
                          half_length, "Gelts/sec");
    case 24:
      return run_multiple(n, rounds, 1, "count sort 2bit long",
                          t_count_sort_2<long>, half_length, "Gelts/sec");
    case 25:
      return run_multiple(n, rounds, ebytes(24, 8), "split3 long",
                          t_split3<long>, half_length);
    case 26:
      return run_multiple(n, rounds, 1, "quicksort long", t_quicksort<long>,
                          half_length, "Gelts/sec");
    case 27:
      return run_multiple(n, rounds, 1, "collect reduce 256 buckets uint",
                          t_collect_reduce_8<uint>, half_length, "Gelts/sec");
    case 28:
      return run_multiple(n, rounds, ebytes(64, 0), "strided read, 128 bytes",
                          t_map_reduce_128, half_length);
    case 29:
      return run_multiple(n, rounds, 1, "collect reduce sparse uint",
                          t_collect_reduce_pair_sparse<uint>, half_length,
                          "Gelts/sec");
    case 30:
      return run_multiple(n, rounds, 1, "remove duplicates",
                          t_remove_duplicates<long>, half_length, "Gelts/sec");
    case 31:
      return run_multiple(n, rounds, 1, "add to bag long", t_bag<long>,
                          half_length, "Gelts/sec");
    case 32:
      return run_multiple(n, rounds, 1, "collect reduce dense uint",
                          t_collect_reduce_pair_dense<uint>, half_length,
                          "Gelts/sec");

    // these are not part of standard suite
    case 50:
      return run_multiple(n, rounds, 1, "histogram reducer",
                          t_histogram_reducer, half_length, "Gelts/sec");
    case 51:
      return run_multiple(n, rounds, ebytes(24, 8), "scan add long seq",
                          t_scan_add_seq<long>, half_length);
    case 52:
      return run_multiple(n, rounds, 1, "range_min long", t_range_min<long>,
                          half_length, "Gelts/sec");
    default:
      assert(false);
      return 0.0;
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  commandLine P(argc, argv,
                "[-n <size>] [-r <rounds>] [-halflen] [-t <testid>]");
  size_t n = P.getOptionLongValue("-n", 100000000);
  int rounds = P.getOptionIntValue("-r", 5);
  int test_num = P.getOptionIntValue("-t", -1);
  bool half_length = P.getOption("-halflen");
  global_check = P.getOption("-check");
  int num_tests = 33;

  cout << "n = " << n << endl;
  cout << "rounds = " << rounds << endl;
  cout << "num threads = " << num_workers() << endl;
  if (half_length)
    cout << "half length on" << endl;
  else
    cout << "half length off" << endl;

  if (test_num == -1)
    for (int i = 0; i < num_tests; i++) pick_test(i, n, rounds, half_length);
  else
    pick_test(test_num, n, rounds, half_length);
  // my_mem_pool.sizes();
}
