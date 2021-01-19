#include "benchmark.h"
#include "parse_command_line.h"

namespace gbbs {
cpu_stats::cpu_stats() {
  ipc = 0;
  total_cycles = 0;
  l2_hit_ratio = 0;
  l3_hit_ratio = 0;
  l2_misses = 0;
  l2_hits = 0;
  l3_misses = 0;
  l3_hits = 0;
  bytes_read = 0;
  bytes_written = 0;
  total_time = 1.0;
  num_rounds = 1;
}

cpu_stats::cpu_stats(
    double ipc,
    size_t total_cycles,
    double l2_hit_ratio,
    double l3_hit_ratio,
    size_t l2_misses,
    size_t l2_hits,
    size_t l3_misses,
    size_t l3_hits,
    size_t bytes_read,
    size_t bytes_written,
    double total_time,
    size_t num_rounds) :
  ipc(ipc), total_cycles(total_cycles), l2_hit_ratio(l2_hit_ratio),
  l3_hit_ratio(l3_hit_ratio), l2_misses(l2_misses), l2_hits(l2_hits),
  l3_misses(l3_misses), l3_hits(l3_hits), bytes_read(bytes_read),
  bytes_written(bytes_written), total_time(total_time), num_rounds(num_rounds) {}

double cpu_stats::get_ipc() {
  return ipc; /* already an average */
}

size_t cpu_stats::get_total_cycles() {
  return total_cycles / num_rounds;
}

double cpu_stats::get_l2_hit_ratio() {
  return l2_hit_ratio; /* already an average */
}

double cpu_stats::get_l3_hit_ratio() {
  return l3_hit_ratio; /* already an average */
}

size_t cpu_stats::get_l2_misses() {
  return l2_misses / num_rounds;
}

size_t cpu_stats::get_l2_hits() {
  return l2_hits / num_rounds;
}

size_t cpu_stats::get_l3_misses() {
  return l3_misses / num_rounds;
}

size_t cpu_stats::get_l3_hits() {
  return l3_hits / num_rounds;
}

double cpu_stats::get_throughput() {
  constexpr size_t GB = 1024*1024*1024;
  return
    ((static_cast<double>(bytes_read + bytes_written) / total_time) /* bytes/sec */
    / GB); /* GB/sec */
}

#ifdef USE_PCM_LIB

cpu_stats get_pcm_stats(
    SystemCounterState& before_state,
    SystemCounterState& after_state,
    double elapsed,
    size_t rounds) {
  double ipc = getIPC(before_state, after_state);
  size_t total_cycles = getCycles(before_state, after_state);
  double l2_hit_ratio = getL2CacheHitRatio(before_state, after_state);
  double l3_hit_ratio = getL3CacheHitRatio(before_state, after_state);
  size_t l2_misses = getL2CacheMisses(before_state, after_state);
  size_t l2_hits = getL2CacheHits(before_state, after_state);
  size_t l3_misses = getL3CacheMisses(before_state, after_state);
  size_t l3_hits = getL3CacheHits(before_state, after_state);

  size_t bytes_read = getBytesReadFromMC(before_state, after_state);
  size_t bytes_written = getBytesWrittenToMC(before_state, after_state);

  return cpu_stats(
      ipc, total_cycles, l2_hit_ratio, l3_hit_ratio,
      l2_misses, l2_hits, l3_misses, l3_hits,
      bytes_read, bytes_written, elapsed, rounds);
}

void print_pcm_stats(SystemCounterState& before_sstate,
                            SystemCounterState& after_sstate, size_t rounds,
                            double elapsed) {
  std::cout << "# Instructions per clock:        "
            << (getIPC(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# Total Cycles:                  "
            << (getCycles(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# ========= Cache misses/hits ========="
            << "\n";
  std::cout << "# L2 Hit ratio:                  "
            << (getL2CacheHitRatio(before_sstate, after_sstate) / rounds)
            << "\n";
  std::cout << "# L3 Hit ratio:                  "
            << (getL3CacheHitRatio(before_sstate, after_sstate) / rounds)
            << "\n";
  std::cout << "# L2 Misses:                     "
            << (getL2CacheMisses(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# L2 Hits:                       "
            << (getL2CacheHits(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# L3 Misses:                     "
            << (getL3CacheMisses(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# L3 Hits:                       "
            << (getL3CacheHits(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "# ========= Bytes read/written ========="
            << "\n";
  auto bytes_read = getBytesReadFromMC(before_sstate, after_sstate) / rounds;
  auto bytes_written =
      getBytesWrittenToMC(before_sstate, after_sstate) / rounds;
  size_t GB = 1024 * 1024 * 1024;
  auto throughput = ((bytes_read + bytes_written) / elapsed) / GB;
  std::cout << "# Bytes read:                    " << bytes_read << "\n";
  std::cout << "# Bytes written:                 " << bytes_written << "\n";
  std::cout << "# Throughput: " << throughput << " GB/s"
            << "\n";
  std::cout << "# ========= Other statistics ========="
            << "\n";
  std::cout << "# Average relative frequency:    "
            << (getActiveRelativeFrequency(before_sstate, after_sstate) /
                rounds)
            << "\n";
}

void pcm_init() {
  auto* m = PCM::getInstance();
  if (m->program() != PCM::Success) {
    std::cout << "# Could not enable program counters"
              << "\n";
    exit(0);
  }
}

#else

void print_pcm_stats(
    size_t before, size_t after, size_t rounds, double elapsed) {}
void pcm_init() {}

#endif
}  // namespace gbbs
