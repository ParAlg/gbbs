#ifndef THIRD_PARTY_GBBS_GBBS_HELPERS_PARALLEL_FOR_WITH_STATUS_H_
#define THIRD_PARTY_GBBS_GBBS_HELPERS_PARALLEL_FOR_WITH_STATUS_H_

#include <cstddef>

#include "absl/status/status.h"
#include "gbbs/bridge.h"
#include "parlay/parallel.h"

namespace gbbs {

// Similar to `parlay::parallel_for`, with the following differences:
// - `f` is expected to return an `absl::Status` instead of `void`.
// - `parallel_for_with_status` returns an `absl::Status` instead of `void`. The
//   returned status is one of the error statuses produced by the calls to `f`,
//   if any, or an Ok status otherwise.
template <typename F>
absl::Status parallel_for_with_status(
    size_t start, size_t end, F&& f,
    // We use `long` for consistency with `parlay::parallel_for`.
    // NOLINTBEGIN(runtime/int)
    // NOLINTBEGIN(google-runtime-int)
    long granularity = 0,
    // NOLINTEND(google-runtime-int)
    // NOLINTEND(runtime/int)
    bool conservative = false) {
  static_assert(std::is_invocable_v<F&, size_t>);
  absl::Status overall_status;
  bool error_found = false;
  parlay::parallel_for(
      start, end,
      [&](size_t i) {
        if (absl::Status status = f(i);
            !status.ok() && gbbs::atomic_compare_and_swap(&error_found, false,
                                                          true)) [[unlikely]] {
          overall_status = status;
        }
      },
      granularity, conservative);
  return overall_status;
}

}  // namespace gbbs

#endif  // THIRD_PARTY_GBBS_GBBS_HELPERS_PARALLEL_FOR_WITH_STATUS_H_
