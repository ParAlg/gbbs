#pragma once

// This library contains utilities for reporting the progress of the execution
// of an algorithm.

#include <cstdint>

#include "absl/base/attributes.h"
#include "absl/functional/any_invocable.h"
#include "absl/functional/function_ref.h"
#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"

namespace gbbs {

// Callable object to be used for reporting the progress of the execution of an
// algorithm, as a number in the closed interval [0.0, 1.0]. Algorithms must
// report non-decreasing values in this range; values closer to 1.0 indicate
// that the algorithm is closer to completion.
//
// Note: beyond the non-decreasing requirement, there's no requirement on the
// concrete meaning of the reported numbers (in particular, the progress numbers
// reported need not be proportional to time elapsed since the start of the
// algorithm). There's also no requirement on how many times, or how often, the
// object is called.
//
// Typically, an algorithm will make at least two calls; the first one reporting
// a value of 0.0 (at the very beginning of the execution of the algorithm), and
// the final one reporting a value of 1.0 (at the very end of the execution of
// the algorithm).
using ReportProgressCallback = absl::AnyInvocable<void(float progress)>;

// Wrapper for a `ReportProgress` instance that's useful for algorithms with the
// following structure:
//
//   - Preprocessing (optional)
//   - One or more iterations
//   - Postprocessing (optional)
//
// This class provides a convenient way to report the progress of the algorithm
// execution. It makes the simplifying assumption that the preprocessing step
// (if any), postprocessing step (if any), and each iteration all represent
// equal amounts of progress.
class IterationProgressReporter {
 public:
  // Creates a new `IterationProgressReporter` instance.
  //
  // Arguments:
  // * `report_progress`: the instance to be used for reporting the progress of
  //   the algorithm execution.
  // * `num_iterations`: the number of iterations (must be strictly positive).
  // * `has_preprocessing`: indicates whether the algorithm has a preprocessing
  //   step.
  // * `has_postprocessing`: indicates whether the algorithm has a
  //   postprocessing step.
  IterationProgressReporter(absl::FunctionRef<void(float)> report_progress
                                ABSL_ATTRIBUTE_LIFETIME_BOUND,
                            int64_t num_iterations,
                            bool has_preprocessing = false,
                            bool has_postprocessing = false)
      : report_progress_(report_progress),
        num_iterations_(num_iterations),
        has_preprocessing_(has_preprocessing),
        has_postprocessing_(has_postprocessing),
        total_steps_(num_iterations + (has_preprocessing ? 1 : 0) +
                     (has_postprocessing ? 1 : 0)) {}

  // Indicate that the preprocessing step is complete. Returns an error status
  // if `has_preprocessing_` is false or if any call of `IterationComplete()`
  // has been made.
  absl::Status PreprocessingComplete();

  // Indicates that the `iteration`-th iteration is complete. The `iteration`
  // value must be in the half-open interval [0, `num_iterations_`) (i.e.,
  // iteration indices are 0-based). At least one call of this method (with
  // argument 0) must be made.
  //
  // Returns an error status if `num_iterations` is non-positive, if
  // preprocessing was indicated but `PreprocessingComplete()` has not been
  // called, or if successive calls of this method do not pass monotonically
  // ascending `iteration` values of 0, 1, 2, ..., `num_iterations` - 1.
  //
  // Example usage:
  //   int64_t n = 100;  // Number of iterations.
  //   gbbs::IterationProgressReporter iter_reporter(report_progress, n);
  //   for (int i = 0; i < n; ++i) {
  //     // ... do work here for the ith iteration ...
  //     ABSL_CHECK_OK(iter_reporter.IterationComplete(i));
  //   }
  //
  // Note: Use of this method is intended for cases within sequential "for"
  // loops where the `iteration` value is increasing sequentially. For use
  // within `parlay::parallel_for()` loops, you must use
  // `ThreadSafeIterationComplete()` instead.
  absl::Status IterationComplete(int64_t iteration);

  // Indicates that one of the iterations has completed. This a convenience
  // method intended for use at the end of a `parlay::parallel_for()` loop,
  // where the iteration counter is not sequentially increasing. In particular,
  // the `iteration` value passed to this method will be generally increasing
  // but not sequentially or even monotonically so (because the parallel threads
  // will complete their work at different rates).
  //
  // Returns an error status if preprocessing was indicated but
  // `PreprocessingComplete()` has not been called, or if this method has been
  // called more than the `num_iterations` value passed to the constructor.
  //
  // Example usage:
  //   int64_t n = 100;  // Number of iterations.
  //   gbbs::IterationProgressReporter iter_reporter(report_progress, n);
  //   parlay::parallel_for(0, n, [&](size_t i) {
  //     // ... do work here for the ith iteration ...
  //     ABSL_CHECK_OK(iter_reporter.ThreadSafeIterationComplete(i));
  //   });
  //
  // Note: This method is not intended for use within a sequential "for" loop
  // because the `iteration` value will be sequentially increasing. Instead, you
  // should use `IterationComplete()` which expects a sequentially increasing
  // `iteration` value.
  absl::Status ThreadSafeIterationComplete(int64_t iteration);

  // Indicates that all iterations are complete. This is equivalent to calling
  // `IterationComplete(num_iterations)`, but can be used in cases where the
  // algorithm detects that it has converged early and not all expected
  // iterations are required
  //
  // Returns an error status if `num_iterations` is non-positive or if
  // preprocessing was indicated but `PreprocessingComplete()` has not been
  // called.
  absl::Status IterationsDoneEarly();

  // Indicate that the postprocessing step is complete. Returns an error status
  // if `has_postprocessing_` is false or if no call of `IterationComplete()`
  // has been made.
  absl::Status PostprocessingComplete();

 private:
  /*const*/ ReportProgressCallback report_progress_;
  const int64_t num_iterations_;
  const bool has_preprocessing_;
  const bool has_postprocessing_;
  const int64_t total_steps_;
  absl::Mutex mu_;  // Guards `num_steps_completed_` in
                    // `ThreadSafeIterationComplete()`.
  int64_t num_steps_completed_ = 0;

  // Invokes the `report_progress_` callback so as to reflect
  // `num_steps_completed_` worth of progress out of `total_steps_`.
  void ReportProgress();
};

}  // namespace gbbs
