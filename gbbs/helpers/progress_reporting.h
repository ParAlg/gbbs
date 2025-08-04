#pragma once

// This library contains utilities for reporting the progress of the execution
// of an algorithm.

#include <utility>

#include "absl/functional/any_invocable.h"
#include "absl/status/status.h"

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
  IterationProgressReporter(ReportProgressCallback report_progress,
                            int num_iterations, bool has_preprocessing = false,
                            bool has_postprocessing = false)
      : report_progress_(std::move(report_progress)),
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
  // iteration indices are 0-based). At least one call of
  // this method (with argument 0) must be made.
  //
  // Returns an error status if `num_iterations` is non-positive, if
  // preprocessing was indicated but `PreprocessingComplete()` has not been
  // called, or if successive calls of this method do not pass monotonically
  // ascending `iteration` values of 1, 2, ..., `num_iterations`.
  absl::Status IterationComplete(int iteration);

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
  const int num_iterations_;
  const bool has_preprocessing_;
  const bool has_postprocessing_;
  const int total_steps_;
  int num_steps_completed_ = 0;

  // Invokes the `report_progress_` callback so as to reflect
  // `num_steps_completed_` worth of progress out of `total_steps_`.
  void ReportProgress();
};

}  // namespace gbbs
