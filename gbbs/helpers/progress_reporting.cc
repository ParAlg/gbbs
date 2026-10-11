#include "gbbs/helpers/progress_reporting.h"

#include <algorithm>
#include <cstdint>

#include "absl/base/attributes.h"
#include "absl/functional/function_ref.h"
#include "absl/log/absl_check.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/synchronization/mutex.h"

namespace gbbs {

IterationProgressReporter::IterationProgressReporter(
    absl::FunctionRef<void(float)> report_progress
        ABSL_ATTRIBUTE_LIFETIME_BOUND,
    int64_t num_iterations, bool has_preprocessing, bool has_postprocessing,
    int64_t max_iteration_callbacks)
    : report_progress_(report_progress),
      num_iterations_(num_iterations),
      has_preprocessing_(has_preprocessing),
      has_postprocessing_(has_postprocessing),
      total_steps_(num_iterations + (has_preprocessing ? 1 : 0) +
                   (has_postprocessing ? 1 : 0)),
      // The callback period is `ceil(num_iterations / max_callbacks)`.
      // The `std::max(1L, ...)` is needed to prevent a callback period of 0
      // in the corner case where `num_iterations` is 0.
      callback_period_(
          std::max(1L, (num_iterations + max_iteration_callbacks - 1) /
                           max_iteration_callbacks)) {
  ABSL_CHECK_GT(max_iteration_callbacks, 0);
}

absl::Status IterationProgressReporter::PreprocessingComplete() {
  if (!has_preprocessing_) {
    return absl::FailedPreconditionError(
        "'PreprocessingComplete()' called on an instance with "
        "'has_preprocessing' of false.");
  }
  if (num_steps_completed_ != 0) {
    return absl::FailedPreconditionError(
        "'PreprocessingComplete()' called more than once.");
  }
  ++num_steps_completed_;
  ReportProgress();
  return absl::OkStatus();
}

bool IterationProgressReporter::PerformCallback(
    const int64_t iteration_num) const {
  const int64_t completed_iterations = iteration_num + 1;
  return completed_iterations % callback_period_ == 0 ||
         completed_iterations == num_iterations_;  // final iteration
}

absl::Status IterationProgressReporter::IterationComplete(int64_t iteration) {
  if (has_preprocessing_ && num_steps_completed_ == 0) {
    return absl::FailedPreconditionError(
        "'PreprocessingComplete()' must be called before any call of "
        "'IterationComplete()'.");
  }
  if (num_iterations_ <= 0) {
    return absl::InvalidArgumentError(
        "'num_iterations' must be strictly positive.");
  }
  if (iteration < 0 || iteration >= num_iterations_) {
    return absl::InvalidArgumentError(
        "'iteration' must be in the half-open interval [0, num_iterations).");
  }

  const int64_t expected_iteration =
      num_steps_completed_ - (has_preprocessing_ ? 1 : 0);
  if (iteration != expected_iteration) {
    return absl::InvalidArgumentError(
        absl::StrCat("'iteration' must be monotonically ascending; expected ",
                     expected_iteration, " but got ", iteration, "."));
  }

  ++num_steps_completed_;
  if (PerformCallback(expected_iteration)) ReportProgress();
  return absl::OkStatus();
}

absl::Status IterationProgressReporter::ThreadSafeIterationComplete(
    const int64_t iteration) {
  // Exit early if not reporting progress for the given `iteration`.
  // This helps avoid lock contention on `mu_`.
  if (!PerformCallback(iteration)) return absl::OkStatus();

  absl::MutexLock lock(mu_);
  if (has_preprocessing_ && num_steps_completed_ == 0) {
    return absl::FailedPreconditionError(
        "'PreprocessingComplete()' must be called before any call of "
        "'ThreadSafeIterationComplete()'.");
  }
  if (iteration >= num_iterations_) {
    return absl::InvalidArgumentError(
        "'ThreadSafeIterationComplete()' (and possibly 'IterationComplete()') "
        "called more than 'num_iterations' times.");
  }
  // Advance `num_steps_completed_` and perform the callback.
  // The `std::min()` is needed to prevent over-reporting in the case where
  // `num_iterations_` is not a multiple of `callback_period_`.
  const int64_t num_iterations_completed =
      num_steps_completed_ - (has_preprocessing_ ? 1 : 0);
  num_steps_completed_ +=
      std::min(callback_period_, num_iterations_ - num_iterations_completed);
  ReportProgress();
  return absl::OkStatus();
}

absl::Status IterationProgressReporter::IterationsDoneEarly() {
  if (has_preprocessing_ && num_steps_completed_ == 0) {
    return absl::FailedPreconditionError(
        "'PreprocessingComplete()' must be called before any call of "
        "'IterationsDoneEarly()'.");
  }
  num_steps_completed_ = (num_iterations_ - 1) + (has_preprocessing_ ? 1 : 0);
  return IterationComplete(num_iterations_ - 1);
}

absl::Status IterationProgressReporter::PostprocessingComplete() {
  if (!has_postprocessing_) {
    return absl::FailedPreconditionError(
        "'PostprocessingComplete()' called on an instance with "
        "'has_postprocessing' of false.");
  }
  if (num_steps_completed_ <= (has_preprocessing_ ? 1 : 0)) {
    return absl::FailedPreconditionError(
        "'PostprocessingComplete()' called before any call of "
        "IterationComplete()' or 'IterationsDoneEarly()'.");
  }
  num_steps_completed_ = total_steps_;
  ReportProgress();
  return absl::OkStatus();
}

void IterationProgressReporter::ReportProgress() {
  report_progress_(static_cast<float>(num_steps_completed_) / total_steps_);
}

}  // namespace gbbs
