#include "gbbs/helpers/progress_reporting.h"

#include "absl/status/status.h"
#include "absl/strings/str_cat.h"

namespace gbbs {

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

absl::Status IterationProgressReporter::IterationComplete(int iteration) {
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

  const int expected_iteration =
      num_steps_completed_ - (has_preprocessing_ ? 1 : 0);
  if (iteration != expected_iteration) {
    return absl::InvalidArgumentError(
        absl::StrCat("'iteration' must be monotonically ascending; expected ",
                     expected_iteration, " but got ", iteration, "."));
  }

  ++num_steps_completed_;
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
