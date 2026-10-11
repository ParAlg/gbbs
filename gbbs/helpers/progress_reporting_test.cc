#include "gbbs/helpers/progress_reporting.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "absl/status/status.h"
#include "gbbs/helpers/progress_reporting_mock.h"
#include "gbbs/helpers/status_macros.h"

namespace gbbs {

namespace {

using ::testing::ElementsAre;

class IterationProgressReporterTest : public ::testing::Test,
                                      public gbbs::ReportProgressMock {};

// Successful complete lifecycle tests without early termination (i.e., where
// all iterations are completed).
using CompleteLifecycleTest = IterationProgressReporterTest;

TEST_F(CompleteLifecycleTest, NoPreprocessingOrPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/4,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.IterationComplete(3));
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/3,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/3,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest, WithPreprocessingAndPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/2,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

// Successful complete lifecycle tests where not all iterations are completed
// and the `IterationsDoneEarly()` method is called.
using CompleteLifecycleEarlyTerminationTest = IterationProgressReporterTest;

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPreprocessingOrPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/8,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.IterationComplete(3));
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 0.25, 0.375, 0.5, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/7,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.IterationComplete(3));
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(0.125, 0.25, 0.375, 0.5, 0.875, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/7,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 0.25, 0.375, 0.5, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest,
       WithPreprocessingAndePostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/6,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_OK(reporter.IterationComplete(2));
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(0.125, 0.25, 0.375, 0.5, 0.875, 1.0));
}

// Successful complete lifecycle tests where NO iterations are completed
// and the `IterationsDoneEarly()` method is called.
using CompleteLifecycleDoneImmediatelyTest = IterationProgressReporterTest;

TEST_F(CompleteLifecycleDoneImmediatelyTest, NoPreprocessingOrPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/8,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(1.0));
}

TEST_F(CompleteLifecycleDoneImmediatelyTest, NoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/7,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.875, 1.0));
}

TEST_F(CompleteLifecycleDoneImmediatelyTest, NoPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/7,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 1.0));
}

TEST_F(CompleteLifecycleDoneImmediatelyTest,
       WithPreprocessingAndPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/6,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 0.875, 1.0));
}

// Tests checking for various error conditions.
using ProcessCompleteTest = IterationProgressReporterTest;
using IterationCompleteTest = IterationProgressReporterTest;
using IterationsDoneEarlyTest = IterationProgressReporterTest;
using PostprocessingCompleteTest = IterationProgressReporterTest;

TEST_F(ProcessCompleteTest, Disallowed) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(reporter.PreprocessingComplete(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       "'PreprocessingComplete()' called on an instance with "
                       "'has_preprocessing' of false."));
}

TEST_F(ProcessCompleteTest, CalledMultipleTimes) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_THAT(reporter.PreprocessingComplete(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       "'PreprocessingComplete()' called more than once."));
}

TEST_F(IterationCompleteTest, PreprocessingCompleteNotCalled) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/true);
  EXPECT_THAT(
      reporter.IterationComplete(0),
      StatusIs(absl::StatusCode::kFailedPrecondition,
               "'PreprocessingComplete()' must be called before any call of "
               "'IterationComplete()'."));
}

TEST_F(IterationCompleteTest, NegativeNumIterations) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/-1,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(
      reporter.IterationComplete(0),
      StatusIs(absl::StatusCode::kInvalidArgument,
               "'num_iterations' must be strictly positive."));
}

TEST_F(IterationCompleteTest, ZeroNumIterations) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/0,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(
      reporter.IterationComplete(0),
      StatusIs(absl::StatusCode::kInvalidArgument,
               "'num_iterations' must be strictly positive."));
}

TEST_F(IterationCompleteTest, IterationTooSmall) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/1,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(reporter.IterationComplete(-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       "'iteration' must be in the half-open interval [0, "
                       "num_iterations)."));
}

TEST_F(IterationCompleteTest, IterationTooLarge) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/2,
                                     /*has_preprocessing=*/false);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_THAT(
      reporter.IterationComplete(2),
      StatusIs(absl::StatusCode::kInvalidArgument,
        "'iteration' must be in the half-open interval [0, "
        "num_iterations)."));
}

TEST_F(IterationCompleteTest, InvalidIterationValueNoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/5,
                                     /*has_preprocessing=*/false);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_THAT(reporter.IterationComplete(3),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       "'iteration' must be monotonically ascending; expected "
                       "2 but got 3."));
}

TEST_F(IterationCompleteTest, InvalidIterationValueWithPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/5,
                                     /*has_preprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_OK(reporter.IterationComplete(1));
  EXPECT_THAT(reporter.IterationComplete(3),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       "'iteration' must be monotonically ascending; expected "
                       "2 but got 3."));
}

TEST_F(IterationsDoneEarlyTest, NegativeNumIterations) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/-1,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(
      reporter.IterationsDoneEarly(),
      StatusIs(absl::StatusCode::kInvalidArgument,
               "'num_iterations' must be strictly positive."));
}

TEST_F(IterationsDoneEarlyTest, ZeroNumIterations) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/0,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(
      reporter.IterationsDoneEarly(),
      StatusIs(absl::StatusCode::kInvalidArgument,
               "'num_iterations' must be strictly positive."));
}

TEST_F(IterationsDoneEarlyTest, PreprocessingCompleteNotCalled) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/true);
  EXPECT_THAT(
      reporter.IterationsDoneEarly(),
      StatusIs(absl::StatusCode::kFailedPrecondition,
               "'PreprocessingComplete()' must be called before any call of "
               "'IterationsDoneEarly()'."));
}

TEST_F(PostprocessingCompleteTest, Disallowed) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.IterationComplete(0));
  EXPECT_THAT(reporter.PostprocessingComplete(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       "'PostprocessingComplete()' called on an instance with "
                       "'has_postprocessing' of false."));
}

TEST_F(PostprocessingCompleteTest, CalledTooEarlyNoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  EXPECT_THAT(reporter.PostprocessingComplete(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       "'PostprocessingComplete()' called before any call of "
                       "IterationComplete()' or 'IterationsDoneEarly()'."));
}

TEST_F(PostprocessingCompleteTest, CalledTooEarlyWithPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback, /*num_iterations=*/1,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  EXPECT_THAT(reporter.PostprocessingComplete(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       "'PostprocessingComplete()' called before any call of "
                       "IterationComplete()' or 'IterationsDoneEarly()'."));
}

}  // namespace

}  // namespace gbbs
