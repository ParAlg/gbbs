#include "gbbs/helpers/progress_reporting.h"

#include <cstddef>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "absl/status/status.h"
#include "gbbs/helpers/progress_reporting_mock.h"
#include "gbbs/helpers/status_macros.h"
#include "parlay/parallel.h"

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
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/4,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  for (int iteration = 0; iteration < 4; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest,
       NoPreprocessingOrPostprocessing_PeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/30,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false,
                                     /*max_iteration_callbacks=*/4);
  for (int iteration = 0; iteration < 30; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(8.0 / 30.0, 16.0 / 30.0, 24.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest,
       NoPreprocessingOrPostprocessing_WithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/4,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  parlay::parallel_for(0, 4, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest,
       NoPreprocessingOrPostprocessing_WithThreadSafetyAndPeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/30,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false,
                                     /*max_iteration_callbacks=*/4);
  parlay::parallel_for(0, 30, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(8.0 / 30.0, 16.0 / 30.0, 24.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/3,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  for (int iteration = 0; iteration < 3; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPreprocessing_PeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/29,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true,
                                     /*max_iteration_callbacks=*/4);
  for (int iteration = 0; iteration < 29; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(8.0 / 30.0, 16.0 / 30.0,
                                                24.0 / 30.0, 29.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPreprocessing_WithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/3,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  parlay::parallel_for(0, 3, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest,
       NoPreprocessing_WithThreadSafetyAndPeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/29,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true,
                                     /*max_iteration_callbacks=*/4);
  parlay::parallel_for(0, 29, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(8.0 / 30.0, 16.0 / 30.0,
                                                24.0 / 30.0, 29.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/3,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 3; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPostprocessing_PeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/29,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false,
                                     /*max_iteration_callbacks=*/4);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 29; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_THAT(GetProgressReports(), ElementsAre(1.0 / 30.0, 9.0 / 30.0,
                                                17.0 / 30.0, 25.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest, NoPostprocessing_WithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/3,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  parlay::parallel_for(0, 3, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest,
       NoPostprocessing_WithThreadSafetyAndPeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/29,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false,
                                     /*max_iteration_callbacks=*/4);
  EXPECT_OK(reporter.PreprocessingComplete());
  parlay::parallel_for(0, 29, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(GetProgressReports(), ElementsAre(1.0 / 30.0, 9.0 / 30.0,
                                                17.0 / 30.0, 25.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest, WithPreprocessingAndPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/2,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 2; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(CompleteLifecycleTest,
       WithPreprocessingAndPostprocessing_PeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/28,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true,
                                     /*max_iteration_callbacks=*/4);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 28; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(1.0 / 30.0, 8.0 / 30.0, 15.0 / 30.0, 22.0 / 30.0,
                          29.0 / 30.0, 1.0));
}

TEST_F(CompleteLifecycleTest,
       WithPreprocessingAndPostprocessing_WithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/2,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  parlay::parallel_for(0, 2, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.25, 0.5, 0.75, 1.0));
}

TEST_F(
    CompleteLifecycleTest,
    WithPreprocessingAndPostprocessing_WithThreadSafetyAndPeriodicCallbacks) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/28,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true,
                                     /*max_iteration_callbacks=*/4);
  EXPECT_OK(reporter.PreprocessingComplete());
  parlay::parallel_for(0, 28, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(1.0 / 30.0, 8.0 / 30.0, 15.0 / 30.0, 22.0 / 30.0,
                          29.0 / 30.0, 1.0));
}

// Tests that `IterationProgressReporter` supports alternating sequential
// and parallel iterations.
TEST_F(CompleteLifecycleTest,
       WithAlternatingSequentialAndParallelIterations) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/8,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 2; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  parlay::parallel_for(2, 6, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  for (int iteration = 6; iteration < 8; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0));
}

// Successful complete lifecycle tests where not all iterations are completed
// and the `IterationsDoneEarly()` method is called.
using CompleteLifecycleEarlyTerminationTest = IterationProgressReporterTest;

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPreprocessingOrPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/8,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/false);
  for (int iteration = 0; iteration < 4; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 0.25, 0.375, 0.5, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPreprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/7,
                                     /*has_preprocessing=*/false,
                                     /*has_postprocessing=*/true);
  for (int iteration = 0; iteration < 4; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_OK(reporter.PostprocessingComplete());
  EXPECT_THAT(GetProgressReports(),
              ElementsAre(0.125, 0.25, 0.375, 0.5, 0.875, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest, NoPostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/7,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/false);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 3; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
  EXPECT_OK(reporter.IterationsDoneEarly());
  EXPECT_THAT(GetProgressReports(), ElementsAre(0.125, 0.25, 0.375, 0.5, 1.0));
}

TEST_F(CompleteLifecycleEarlyTerminationTest,
       WithPreprocessingAndePostprocessing) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/6,
                                     /*has_preprocessing=*/true,
                                     /*has_postprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  for (int iteration = 0; iteration < 3; ++iteration) {
    EXPECT_OK(reporter.IterationComplete(iteration));
  }
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

TEST_F(IterationCompleteTest, PreprocessingCompleteNotCalledWithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/1,
                                     /*has_preprocessing=*/true);
  EXPECT_THAT(
      reporter.ThreadSafeIterationComplete(0),
      StatusIs(absl::StatusCode::kFailedPrecondition,
               "'PreprocessingComplete()' must be called before any call of "
               "'ThreadSafeIterationComplete()'."));
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
  EXPECT_THAT(reporter.IterationComplete(0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       "'num_iterations' must be strictly positive."));
}

TEST_F(IterationCompleteTest, ZeroNumIterationsWithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/0,
                                     /*has_preprocessing=*/false);
  EXPECT_THAT(
      reporter.ThreadSafeIterationComplete(0),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          "'ThreadSafeIterationComplete()' (and possibly "
          "'IterationComplete()') called more than 'num_iterations' times."));
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

TEST_F(IterationCompleteTest, IterationTooLargeWithThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/2,
                                     /*has_preprocessing=*/false);
  parlay::parallel_for(0, 2, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(
      reporter.ThreadSafeIterationComplete(2),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          "'ThreadSafeIterationComplete()' (and possibly "
          "'IterationComplete()') called more than 'num_iterations' times."));
}

TEST_F(IterationCompleteTest,
       IterationTooLargeWithPreprocessingAndThreadSafety) {
  ReportProgressCallback mock_report_progress_callback =
      MockReportProgressCallback();
  IterationProgressReporter reporter(mock_report_progress_callback,
                                     /*num_iterations=*/2,
                                     /*has_preprocessing=*/true);
  EXPECT_OK(reporter.PreprocessingComplete());
  parlay::parallel_for(0, 2, [&](size_t i) {
    EXPECT_OK(reporter.ThreadSafeIterationComplete(i));
  });
  EXPECT_THAT(
      reporter.ThreadSafeIterationComplete(2),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          "'ThreadSafeIterationComplete()' (and possibly "
          "'IterationComplete()') called more than 'num_iterations' times."));
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
