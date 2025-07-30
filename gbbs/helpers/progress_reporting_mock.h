#pragma once

#include <vector>

#include "absl/base/attributes.h"
#include "gbbs/helpers/progress_reporting.h"

namespace gbbs {

// Test-only class for mocking the progress reporting of an algorithm. This
// is intended to be a mix-in class to be combined into a test fixture, like
// this:
//
// class MyAlgorithmTest : public testing::Test, public ReportProgressMock {
//   // ...
// };
//
// TEST_F(MyAlgorithmTest, DoesFoo) {
//   // Invoke method under test.
//   EXPECT_THAT(clusterer_.Cluster(..., MockReportProgress()),
//               IsOkAndHolds(...));
//   // Check that the progress reports were as expected.
//   EXPECT_THAT(GetProgressReports(), ElementsAre(0.0, 1.0));
// }
//
// Note that because `MockReportProgress()` clears any previous progress
// reports, multiple `MockReportProgress()` ... `GetProgressReports()` call
// pairs can be used in a single test.
class ReportProgressMock {
 public:
  // Returns a new ReportProgress `AnyInvocable` that records its successive
  // arguments in `progress_reports_`. Any previously recorded progress reports
  // are cleared.
  ReportProgressT MockReportProgress() ABSL_ATTRIBUTE_LIFETIME_BOUND {
    progress_reports_.clear();
    return [this](float progress) {
      this->progress_reports_.push_back(progress);
    };
  }

  // Returns the arguments to the sequence of progress report calls.
  std::vector<float> GetProgressReports() const { return progress_reports_; }

 private:
  std::vector<float> progress_reports_;
};

}  // namespace gbbs
