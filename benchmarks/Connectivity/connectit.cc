#include "benchmarks/Connectivity/connectit.h"

#include "pbbslib/assert.h"

// Creates a case on `case_label` and returns `case_label` as a string in the
// case.
#define RETURN_CASE_LABEL(case_label) \
  case case_label: \
    return "case_label";

namespace connectit {

std::string find_to_string(FindOption find_option) {
  switch (find_option) {
    RETURN_CASE_LABEL(find_compress);
    RETURN_CASE_LABEL(find_naive);
    RETURN_CASE_LABEL(find_split);
    RETURN_CASE_LABEL(find_halve);
    RETURN_CASE_LABEL(find_atomic_split);
    RETURN_CASE_LABEL(find_atomic_halve);
  }
  ABORT_INVALID_ENUM(FindOption, find_option);
}

} // namespace connectit
