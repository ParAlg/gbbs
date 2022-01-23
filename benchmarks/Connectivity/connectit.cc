#include "benchmarks/Connectivity/connectit.h"

#include "gbbs/helpers/assert.h"

namespace gbbs {

// Creates a case on `case_label` and returns `case_label` as a string in the
// case.
#define RETURN_CASE_LABEL(case_label) \
  case case_label:                    \
    return #case_label;

namespace connectit {

std::string find_to_string(const FindOption find_option) {
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

std::string splice_to_string(const SpliceOption splice_option) {
  switch (splice_option) {
    RETURN_CASE_LABEL(split_atomic_one);
    RETURN_CASE_LABEL(halve_atomic_one);
    RETURN_CASE_LABEL(splice_simple);
    RETURN_CASE_LABEL(splice_atomic);
  }
  ABORT_INVALID_ENUM(SpliceOption, splice_option);
}

std::string unite_to_string(const UniteOption unite_option) {
  switch (unite_option) {
    RETURN_CASE_LABEL(unite);
    RETURN_CASE_LABEL(unite_early);
    RETURN_CASE_LABEL(unite_nd);
    RETURN_CASE_LABEL(unite_rem_lock);
    RETURN_CASE_LABEL(unite_rem_cas);
  }
  ABORT_INVALID_ENUM(UniteOption, unite_option);
}

std::string sampling_to_string(const SamplingOption sampling_option) {
  switch (sampling_option) {
    RETURN_CASE_LABEL(sample_kout);
    RETURN_CASE_LABEL(sample_bfs);
    RETURN_CASE_LABEL(sample_ldd);
    RETURN_CASE_LABEL(no_sampling);
  }
  ABORT_INVALID_ENUM(SamplingOption, sampling_option);
}

std::string jayanti_find_to_string(const JayantiFindOption find_option) {
  switch (find_option) {
    RETURN_CASE_LABEL(find_twotrysplit);
    RETURN_CASE_LABEL(find_simple);
  }
  ABORT_INVALID_ENUM(JayantiFindOption, find_option);
}

std::string connect_to_string(const LiuTarjanConnectOption connect_option) {
  switch (connect_option) {
    RETURN_CASE_LABEL(simple_connect);
    RETURN_CASE_LABEL(parent_connect);
    RETURN_CASE_LABEL(extended_connect);
  }
  ABORT_INVALID_ENUM(LiuTarjanConnectOption, connect_option);
}

std::string update_to_string(const LiuTarjanUpdateOption update_option) {
  switch (update_option) {
    RETURN_CASE_LABEL(simple_update);
    RETURN_CASE_LABEL(root_update);
  }
  ABORT_INVALID_ENUM(LiuTarjanUpdateOption, update_option);
}

std::string shortcut_to_string(const LiuTarjanShortcutOption shortcut_option) {
  switch (shortcut_option) {
    RETURN_CASE_LABEL(shortcut);
    RETURN_CASE_LABEL(full_shortcut);
  }
  ABORT_INVALID_ENUM(LiuTarjanShortcutOption, shortcut_option);
}

std::string alter_to_string(const LiuTarjanAlterOption alter_option) {
  switch (alter_option) {
    RETURN_CASE_LABEL(alter);
    RETURN_CASE_LABEL(no_alter);
  }
  ABORT_INVALID_ENUM(LiuTarjanAlterOption, alter_option);
}

std::string uf_options_to_string(const SamplingOption sampling_option,
                                 const FindOption find_option,
                                 const UniteOption unite_option) {
  return "uf; sample=" + sampling_to_string(sampling_option) + "; unite=" +
         unite_to_string(unite_option) + "; find=" +
         find_to_string(find_option);
}

std::string uf_options_to_string(const SamplingOption sampling_option,
                                 const FindOption find_option,
                                 const UniteOption unite_option,
                                 const SpliceOption splice_option) {
  return "uf; sample=" + sampling_to_string(sampling_option) + "; unite=" +
         unite_to_string(unite_option) + "; find=" +
         find_to_string(find_option) + +"; splice=" +
         splice_to_string(splice_option);
}

std::string jayanti_options_to_string(const SamplingOption sampling_option,
                                      const JayantiFindOption find_option) {
  return "jayanti; sample=" + sampling_to_string(sampling_option) + "; find=" +
         jayanti_find_to_string(find_option);
}

std::string liu_tarjan_options_to_string(
    const SamplingOption sampling_option,
    const LiuTarjanConnectOption connect_option,
    const LiuTarjanUpdateOption update_option,
    const LiuTarjanShortcutOption shortcut_option,
    const LiuTarjanAlterOption alter_option) {
  return "liu_tarjan; sample=" + sampling_to_string(sampling_option) +
         "; connect=" + connect_to_string(connect_option) + "; update=" +
         update_to_string(update_option) + "; shortcut=" +
         shortcut_to_string(shortcut_option) + "; alter=" +
         alter_to_string(alter_option);
}

}  // namespace connectit
}  // namespace gbbs
