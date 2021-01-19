#pragma once

#include <string>

namespace pbbs {

struct timer {
  double total_time;
  double last_time;
  bool on;
  std::string name;

  timer(std::string name = "# PBBS time", bool _start = true);

  double get_time() const;
  void start();
  double stop();
  void reset();
  double get_total() const;
  double get_next();
  void report(double time, const std::string& str) const;
  void total();
  void reportTotal(const std::string& str) const;
  void next(const std::string& str);
};

}  // namespace pbbs
