#pragma once

template <class T>
struct maybe {
  T value;
  bool valid;

  maybe(T v, bool u) : value(v) { valid = u; }
  maybe(T v) : value(v) { valid = true; }
  maybe() { valid = false; }

  bool operator!() const { return !valid; }
  operator bool() const { return valid; };
  T& operator*() { return value; }
};
