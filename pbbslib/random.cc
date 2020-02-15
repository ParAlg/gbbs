#include "random.h"

#include "utilities.h"

namespace pbbs {

random::random(size_t seed) : state(seed) {};
random::random() : state(0) {};
random random::fork(uint64_t i) const {
  return random(hash64(hash64(i+state))); }
random random::next() const { return fork(0);}
size_t random::ith_rand(uint64_t i) const {
  return hash64(i+state);}
size_t random::operator[] (size_t i) const {return ith_rand(i);}
size_t random::rand() { return ith_rand(0);}

}  // namespace pbbs
