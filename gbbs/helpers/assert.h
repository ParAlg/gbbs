// Assertion macros.
//
// These are implemented as macros instead of functions so that `__FILE__` and
// `__LINE__` may be used to print the location of a failing assertion.
#pragma once

#include <exception>
#include <iostream>

// Prints out message and terminates the program unconditionally.
//
// Usage examples
// --------------
// If the program 'foo/bar.cc' hits
//     ABORT("A bad thing happened");
// on line 45, then the program will terminate at runtime with the message
//     foo/bar.cc:45: Abort: A bad thing happened
// Moreover, to some extent, the input message can be treated like a stream:
//     int x = 0;
//     ABORT("Unexpected value of x: " << x);
#define ABORT(message)                                                 \
  do {                                                                 \
    std::cerr << __FILE__ << ":" << __LINE__ << ": Abort: " << message \
              << std::endl;                                            \
    std::terminate();                                                  \
  } while (false)

// Prints out enum value and terminates the program.
//
// This is intended for use after a switch statement on an enum that handles all
// possible enum cases by returning. Using this macro will stop GCC's
// -Wreturn-type flag from issuing a warning about this situation.
// Example:
//     enum class OneOrTwo { kOne, kTwo };
//
//     int ReturnOneOrTwo(OneOrTwo one_or_two) {
//       switch (one_or_two) {
//         case OneOrTwo::kOne:
//           return 1;
//         case OneOrTwo::kTwo:
//           return 2;
//       }
//       // If we didn't put this ABORT statement below, GCC would complain:
//       //   warning: control reaches end of non-void function [-Wreturn-type]
//       // We shouldn't reach this point unless someone makes a weird function
//       // call like `ReturnOneOrTwo(static_cast<OneOrTwo>(3))`.
//       ABORT_INVALID_ENUM(OneOrTwo, one_or_two)
//     }
#define ABORT_INVALID_ENUM(EnumType, enum_value)                       \
  ABORT("Unexpected " #EnumType " value: "                             \
        << static_cast<typename std::underlying_type<EnumType>::type>( \
            enum_value));

// Asserts on a condition, printing an error and terminating if the condition is
// false.
//
// This is more intended for debugging and validating internal programming logic
// than for detecting user errors. We may choose to replace `ASSERT` with a
// no-op in optimized builds.
//
// This is overloaded and can take either one or two arguments.
//
// Arguments
// ---------
// condition: bool
//   The condition on which to assert.
// (optional) message: string
//   An explanatory message to print out when the condition fails.
//
// Usage examples
// --------------
// If the program 'foo/bar.cc' hits
//     ASSERT(0 == 1, "Zero should be equal to one");
// on line 45, then the program will terminate at runtime with the message
//     foo/bar.cc:45: Failed assertion `0 == 1`: Zero should be equal to one
// The input message may also be omitted:
//     ASSERT(0 == 1);
// Moreover, to some extent, the input message can be treated like a stream:
//     int x = 0;
//     ABORT(x > 0, "x must be positive, was " << x << " instead");
#define ASSERT(...) _GET_MACRO(__VA_ARGS__, _ASSERT2, _ASSERT1)(__VA_ARGS__)

//////////////
// Internal //
//////////////

// Trickery from https://stackoverflow.com/a/11763277/4865149 to support
// overloading macros by number of arguments.
// If we need to do this kind of overloading more often, consider
// importing Boost and replacing this macro with BOOST_PP_OVERLOAD.
#define _GET_MACRO(_1, _2, NAME, ...) NAME

// Implementation of ASSERT with one argument.
#define _ASSERT1(condition)                                            \
  do {                                                                 \
    if (!(condition)) {                                                \
      std::cerr << __FILE__ << ":" << __LINE__                         \
                << ": Failed assertion `" #condition "`" << std::endl; \
      std::terminate();                                                \
    }                                                                  \
  } while (false)

// Implementation of ASSERT with two arguments.
#define _ASSERT2(condition, message)                                  \
  do {                                                                \
    if (!(condition)) {                                               \
      std::cerr << __FILE__ << ":" << __LINE__                        \
                << ": Failed assertion `" #condition "`: " << message \
                << std::endl;                                         \
      std::terminate();                                               \
    }                                                                 \
  } while (false)
