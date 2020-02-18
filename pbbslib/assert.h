// Assertion macros.
//
// These are implemented as macros instead of functions so that `__FILE__` and
// `__LINE__` may be used to determine the location of a failing assertion.
#pragma once

#include <exception>
#include <iostream>

// Prints out message and terminates the program unconditionally.
#define ABORT(message) \
  do { \
    std::cerr << __FILE__ << ":" << __LINE__ << ": Abort: "  \
         << message << std::endl; \
    std::terminate(); \
  } while (false)

// Asserts on a condition, printing an error and terminating if the condition is
// false.
//
// This is overloaded and can take either one or two arguments.
//
// Arguments
// ---------
//   condition: bool
//     The condition on which to assert.
//   (optional) message: string
//     An explanatory message to print out when the condition fails.
#define ASSERT(...) _GET_MACRO(__VA_ARGS__, _ASSERT2, _ASSERT1)(__VA_ARGS__)


//////////////
// Internal //
//////////////

// Trickery from https://stackoverflow.com/a/11763277/4865149 to support
// overloading macros by number of arguments.
// If we need to do this kind of overloading more often, consider
// importing Boost and replacing this macro with BOOST_PP_OVERLOAD.
#define _GET_MACRO(_1,_2,NAME,...) NAME

// Implementation of ASSERT with one argument.
#define _ASSERT1(condition) \
  do { \
    if (!(condition)) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": Failed assertion `"  \
          #condition "`" << std::endl; \
      std::terminate(); \
    } \
  } while (false)

// Implementation of ASSERT with two arguments.
#define _ASSERT2(condition, message) \
  do { \
    if (!(condition)) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": Failed assertion `"  \
          #condition "`: " << message << std::endl; \
      std::terminate(); \
    } \
  } while (false)
