// Error handling across the C++/R boundary.
//
// Rf_error() longjmps. Called while C++ objects with destructors are alive
// (std::vector, std::string, an active catch handler and its exception
// object) it skips their destructors and leaks them, and leaving a catch
// handler by longjmp also leaves the C++ runtime's caught-exception state
// inconsistent. Native code therefore reports failures by throwing
// (compressor_fail()), and each .Call entry point runs its body under
// compressor_guard(), which copies the message into a plain buffer, lets the
// handler and every C++ object of the body be destroyed, and only then
// raises the R error.
#ifndef COMPRESSOR_NATIVE_ERROR_H
#define COMPRESSOR_NATIVE_ERROR_H

#include <R.h>
#include <Rinternals.h>

#include <cstdarg>
#include <cstdio>
#include <exception>
#include <stdexcept>

[[noreturn]] inline void compressor_fail(const char* format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  std::vsnprintf(buffer, sizeof(buffer), format, args);
  va_end(args);
  throw std::runtime_error(buffer);
}

template <typename F>
SEXP compressor_guard(F&& body, const char* fallback) {
  char message[1024];
  message[0] = '\0';
  try {
    return body();
  } catch (const std::exception& failure) {
    std::snprintf(message, sizeof(message), "%s", failure.what());
  } catch (...) {
    std::snprintf(message, sizeof(message), "%s", fallback);
  }
  // Outside the handler: no C++ object of `body` and no exception is alive.
  Rf_error("%s", message);
  return R_NilValue;  // not reached
}

#endif
