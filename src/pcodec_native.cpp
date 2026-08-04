#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#ifdef COMPRESSOR_NATIVE_PCODEC
#include "pcodec_native.h"

namespace {
constexpr unsigned char kPcoTypeU32 = 1;
constexpr unsigned char kPcoTypeU16 = 7;
constexpr unsigned char kPcoTypeU8 = 10;

void check_status(CompressorPcoError status, const char* operation) {
  if (status != COMPRESSOR_PCO_SUCCESS) {
    Rf_error("native Pcodec %s failed (error code %d)", operation,
             static_cast<int>(status));
  }
}

template <typename T>
SEXP compress_integer(SEXP input, SEXP level, SEXP page_n, unsigned char dtype) {
  if (TYPEOF(input) != INTSXP) Rf_error("native Pcodec expects an integer vector");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  std::vector<T> values(n);
  const int* source = INTEGER(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (source[i] == NA_INTEGER || source[i] < 0 ||
        static_cast<std::uint64_t>(source[i]) >
          static_cast<std::uint64_t>(std::numeric_limits<T>::max())) {
      Rf_error("native Pcodec input contains an out-of-range integer");
    }
    values[i] = static_cast<T>(source[i]);
  }
  CompressorPcoChunkConfig config{
    static_cast<unsigned int>(Rf_asInteger(level)),
    static_cast<std::size_t>(std::max(0, Rf_asInteger(page_n)))
  };
  const std::size_t capacity = compressor_pco_guarantee_file_size(n, dtype);
  if (!capacity) Rf_error("native Pcodec returned zero output capacity");
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(capacity)));
  std::size_t written = 0;
  check_status(compressor_pco_compress_into(
      values.data(), n, dtype, &config, RAW(output), capacity, &written),
    "compression");
  if (written > capacity) Rf_error("native Pcodec wrote beyond its output buffer");
  if (written != capacity) {
    SEXP trimmed = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(written)));
    if (written) std::memcpy(RAW(trimmed), RAW(output), written);
    UNPROTECT(2);
    return trimmed;
  }
  UNPROTECT(1);
  return output;
}

template <typename T>
SEXP decompress_integer(SEXP compressed, SEXP n, unsigned char dtype) {
  if (TYPEOF(compressed) != RAWSXP) Rf_error("compressed must be a raw vector");
  const R_xlen_t expected = Rf_asInteger(n);
  if (expected < 0) Rf_error("n must be non-negative");
  std::vector<T> values(static_cast<std::size_t>(expected));
  std::size_t written = 0;
  check_status(compressor_pco_decompress_into(
      RAW(compressed), static_cast<std::size_t>(XLENGTH(compressed)), dtype,
      values.data(), static_cast<std::size_t>(expected), &written),
    "decompression");
  if (written != static_cast<std::size_t>(expected)) {
    Rf_error("native Pcodec returned %zu values; expected %td", written, expected);
  }
  SEXP output = PROTECT(Rf_allocVector(INTSXP, expected));
  int* target = INTEGER(output);
  for (std::size_t i = 0; i < written; ++i) {
    if (static_cast<std::uint64_t>(values[i]) >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      UNPROTECT(1);
      Rf_error("native Pcodec value does not fit R's integer type");
    }
    target[i] = static_cast<int>(values[i]);
  }
  UNPROTECT(1);
  return output;
}

SEXP compress_numeric_u32(SEXP input, SEXP level, SEXP page_n) {
  if (TYPEOF(input) != REALSXP) Rf_error("native Pcodec uint32 input must be numeric");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  std::vector<std::uint32_t> values(n);
  const double* source = REAL(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (!R_FINITE(source[i]) || source[i] < 0.0 ||
        source[i] != std::floor(source[i]) ||
        source[i] > 4294967295.0) {
      Rf_error("native Pcodec uint32 input contains an invalid value");
    }
    values[i] = static_cast<std::uint32_t>(source[i]);
  }
  CompressorPcoChunkConfig config{
    static_cast<unsigned int>(Rf_asInteger(level)),
    static_cast<std::size_t>(std::max(0, Rf_asInteger(page_n)))
  };
  const std::size_t capacity = compressor_pco_guarantee_file_size(n, kPcoTypeU32);
  if (!capacity) Rf_error("native Pcodec returned zero output capacity");
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(capacity)));
  std::size_t written = 0;
  check_status(compressor_pco_compress_into(
      values.data(), n, kPcoTypeU32, &config, RAW(output), capacity, &written),
    "uint32 compression");
  if (written > capacity) Rf_error("native Pcodec wrote beyond its output buffer");
  if (written != capacity) {
    SEXP trimmed = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(written)));
    if (written) std::memcpy(RAW(trimmed), RAW(output), written);
    UNPROTECT(2);
    return trimmed;
  }
  UNPROTECT(1);
  return output;
}

SEXP decompress_numeric_u32(SEXP compressed, SEXP n) {
  if (TYPEOF(compressed) != RAWSXP) Rf_error("compressed must be a raw vector");
  const R_xlen_t expected = Rf_asInteger(n);
  if (expected < 0) Rf_error("n must be non-negative");
  std::vector<std::uint32_t> values(static_cast<std::size_t>(expected));
  std::size_t written = 0;
  check_status(compressor_pco_decompress_into(
      RAW(compressed), static_cast<std::size_t>(XLENGTH(compressed)), kPcoTypeU32,
      values.data(), static_cast<std::size_t>(expected), &written),
    "uint32 decompression");
  if (written != static_cast<std::size_t>(expected)) {
    Rf_error("native Pcodec returned %zu values; expected %td", written, expected);
  }
  SEXP output = PROTECT(Rf_allocVector(REALSXP, expected));
  double* target = REAL(output);
  for (std::size_t i = 0; i < written; ++i) target[i] = values[i];
  UNPROTECT(1);
  return output;
}
}

extern "C" SEXP compressor_pcodec_native_available() {
  return Rf_ScalarLogical(TRUE);
}

extern "C" SEXP compressor_pcodec_compress_u8(SEXP input, SEXP level, SEXP page_n) {
  return compress_integer<std::uint8_t>(input, level, page_n, kPcoTypeU8);
}

extern "C" SEXP compressor_pcodec_compress_u16(SEXP input, SEXP level, SEXP page_n) {
  return compress_integer<std::uint16_t>(input, level, page_n, kPcoTypeU16);
}

extern "C" SEXP compressor_pcodec_compress_u32(SEXP input, SEXP level, SEXP page_n) {
  return compress_numeric_u32(input, level, page_n);
}

extern "C" SEXP compressor_pcodec_decompress_u8(SEXP compressed, SEXP n) {
  return decompress_integer<std::uint8_t>(compressed, n, kPcoTypeU8);
}

extern "C" SEXP compressor_pcodec_decompress_u16(SEXP compressed, SEXP n) {
  return decompress_integer<std::uint16_t>(compressed, n, kPcoTypeU16);
}

extern "C" SEXP compressor_pcodec_decompress_u32(SEXP compressed, SEXP n) {
  return decompress_numeric_u32(compressed, n);
}

#else

extern "C" SEXP compressor_pcodec_native_available() {
  return Rf_ScalarLogical(FALSE);
}

extern "C" SEXP compressor_pcodec_compress_u8(SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_compress_u16(SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_compress_u32(SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_decompress_u8(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_decompress_u16(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_decompress_u32(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}

#endif
