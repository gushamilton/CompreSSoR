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
#include <atomic>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef COMPRESSOR_NATIVE_PCODEC
#include "pcodec_native.h"
#include "native_error.h"

namespace {
constexpr unsigned char kPcoTypeU32 = 1;
constexpr unsigned char kPcoTypeU16 = 7;
constexpr unsigned char kPcoTypeU8 = 10;

void check_status(CompressorPcoError status, const char* operation) {
  if (status != COMPRESSOR_PCO_SUCCESS) {
    compressor_fail("native Pcodec %s failed (error code %d)", operation,
             static_cast<int>(status));
  }
}

template <typename T>
SEXP compress_integer(SEXP input, SEXP level, SEXP page_n, unsigned char dtype) {
  if (TYPEOF(input) != INTSXP) compressor_fail("native Pcodec expects an integer vector");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  if (n == 0) return Rf_allocVector(RAWSXP, 0);
  std::vector<T> values(n);
  const int* source = INTEGER(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (source[i] == NA_INTEGER || source[i] < 0 ||
        static_cast<std::uint64_t>(source[i]) >
          static_cast<std::uint64_t>(std::numeric_limits<T>::max())) {
      compressor_fail("native Pcodec input contains an out-of-range integer");
    }
    values[i] = static_cast<T>(source[i]);
  }
  CompressorPcoChunkConfig config{
    static_cast<unsigned int>(Rf_asInteger(level)),
    static_cast<std::size_t>(std::max(0, Rf_asInteger(page_n)))
  };
  const std::size_t capacity = compressor_pco_guarantee_file_size(n, dtype);
  if (!capacity) compressor_fail("native Pcodec returned zero output capacity");
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(capacity)));
  std::size_t written = 0;
  check_status(compressor_pco_compress_into(
      values.data(), n, dtype, &config, RAW(output), capacity, &written),
    "compression");
  if (written > capacity) compressor_fail("native Pcodec wrote beyond its output buffer");
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
  if (TYPEOF(compressed) != RAWSXP) compressor_fail("compressed must be a raw vector");
  const R_xlen_t expected = Rf_asInteger(n);
  if (expected < 0) compressor_fail("n must be non-negative");
  if (expected == 0) return Rf_allocVector(INTSXP, 0);
  std::vector<T> values(static_cast<std::size_t>(expected));
  std::size_t written = 0;
  check_status(compressor_pco_decompress_into(
      RAW(compressed), static_cast<std::size_t>(XLENGTH(compressed)), dtype,
      values.data(), static_cast<std::size_t>(expected), &written),
    "decompression");
  if (written != static_cast<std::size_t>(expected)) {
    compressor_fail("native Pcodec returned %zu values; expected %td", written, expected);
  }
  SEXP output = PROTECT(Rf_allocVector(INTSXP, expected));
  int* target = INTEGER(output);
  for (std::size_t i = 0; i < written; ++i) {
    if (static_cast<std::uint64_t>(values[i]) >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      UNPROTECT(1);
      compressor_fail("native Pcodec value does not fit R's integer type");
    }
    target[i] = static_cast<int>(values[i]);
  }
  UNPROTECT(1);
  return output;
}


// Run fn(task) for task in [0, tasks) on up to `threads` std::threads.  The
// callable must not touch the R API.  The first exception is rethrown on the
// calling thread after every worker has joined.
template <typename F>
void run_parallel_tasks(std::size_t tasks, std::size_t threads, F fn) {
  if (!tasks) return;
  const std::size_t workers = std::max<std::size_t>(1, std::min(threads, tasks));
  if (workers == 1) {
    for (std::size_t task = 0; task < tasks; ++task) fn(task);
    return;
  }
  std::atomic<std::size_t> next(0);
  std::mutex error_mutex;
  std::exception_ptr first_error;
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (std::size_t worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&]() {
      try {
        while (true) {
          const std::size_t task = next.fetch_add(1, std::memory_order_relaxed);
          if (task >= tasks) break;
          fn(task);
        }
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!first_error) first_error = std::current_exception();
      }
    });
  }
  for (std::thread& worker : pool) worker.join();
  if (first_error) std::rethrow_exception(first_error);
}

std::size_t thread_count_arg(SEXP threads) {
  const int value = Rf_asInteger(threads);
  return value == NA_INTEGER || value < 1 ? 1 : static_cast<std::size_t>(value);
}

SEXP raw_list(const std::vector<std::vector<unsigned char>>& blobs) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(blobs.size())));
  for (std::size_t i = 0; i < blobs.size(); ++i) {
    SEXP blob = Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(blobs[i].size()));
    SET_VECTOR_ELT(out, static_cast<R_xlen_t>(i), blob);
    if (!blobs[i].empty()) std::memcpy(RAW(blob), blobs[i].data(), blobs[i].size());
  }
  UNPROTECT(1);
  return out;
}

// Compress contiguous `block_rows`-row frames of one stream.  Each frame is
// compressed exactly as compress_integer()/compress_numeric_u32() would
// compress it alone, so the frames are byte-identical to the per-block R path
// regardless of the thread count.
template <typename T>
SEXP compress_blocks(const std::vector<T>& values, std::size_t block_rows,
                     unsigned char dtype, SEXP level, SEXP page_n, SEXP threads) {
  const std::size_t n = values.size();
  const std::size_t blocks = n ? (n + block_rows - 1) / block_rows : 0;
  const CompressorPcoChunkConfig config{
    static_cast<unsigned int>(Rf_asInteger(level)),
    static_cast<std::size_t>(std::max(0, Rf_asInteger(page_n)))
  };
  std::vector<std::vector<unsigned char>> blobs(blocks);
  std::vector<std::size_t> capacity(blocks);
  for (std::size_t block = 0; block < blocks; ++block) {
    const std::size_t count = std::min(n, (block + 1) * block_rows) - block * block_rows;
    capacity[block] = compressor_pco_guarantee_file_size(count, dtype);
    if (!capacity[block]) compressor_fail("native Pcodec returned zero output capacity");
  }
  std::vector<CompressorPcoError> status(blocks, COMPRESSOR_PCO_SUCCESS);
  std::vector<int> overflow(blocks, 0);
  try {
    run_parallel_tasks(blocks, thread_count_arg(threads), [&](std::size_t block) {
      const std::size_t start = block * block_rows;
      const std::size_t count = std::min(n, start + block_rows) - start;
      std::vector<unsigned char>& blob = blobs[block];
      blob.resize(capacity[block]);
      std::size_t written = 0;
      status[block] = compressor_pco_compress_into(
        values.data() + start, count, dtype, &config, blob.data(), capacity[block], &written);
      if (written > capacity[block]) overflow[block] = 1;
      blob.resize(std::min(written, capacity[block]));
    });
  } catch (const std::exception& failure) {
    throw std::runtime_error(std::string("native Pcodec block compression failed: ") +
                             failure.what());
  }
  for (std::size_t block = 0; block < blocks; ++block) {
    check_status(status[block], dtype == kPcoTypeU32 ? "uint32 compression" : "compression");
    if (overflow[block]) compressor_fail("native Pcodec wrote beyond its output buffer");
  }
  return raw_list(blobs);
}

template <typename T>
std::vector<T> integer_stream_values(SEXP input) {
  if (TYPEOF(input) != INTSXP) compressor_fail("native Pcodec expects an integer vector");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  std::vector<T> values(n);
  const int* source = INTEGER(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (source[i] == NA_INTEGER || source[i] < 0 ||
        static_cast<std::uint64_t>(source[i]) >
          static_cast<std::uint64_t>(std::numeric_limits<T>::max())) {
      compressor_fail("native Pcodec input contains an out-of-range integer");
    }
    values[i] = static_cast<T>(source[i]);
  }
  return values;
}

std::vector<std::uint32_t> numeric_u32_stream_values(SEXP input) {
  if (TYPEOF(input) != REALSXP) compressor_fail("native Pcodec uint32 input must be numeric");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  std::vector<std::uint32_t> values(n);
  const double* source = REAL(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (!R_FINITE(source[i]) || source[i] < 0.0 ||
        source[i] != std::floor(source[i]) ||
        source[i] > 4294967295.0) {
      compressor_fail("native Pcodec uint32 input contains an invalid value");
    }
    values[i] = static_cast<std::uint32_t>(source[i]);
  }
  return values;
}

void append_le32(std::vector<unsigned char>& out, std::size_t at, std::uint32_t bits) {
  out[at] = static_cast<unsigned char>(bits & 0xffu);
  out[at + 1] = static_cast<unsigned char>((bits >> 8) & 0xffu);
  out[at + 2] = static_cast<unsigned char>((bits >> 16) & 0xffu);
  out[at + 3] = static_cast<unsigned char>((bits >> 24) & 0xffu);
}

std::uint32_t float_bits(double value) {
  // writeBin(size = 4) narrows with a C cast to float.
  const float narrowed = static_cast<float>(value);
  std::uint32_t bits;
  std::memcpy(&bits, &narrowed, sizeof(bits));
  return bits;
}

SEXP compress_numeric_u32(SEXP input, SEXP level, SEXP page_n) {
  if (TYPEOF(input) != REALSXP) compressor_fail("native Pcodec uint32 input must be numeric");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  if (n == 0) return Rf_allocVector(RAWSXP, 0);
  std::vector<std::uint32_t> values(n);
  const double* source = REAL(input);
  for (std::size_t i = 0; i < n; ++i) {
    if (!R_FINITE(source[i]) || source[i] < 0.0 ||
        source[i] != std::floor(source[i]) ||
        source[i] > 4294967295.0) {
      compressor_fail("native Pcodec uint32 input contains an invalid value");
    }
    values[i] = static_cast<std::uint32_t>(source[i]);
  }
  CompressorPcoChunkConfig config{
    static_cast<unsigned int>(Rf_asInteger(level)),
    static_cast<std::size_t>(std::max(0, Rf_asInteger(page_n)))
  };
  const std::size_t capacity = compressor_pco_guarantee_file_size(n, kPcoTypeU32);
  if (!capacity) compressor_fail("native Pcodec returned zero output capacity");
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(capacity)));
  std::size_t written = 0;
  check_status(compressor_pco_compress_into(
      values.data(), n, kPcoTypeU32, &config, RAW(output), capacity, &written),
    "uint32 compression");
  if (written > capacity) compressor_fail("native Pcodec wrote beyond its output buffer");
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
  if (TYPEOF(compressed) != RAWSXP) compressor_fail("compressed must be a raw vector");
  const R_xlen_t expected = Rf_asInteger(n);
  if (expected < 0) compressor_fail("n must be non-negative");
  if (expected == 0) return Rf_allocVector(REALSXP, 0);
  std::vector<std::uint32_t> values(static_cast<std::size_t>(expected));
  std::size_t written = 0;
  check_status(compressor_pco_decompress_into(
      RAW(compressed), static_cast<std::size_t>(XLENGTH(compressed)), kPcoTypeU32,
      values.data(), static_cast<std::size_t>(expected), &written),
    "uint32 decompression");
  if (written != static_cast<std::size_t>(expected)) {
    compressor_fail("native Pcodec returned %zu values; expected %td", written, expected);
  }
  SEXP output = PROTECT(Rf_allocVector(REALSXP, expected));
  double* target = REAL(output);
  for (std::size_t i = 0; i < written; ++i) target[i] = values[i];
  UNPROTECT(1);
  return output;
}

SEXP zstd_compress(SEXP input, SEXP level) {
  if (TYPEOF(input) != RAWSXP) compressor_fail("Zstandard input must be a raw vector");
  const std::size_t n = static_cast<std::size_t>(XLENGTH(input));
  if (n == 0) return Rf_allocVector(RAWSXP, 0);
  const int compression_level = Rf_asInteger(level);
  const std::size_t capacity = n + (n / 8) + 131072;
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(capacity)));
  std::size_t written = 0;
  check_status(compressor_zstd_compress_into(
      RAW(input), n, compression_level, RAW(output), capacity, &written),
    "Zstandard compression");
  if (written != capacity) {
    SEXP trimmed = PROTECT(Rf_allocVector(RAWSXP, static_cast<R_xlen_t>(written)));
    if (written) std::memcpy(RAW(trimmed), RAW(output), written);
    UNPROTECT(2);
    return trimmed;
  }
  UNPROTECT(1);
  return output;
}

SEXP zstd_decompress(SEXP input, SEXP expected) {
  if (TYPEOF(input) != RAWSXP) compressor_fail("Zstandard input must be a raw vector");
  const R_xlen_t expected_length = Rf_asInteger(expected);
  if (expected_length < 0) compressor_fail("expected decompressed length must be non-negative");
  if (expected_length == 0) return Rf_allocVector(RAWSXP, 0);
  SEXP output = PROTECT(Rf_allocVector(RAWSXP, expected_length));
  std::size_t written = 0;
  check_status(compressor_zstd_decompress_into(
      RAW(input), static_cast<std::size_t>(XLENGTH(input)), RAW(output),
      static_cast<std::size_t>(expected_length), &written),
    "Zstandard decompression");
  if (written != static_cast<std::size_t>(expected_length)) {
    UNPROTECT(1);
    compressor_fail("Zstandard returned %zu bytes; expected %td", written, expected_length);
  }
  UNPROTECT(1);
  return output;
}
}

extern "C" SEXP compressor_pcodec_native_available() {
  return Rf_ScalarLogical(TRUE);
}

extern "C" SEXP compressor_pcodec_compress_u8(SEXP input, SEXP level, SEXP page_n) {
  return compressor_guard([&]() -> SEXP {
    return compress_integer<std::uint8_t>(input, level, page_n, kPcoTypeU8);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_compress_u16(SEXP input, SEXP level, SEXP page_n) {
  return compressor_guard([&]() -> SEXP {
    return compress_integer<std::uint16_t>(input, level, page_n, kPcoTypeU16);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_compress_u32(SEXP input, SEXP level, SEXP page_n) {
  return compressor_guard([&]() -> SEXP {
    return compress_numeric_u32(input, level, page_n);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_decompress_u8(SEXP compressed, SEXP n) {
  return compressor_guard([&]() -> SEXP {
    return decompress_integer<std::uint8_t>(compressed, n, kPcoTypeU8);
  }, "native Pcodec error");
}

// Zero-based row IDs (row_start + i) of the nonzero values of one aligned
// binary u8 flag block; any value above 1 is an error.  Avoids returning the
// block's full integer vector to R only to scan it again there.
extern "C" SEXP compressor_pcodec_flag_rows_u8(SEXP compressed, SEXP n, SEXP row_start) {
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(compressed) != RAWSXP) compressor_fail("compressed must be a raw vector");
    const R_xlen_t expected = Rf_asInteger(n);
    if (expected < 0) compressor_fail("n must be non-negative");
    const double start = Rf_asReal(row_start);
    if (!(start >= 0) || start + static_cast<double>(expected) >
                             static_cast<double>(std::numeric_limits<int>::max())) {
      compressor_fail("flag block row range does not fit R's integer type");
    }
    if (expected == 0) return Rf_allocVector(INTSXP, 0);
    std::vector<std::uint8_t> values(static_cast<std::size_t>(expected));
    std::size_t written = 0;
    check_status(compressor_pco_decompress_into(
        RAW(compressed), static_cast<std::size_t>(XLENGTH(compressed)), kPcoTypeU8,
        values.data(), static_cast<std::size_t>(expected), &written),
      "decompression");
    if (written != static_cast<std::size_t>(expected)) {
      compressor_fail("native Pcodec returned %zu values; expected %td", written, expected);
    }
    R_xlen_t hits = 0;
    for (std::size_t i = 0; i < written; ++i) {
      if (values[i] > 1u) {
        compressor_fail("native Pcodec p-value flag payload is not a binary row-aligned stream");
      }
      hits += values[i];
    }
    SEXP output = PROTECT(Rf_allocVector(INTSXP, hits));
    int* target = INTEGER(output);
    const int base = static_cast<int>(start);
    R_xlen_t k = 0;
    for (std::size_t i = 0; i < written; ++i) {
      if (values[i]) target[k++] = base + static_cast<int>(i);
    }
    UNPROTECT(1);
    return output;
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_decompress_u16(SEXP compressed, SEXP n) {
  return compressor_guard([&]() -> SEXP {
    return decompress_integer<std::uint16_t>(compressed, n, kPcoTypeU16);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_decompress_u32(SEXP compressed, SEXP n) {
  return compressor_guard([&]() -> SEXP {
    return decompress_numeric_u32(compressed, n);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_zstd_compress(SEXP input, SEXP level) {
  return compressor_guard([&]() -> SEXP {
    return zstd_compress(input, level);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_zstd_decompress(SEXP input, SEXP expected) {
  return compressor_guard([&]() -> SEXP {
    return zstd_decompress(input, expected);
  }, "native Pcodec error");
}

extern "C" SEXP compressor_pcodec_compress_blocks(SEXP input, SEXP dtype, SEXP block_rows,
                                                  SEXP level, SEXP page_n, SEXP threads) {
  return compressor_guard([&]() -> SEXP {
    const int rows = Rf_asInteger(block_rows);
    if (rows == NA_INTEGER || rows < 1) compressor_fail("native Pcodec stream block_rows must be positive");
    const char* type = CHAR(STRING_ELT(dtype, 0));
    const std::size_t block = static_cast<std::size_t>(rows);
    if (std::strcmp(type, "u8") == 0) {
      return compress_blocks<std::uint8_t>(integer_stream_values<std::uint8_t>(input), block,
                                           kPcoTypeU8, level, page_n, threads);
    }
    if (std::strcmp(type, "u16") == 0) {
      return compress_blocks<std::uint16_t>(integer_stream_values<std::uint16_t>(input), block,
                                            kPcoTypeU16, level, page_n, threads);
    }
    if (std::strcmp(type, "u32") == 0) {
      return compress_blocks<std::uint32_t>(numeric_u32_stream_values(input), block,
                                            kPcoTypeU32, level, page_n, threads);
    }
    compressor_fail("unsupported native Pcodec dtype: %s", type);
    return R_NilValue;
  }, "native Pcodec error");
}

// Build and Zstandard-compress the per-value-block exception frames.  A frame
// holds the block's records column by column, exactly as
// pcodec_native_exception_bytes() writes them: int32 rows, float32 z,
// float32 log2se, float32 eaf, int8 flags, all little-endian.  `block_stops`
// are the exclusive 0-based row ends of the value blocks; `row` must be
// strictly increasing (checked by the caller).  Returns
// list(blobs, raw_lengths, counts).
extern "C" SEXP compressor_exception_blocks(SEXP row, SEXP z, SEXP log2se, SEXP eaf,
                                            SEXP flags, SEXP block_stops, SEXP level,
                                            SEXP threads) {
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(row) != INTSXP || TYPEOF(flags) != INTSXP || TYPEOF(z) != REALSXP ||
        TYPEOF(log2se) != REALSXP || TYPEOF(eaf) != REALSXP || TYPEOF(block_stops) != REALSXP) {
      compressor_fail("native exception frame inputs have the wrong types");
    }
    const std::size_t n = static_cast<std::size_t>(XLENGTH(row));
    if (static_cast<std::size_t>(XLENGTH(z)) != n || static_cast<std::size_t>(XLENGTH(log2se)) != n ||
        static_cast<std::size_t>(XLENGTH(eaf)) != n || static_cast<std::size_t>(XLENGTH(flags)) != n) {
      compressor_fail("native exception frame inputs have unequal lengths");
    }
    const std::size_t blocks = static_cast<std::size_t>(XLENGTH(block_stops));
    const int* rows = INTEGER(row);
    const double* stops = REAL(block_stops);
    // Member range [first[b], first[b + 1]) of each block, by a linear sweep.
    std::vector<std::size_t> first(blocks + 1, n);
    std::size_t at = 0;
    for (std::size_t block = 0; block < blocks; ++block) {
      first[block] = at;
      while (at < n && static_cast<double>(rows[at]) < stops[block]) ++at;
    }
    first[blocks] = at;
    if (at != n) compressor_fail("native Pcodec exception row is outside the value blocks");
    const double* zv = REAL(z);
    const double* sv = REAL(log2se);
    const double* ev = REAL(eaf);
    const int* fv = INTEGER(flags);
    const int compression_level = Rf_asInteger(level);
    std::vector<std::vector<unsigned char>> blobs(blocks);
    std::vector<int> raw_lengths(blocks, 0);
    std::vector<int> counts(blocks, 0);
    std::vector<CompressorPcoError> status(blocks, COMPRESSOR_PCO_SUCCESS);
    try {
      run_parallel_tasks(blocks, thread_count_arg(threads), [&](std::size_t block) {
        const std::size_t lo = first[block];
        const std::size_t count = first[block + 1] - lo;
        counts[block] = static_cast<int>(count);
        if (!count) return;
        std::vector<unsigned char> frame(count * 17u);
        for (std::size_t i = 0; i < count; ++i) {
          append_le32(frame, 4u * i, static_cast<std::uint32_t>(rows[lo + i]));
          append_le32(frame, 4u * (count + i), float_bits(zv[lo + i]));
          append_le32(frame, 4u * (2u * count + i), float_bits(sv[lo + i]));
          append_le32(frame, 4u * (3u * count + i), float_bits(ev[lo + i]));
          frame[16u * count + i] = static_cast<unsigned char>(static_cast<signed char>(fv[lo + i]));
        }
        raw_lengths[block] = static_cast<int>(frame.size());
        const std::size_t capacity = frame.size() + (frame.size() / 8) + 131072;
        std::vector<unsigned char>& blob = blobs[block];
        blob.resize(capacity);
        std::size_t written = 0;
        status[block] = compressor_zstd_compress_into(
          frame.data(), frame.size(), compression_level, blob.data(), capacity, &written);
        blob.resize(std::min(written, capacity));
      });
    } catch (const std::exception& failure) {
      throw std::runtime_error(std::string("native exception frame compression failed: ") +
                               failure.what());
    }
    for (std::size_t block = 0; block < blocks; ++block) {
      check_status(status[block], "Zstandard compression");
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, raw_list(blobs));
    SEXP lengths = Rf_allocVector(INTSXP, static_cast<R_xlen_t>(blocks));
    SET_VECTOR_ELT(out, 1, lengths);
    SEXP member_counts = Rf_allocVector(INTSXP, static_cast<R_xlen_t>(blocks));
    SET_VECTOR_ELT(out, 2, member_counts);
    for (std::size_t block = 0; block < blocks; ++block) {
      INTEGER(lengths)[block] = raw_lengths[block];
      INTEGER(member_counts)[block] = counts[block];
    }
    UNPROTECT(1);
    return out;
  }, "native Pcodec error");
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
extern "C" SEXP compressor_pcodec_flag_rows_u8(SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_decompress_u16(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_decompress_u32(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_zstd_compress(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_zstd_decompress(SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_pcodec_compress_blocks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}
extern "C" SEXP compressor_exception_blocks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) {
  Rf_error("native Pcodec is not available in this build");
}

#endif
