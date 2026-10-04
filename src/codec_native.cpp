#include <algorithm>
#include <atomic>
#include <climits>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <R.h>
#include <Rinternals.h>

#include "native_error.h"

// Decoded values must not depend on the compiler fusing a * b + c into one
// rounding (FMA): every reader and platform must produce the same bits.
// Makevars also passes -ffp-contract=off where the compiler accepts it.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

#ifdef COMPRESSOR_NATIVE_PCODEC
#include "pcodec_native.h"
#endif

// Rinternals.h defines `length` as a macro, which collides with libc++'s
// locale headers when they are included after it. This translation unit uses
// XLENGTH/Rf_length explicitly, so the macro is not needed below.
#ifdef length
#undef length
#endif

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInvSqrtTwo = 0.70710678118654752440;
constexpr double kPi = 3.14159265358979323846;

double clamp_probability_eaf(double value) {
  if (!R_FINITE(value)) return 0.5;
  return std::min(1.0 - 1e-12, std::max(1e-12, value));
}

SEXP scalar_name(const char* value) {
  return Rf_mkCharCE(value, CE_UTF8);
}

bool requested_has(SEXP requested, const char* target) {
  const R_xlen_t n = XLENGTH(requested);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP value = STRING_ELT(requested, i);
    if (value == NA_STRING) {
      throw std::runtime_error("native Pcodec requested streams contain NA");
    }
    if (std::string(CHAR(value)) == target) return true;
  }
  return false;
}

double native_scalar_integer(SEXP value, bool allow_zero, const char* label) {
  if ((TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP) ||
      XLENGTH(value) != 1) {
    throw std::runtime_error(std::string(label) + " must be one integer");
  }
  const double number = Rf_asReal(value);
  if (!R_FINITE(number) || number != std::floor(number) ||
      (allow_zero ? number < 0.0 : number < 1.0)) {
    throw std::runtime_error(std::string(label) +
                             (allow_zero ? " must be non-negative" :
                                           " must be positive"));
  }
  return number;
}

#ifdef COMPRESSOR_NATIVE_PCODEC

constexpr unsigned char kPcoTypeU8 = 10;
constexpr unsigned char kPcoTypeU16 = 7;
constexpr unsigned char kPcoTypeU32 = 1;

struct NativePcodecBlock {
  std::uint64_t offset;
  std::uint64_t length;
  std::uint64_t values;
  std::uint64_t row_start;
  std::uint64_t row_stop;
  std::uint64_t first_position;
};

struct NativePcodecExceptionBlock {
  std::uint64_t offset;
  std::uint64_t length;
  std::uint64_t count;
  std::uint64_t raw_length;
};

std::uint64_t matrix_uint64(SEXP matrix, R_xlen_t row, int column,
                            const char* label) {
  if (TYPEOF(matrix) != REALSXP || !Rf_isMatrix(matrix)) {
    throw std::runtime_error(std::string(label) + " must be a numeric matrix");
  }
  const SEXP dimensions = Rf_getAttrib(matrix, R_DimSymbol);
  if (XLENGTH(dimensions) != 2 || INTEGER(dimensions)[1] <= column) {
    throw std::runtime_error(std::string(label) + " has too few columns");
  }
  const R_xlen_t rows = INTEGER(dimensions)[0];
  if (row < 0 || row >= rows) throw std::runtime_error("native Pcodec block row is invalid");
  const double value = REAL(matrix)[row + rows * column];
  if (!R_FINITE(value) || value < 0.0 || value != std::floor(value) ||
      value > static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
    throw std::runtime_error(std::string(label) + " contains an invalid integer");
  }
  return static_cast<std::uint64_t>(value);
}

std::vector<NativePcodecBlock> read_native_blocks(SEXP matrix,
                                                  R_xlen_t n,
                                                  bool has_first_position) {
  if (TYPEOF(matrix) != REALSXP || !Rf_isMatrix(matrix)) {
    throw std::runtime_error("native Pcodec blocks must be a numeric matrix");
  }
  const SEXP dimensions = Rf_getAttrib(matrix, R_DimSymbol);
  if (XLENGTH(dimensions) != 2 || INTEGER(dimensions)[1] < (has_first_position ? 6 : 5)) {
    throw std::runtime_error("native Pcodec block matrix has the wrong shape");
  }
  const R_xlen_t rows = INTEGER(dimensions)[0];
  std::vector<NativePcodecBlock> blocks;
  blocks.reserve(static_cast<std::size_t>(rows));
  std::uint64_t expected_row = 0;
  for (R_xlen_t row = 0; row < rows; ++row) {
    NativePcodecBlock block{
      matrix_uint64(matrix, row, 0, "native Pcodec block offset"),
      matrix_uint64(matrix, row, 1, "native Pcodec block length"),
      matrix_uint64(matrix, row, 2, "native Pcodec block value count"),
      matrix_uint64(matrix, row, 3, "native Pcodec block row start"),
      matrix_uint64(matrix, row, 4, "native Pcodec block row stop"),
      has_first_position ? matrix_uint64(matrix, row, 5,
                                         "native Pcodec block first position") : 0
    };
    if (block.row_start != expected_row || block.row_stop < block.row_start ||
        block.row_stop - block.row_start != block.values || block.row_stop > n) {
      throw std::runtime_error("native Pcodec blocks do not form a contiguous row index");
    }
    if (has_first_position &&
        block.first_position > std::numeric_limits<std::uint32_t>::max()) {
      throw std::runtime_error("native Pcodec block first position exceeds uint32");
    }
    expected_row = block.row_stop;
    blocks.push_back(block);
  }
  if (expected_row != n) {
    throw std::runtime_error("native Pcodec blocks do not cover the declared row count");
  }
  return blocks;
}

void require_block_row_parity(const std::vector<NativePcodecBlock>& left,
                              const std::vector<NativePcodecBlock>& right,
                              const char* label) {
  if (left.size() != right.size()) {
    throw std::runtime_error(std::string(label) + " block counts differ");
  }
  for (std::size_t i = 0; i < left.size(); ++i) {
    if (left[i].row_start != right[i].row_start ||
        left[i].row_stop != right[i].row_stop ||
        left[i].values != right[i].values) {
      throw std::runtime_error(std::string(label) + " blocks are not row-aligned");
    }
  }
}

std::vector<NativePcodecExceptionBlock> read_native_exception_blocks(SEXP matrix) {
  if (TYPEOF(matrix) != REALSXP || !Rf_isMatrix(matrix)) {
    throw std::runtime_error("native Pcodec exception blocks must be a numeric matrix");
  }
  const SEXP dimensions = Rf_getAttrib(matrix, R_DimSymbol);
  if (XLENGTH(dimensions) != 2 || INTEGER(dimensions)[1] < 4) {
    throw std::runtime_error("native Pcodec exception block matrix has the wrong shape");
  }
  const R_xlen_t rows = INTEGER(dimensions)[0];
  std::vector<NativePcodecExceptionBlock> blocks;
  blocks.reserve(static_cast<std::size_t>(rows));
  for (R_xlen_t row = 0; row < rows; ++row) {
    blocks.push_back({
      matrix_uint64(matrix, row, 0, "native Pcodec exception offset"),
      matrix_uint64(matrix, row, 1, "native Pcodec exception length"),
      matrix_uint64(matrix, row, 2, "native Pcodec exception count"),
      matrix_uint64(matrix, row, 3, "native Pcodec exception raw length")
    });
  }
  return blocks;
}

std::string native_path(SEXP files, int index, const char* label) {
  if (TYPEOF(files) != STRSXP || XLENGTH(files) <= index ||
      STRING_ELT(files, index) == NA_STRING) {
    throw std::runtime_error(std::string("native Pcodec ") + label + " path is invalid");
  }
  return CHAR(STRING_ELT(files, index));
}

class NativePcodecFile {
 public:
  explicit NativePcodecFile(std::string path) : path_(std::move(path)) {
    stream_.open(path_, std::ios::binary);
    if (!stream_) throw std::runtime_error("cannot open native Pcodec stream: " + path_);
    stream_.seekg(0, std::ios::end);
    const std::streamoff file_size = stream_.tellg();
    if (file_size < 0) throw std::runtime_error("cannot determine native Pcodec stream size: " + path_);
    file_size_ = static_cast<std::uint64_t>(file_size);
  }

  std::vector<std::uint8_t> read(std::uint64_t offset, std::uint64_t length) {
    if (length > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
        offset > file_size_ || length > file_size_ - offset) {
      throw std::runtime_error("native Pcodec stream block is outside the file: " + path_);
    }
    std::vector<std::uint8_t> blob(static_cast<std::size_t>(length));
    stream_.clear();
    stream_.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    if (length > 0 && !stream_.read(reinterpret_cast<char*>(blob.data()),
                                    static_cast<std::streamsize>(length))) {
      throw std::runtime_error("native Pcodec stream block is truncated: " + path_);
    }
    return blob;
  }

  const std::string& path() const { return path_; }

 private:
  std::string path_;
  std::ifstream stream_;
  std::uint64_t file_size_ = 0;
};

template <typename T>
std::vector<T> native_decompress_block(NativePcodecFile& file,
                                       const NativePcodecBlock& block,
                                       unsigned char dtype) {
  if (block.values > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error("native Pcodec block has too many values");
  }
  std::vector<std::uint8_t> blob = file.read(block.offset, block.length);
  std::vector<T> values(static_cast<std::size_t>(block.values));
  std::size_t written = 0;
  const CompressorPcoError status = compressor_pco_decompress_into(
    blob.data(), blob.size(), dtype, values.data(), values.size(), &written);
  if (status != COMPRESSOR_PCO_SUCCESS || written != values.size()) {
    throw std::runtime_error("native Pcodec block decompression failed: " + file.path());
  }
  return values;
}

// Pcodec blocks are independent.  Give each worker its own ifstream rather
// than sharing a seekable stream between threads; this keeps the reader's
// hot path free of locks and also makes the function safe when a store is
// read concurrently by more than one R call.
template <typename T>
void native_decompress_blocks_parallel(
    const std::string& path,
    const std::vector<NativePcodecBlock>& blocks,
    unsigned char dtype,
    std::size_t requested_threads,
    std::vector<T>& output) {
  if (blocks.empty()) return;
  const std::size_t worker_count = std::max<std::size_t>(
    1, std::min<std::size_t>(requested_threads, blocks.size()));
  if (output.size() < static_cast<std::size_t>(blocks.back().row_stop)) {
    throw std::runtime_error("native Pcodec output vector is too small");
  }
  if (worker_count == 1) {
    NativePcodecFile file(path);
    for (const NativePcodecBlock& block : blocks) {
      std::vector<T> values = native_decompress_block<T>(file, block, dtype);
      std::copy(values.begin(), values.end(),
                output.begin() + static_cast<std::ptrdiff_t>(block.row_start));
    }
    return;
  }

  std::atomic<std::size_t> next_block(0);
  std::mutex error_mutex;
  std::exception_ptr first_error;
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    workers.emplace_back([&]() {
      try {
        NativePcodecFile file(path);
        while (true) {
          const std::size_t block_id = next_block.fetch_add(
            1, std::memory_order_relaxed);
          if (block_id >= blocks.size()) break;
          const NativePcodecBlock& block = blocks[block_id];
          std::vector<T> values = native_decompress_block<T>(file, block, dtype);
          std::copy(values.begin(), values.end(),
                    output.begin() + static_cast<std::ptrdiff_t>(block.row_start));
        }
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!first_error) first_error = std::current_exception();
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  if (first_error) std::rethrow_exception(first_error);
}

std::uint32_t native_read_u32_le(const std::uint8_t* bytes) {
  return static_cast<std::uint32_t>(bytes[0]) |
    (static_cast<std::uint32_t>(bytes[1]) << 8) |
    (static_cast<std::uint32_t>(bytes[2]) << 16) |
    (static_cast<std::uint32_t>(bytes[3]) << 24);
}

float native_read_float_le(const std::uint8_t* bytes) {
  const std::uint32_t bits = native_read_u32_le(bytes);
  float value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void native_append_exception_block(
    NativePcodecFile& file, const NativePcodecExceptionBlock& block,
    const NativePcodecBlock& value_block, std::uint64_t row_count,
    std::int64_t* previous_row,
    const std::string& codec, std::vector<int>& rows, std::vector<double>& z,
    std::vector<double>& log2se, std::vector<double>& eaf,
    std::vector<int>& flags) {
  if (block.count == 0) {
    if (block.raw_length != 0 || block.length != 0) {
      throw std::runtime_error("native Pcodec empty exception block has a payload");
    }
    return;
  }
  if (block.count > value_block.values) {
    throw std::runtime_error("native Pcodec exception block has too many records");
  }
  if (block.count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() / 17)) {
    throw std::runtime_error("native Pcodec exception block is too large");
  }
  const std::uint64_t expected_raw = block.count * 17;
  if (block.raw_length != expected_raw) {
    throw std::runtime_error("native Pcodec exception block has an invalid raw length");
  }
  std::vector<std::uint8_t> compressed = file.read(block.offset, block.length);
  std::vector<std::uint8_t> raw(static_cast<std::size_t>(block.raw_length));
  if (codec == "zstd") {
    std::size_t written = 0;
    const CompressorPcoError status = compressor_zstd_decompress_into(
      compressed.data(), compressed.size(), raw.data(), raw.size(), &written);
    if (status != COMPRESSOR_PCO_SUCCESS || written != raw.size()) {
      throw std::runtime_error("native Pcodec exception decompression failed");
    }
  } else if (codec == "raw") {
    if (compressed.size() != raw.size()) {
      throw std::runtime_error("native Pcodec raw exception block is truncated");
    }
    raw.swap(compressed);
  } else {
    throw std::runtime_error("unsupported native Pcodec exception codec: " + codec);
  }
  for (std::uint64_t i = 0; i < block.count; ++i) {
    const std::size_t index = static_cast<std::size_t>(i);
    const std::size_t count = static_cast<std::size_t>(block.count);
    const std::uint32_t row = native_read_u32_le(raw.data() + index * 4);
    if (static_cast<std::uint64_t>(row) < value_block.row_start ||
        static_cast<std::uint64_t>(row) >= value_block.row_stop ||
        static_cast<std::uint64_t>(row) >= row_count ||
        (previous_row != nullptr && static_cast<std::int64_t>(row) <= *previous_row) ||
        row > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("native Pcodec exception row does not fit R integer");
    }
    if (previous_row != nullptr) *previous_row = static_cast<std::int64_t>(row);
    rows.push_back(static_cast<int>(row));
    z.push_back(static_cast<double>(native_read_float_le(raw.data() + count * 4 + index * 4)));
    log2se.push_back(static_cast<double>(native_read_float_le(raw.data() + count * 8 + index * 4)));
    eaf.push_back(static_cast<double>(native_read_float_le(raw.data() + count * 12 + index * 4)));
    const int record_flags = static_cast<int>(raw[count * 16 + index]);
    if (record_flags <= 0 || (record_flags & ~7) != 0) {
      throw std::runtime_error("native Pcodec exception flags are invalid");
    }
    flags.push_back(record_flags);
  }
}

#endif

}  // namespace

extern "C" SEXP compressor_decode_native(
    SEXP z_code,
    SEXP se_code,
    SEXP eaf_code,
    SEXP z_min,
    SEXP z_max,
    SEXP z_count,
    SEXP se_count,
    SEXP eaf_count,
    SEXP z_bits,
    SEXP se_bits,
    SEXP eaf_bits,
    SEXP block_rows,
    SEXP centres,
    SEXP exception_row,
    SEXP exception_z,
    SEXP exception_se,
    SEXP exception_eaf,
    SEXP exception_flags,
    SEXP include_beta,
    SEXP include_p,
    SEXP se_residual_min,
    SEXP se_residual_max,
    SEXP centre_index) {
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(z_code) != INTSXP || TYPEOF(se_code) != INTSXP ||
        TYPEOF(eaf_code) != INTSXP) {
      compressor_fail("native decoder requires integer code vectors");
    }
    const R_xlen_t n = XLENGTH(z_code);
    if (XLENGTH(se_code) != n || XLENGTH(eaf_code) != n) {
      compressor_fail("native decoder code vectors must have equal lengths");
    }
  
    const int z_bits_value = Rf_asInteger(z_bits);
    const int se_bits_value = Rf_asInteger(se_bits);
    const int eaf_bits_value = Rf_asInteger(eaf_bits);
    const int z_count_value = Rf_asInteger(z_count);
    const int se_count_value = Rf_asInteger(se_count);
    const int eaf_count_value = Rf_asInteger(eaf_count);
    const double se_residual_min_value = Rf_asReal(se_residual_min);
    const double se_residual_max_value = Rf_asReal(se_residual_max);
    const double block_rows_number = Rf_asReal(block_rows);
    if (!R_FINITE(block_rows_number) || block_rows_number < 1.0 ||
        block_rows_number != std::floor(block_rows_number) ||
        block_rows_number > static_cast<double>(std::numeric_limits<int>::max())) {
      compressor_fail("native decoder block_rows must be a positive integer");
    }
    const int block_rows_value = static_cast<int>(block_rows_number);
    if (z_bits_value <= 0 || se_bits_value <= 0 || eaf_bits_value <= 0 ||
        z_count_value <= 0 || se_count_value <= 0 || eaf_count_value <= 0) {
      compressor_fail("native decoder metadata contains invalid code domains");
    }
    if (!R_FINITE(se_residual_min_value) || !R_FINITE(se_residual_max_value) ||
        se_residual_max_value <= se_residual_min_value) {
      compressor_fail("native decoder metadata contains an invalid SE residual range");
    }
    const double z_min_value = Rf_asReal(z_min);
    const double z_max_value = Rf_asReal(z_max);
    if (!R_FINITE(z_min_value) || !R_FINITE(z_max_value) ||
        z_max_value <= z_min_value) {
      compressor_fail("native decoder metadata contains an invalid Z range");
    }
    if (z_bits_value > 16 || se_bits_value > 16 || eaf_bits_value > 16) {
      compressor_fail("native decoder code domains exceed the supported 16-bit limit");
    }
  
    const std::size_t z_width = static_cast<std::size_t>(1) << z_bits_value;
    const std::size_t se_width = static_cast<std::size_t>(1) << se_bits_value;
    const std::size_t eaf_width = static_cast<std::size_t>(1) << eaf_bits_value;
    if (z_width > 65536 || se_width > 65536 || eaf_width > 65536) {
      compressor_fail("native decoder code domains are too large");
    }
    if (static_cast<std::size_t>(z_count_value) + 2 > z_width ||
        static_cast<std::size_t>(se_count_value) + 2 > se_width ||
        static_cast<std::size_t>(eaf_count_value) + 1 > eaf_width) {
      compressor_fail("native decoder metadata does not fit its code domains");
    }
    if (static_cast<std::size_t>(se_count_value) * eaf_width > 16u * 1024u * 1024u) {
      compressor_fail("native decoder SE lookup table is too large");
    }
  
    std::vector<double> z_table(z_width, kNaN);
    const double z_step = (z_max_value - z_min_value) /
                          static_cast<double>(z_count_value);
    for (int code = 0; code < z_count_value &&
                         static_cast<std::size_t>(code) < z_width; ++code) {
      z_table[static_cast<std::size_t>(code)] = z_min_value +
        (static_cast<double>(code) + 0.5) * z_step;
    }
    std::vector<double> p_table(z_width, kNaN);
    for (int code = 0; code < z_count_value &&
                         static_cast<std::size_t>(code) < z_width; ++code) {
      p_table[static_cast<std::size_t>(code)] =
        std::erfc(std::abs(z_table[static_cast<std::size_t>(code)]) * kInvSqrtTwo);
    }
  
    std::vector<double> eaf_table(eaf_width, kNaN);
    // q_encode uses codes 0..eaf_count inclusive and reserves the final code
    // (2^bits-1) for an exact/missing exception.
    for (int code = 0; code <= eaf_count_value &&
                         static_cast<std::size_t>(code) < eaf_width; ++code) {
      eaf_table[static_cast<std::size_t>(code)] =
        std::pow(std::sin((static_cast<double>(code) /
                           static_cast<double>(eaf_count_value)) * kPi / 2.0), 2.0);
    }
    const double fallback_eaf = 0.5;
    std::vector<double> se_table(static_cast<std::size_t>(se_count_value) * eaf_width, kNaN);
    for (int se_value = 0; se_value < se_count_value; ++se_value) {
      const double residual = se_residual_min_value +
        (static_cast<double>(se_value) + 0.5) *
          ((se_residual_max_value - se_residual_min_value) / se_count_value);
      for (std::size_t eaf_value = 0; eaf_value < eaf_width; ++eaf_value) {
        const double decoded_eaf = R_FINITE(eaf_table[eaf_value])
          ? eaf_table[eaf_value] : fallback_eaf;
        const double safe_eaf = clamp_probability_eaf(decoded_eaf);
        const double correction = -0.5 *
          std::log2(2.0 * safe_eaf * (1.0 - safe_eaf));
        se_table[static_cast<std::size_t>(se_value) * eaf_width + eaf_value] =
          std::exp2(residual + correction);
      }
    }
  
    const R_xlen_t centre_count = XLENGTH(centres);
    const double* centre_values = TYPEOF(centres) == REALSXP ? REAL(centres) : nullptr;
    if (centre_count > 0 && centre_values == nullptr) {
      compressor_fail("native decoder block centres must be numeric");
    }
    // The block centre is constant for every row in a block.  Materialise its
    // multiplicative factor once; calling exp2() in the row loop made the
    // R-facing decoder materially slower than the historical native reader.
    std::vector<double> centre_factors(static_cast<std::size_t>(centre_count));
    for (R_xlen_t i = 0; i < centre_count; ++i) {
      if (!R_FINITE(centre_values[i])) {
        compressor_fail("native decoder block centres contain a non-finite value");
      }
      centre_factors[static_cast<std::size_t>(i)] = std::exp2(centre_values[i]);
    }
  
    // Optional zero-based SE-centre block of every row, for rows that are not
    // a contiguous run from row 0 (selective, region and candidate reads).
    // Without it the centre is row / block_rows, as for a full read. Either
    // way the arithmetic is the same, so every reader decodes a row to the
    // same bits.
    const int* centre_ids = nullptr;
    if (centre_index != R_NilValue) {
      if (TYPEOF(centre_index) != INTSXP || XLENGTH(centre_index) != n) {
        compressor_fail("native decoder centre index must be one integer per row");
      }
      centre_ids = INTEGER(centre_index);
    }
    const bool want_beta = Rf_asLogical(include_beta) == TRUE;
    const bool want_p = Rf_asLogical(include_p) == TRUE;
    const int output_count = 3 + (want_beta ? 1 : 0) + (want_p ? 1 : 0);
    int protect_count = 0;
    SEXP output = PROTECT(Rf_allocVector(VECSXP, output_count));
    ++protect_count;
    SEXP z = PROTECT(Rf_allocVector(REALSXP, n));
    ++protect_count;
    SEXP se = PROTECT(Rf_allocVector(REALSXP, n));
    ++protect_count;
    SEXP eaf = PROTECT(Rf_allocVector(REALSXP, n));
    ++protect_count;
    SEXP beta = R_NilValue;
    SEXP p = R_NilValue;
    if (want_beta) {
      beta = PROTECT(Rf_allocVector(REALSXP, n));
      ++protect_count;
    }
    if (want_p) {
      p = PROTECT(Rf_allocVector(REALSXP, n));
      ++protect_count;
    }
  
    double* z_out = REAL(z);
    double* se_out = REAL(se);
    double* eaf_out = REAL(eaf);
    double* beta_out = want_beta ? REAL(beta) : nullptr;
    double* p_out = want_p ? REAL(p) : nullptr;
    const int* z_codes = INTEGER(z_code);
    const int* se_codes = INTEGER(se_code);
    const int* eaf_codes = INTEGER(eaf_code);
    for (R_xlen_t row = 0; row < n; ++row) {
      const int z_value = z_codes[row];
      const int se_value = se_codes[row];
      const int eaf_value = eaf_codes[row];
      const bool z_ok = z_value >= 0 && static_cast<std::size_t>(z_value) < z_width &&
                        z_value < z_count_value;
      const bool se_ok = se_value >= 0 && static_cast<std::size_t>(se_value) < se_width &&
                         se_value < se_count_value;
      const bool eaf_ok = eaf_value >= 0 && static_cast<std::size_t>(eaf_value) < eaf_width &&
                          R_FINITE(eaf_table[static_cast<std::size_t>(eaf_value)]);
      const std::size_t safe_eaf_value =
        eaf_value >= 0 && static_cast<std::size_t>(eaf_value) < eaf_width
          ? static_cast<std::size_t>(eaf_value) : eaf_width - 1;
      z_out[row] = z_ok ? z_table[static_cast<std::size_t>(z_value)] : kNaN;
      eaf_out[row] = eaf_ok ? eaf_table[static_cast<std::size_t>(eaf_value)] : kNaN;
      if (se_ok && centre_count > 0) {
        R_xlen_t block = centre_ids ? static_cast<R_xlen_t>(centre_ids[row])
                                    : row / block_rows_value;
        if (block < 0) block = 0;  // includes NA_INTEGER
        const R_xlen_t centre = std::min<R_xlen_t>(block, centre_count - 1);
        se_out[row] = se_table[static_cast<std::size_t>(se_value) * eaf_width +
                                safe_eaf_value] *
                      centre_factors[static_cast<std::size_t>(centre)];
      } else {
        se_out[row] = kNaN;
      }
      if (want_beta) beta_out[row] = z_out[row] * se_out[row];
      if (want_p) p_out[row] = z_ok
        ? p_table[static_cast<std::size_t>(z_value)] : kNaN;
    }
  
    const R_xlen_t exception_count = XLENGTH(exception_row);
    if (TYPEOF(exception_row) != INTSXP || TYPEOF(exception_z) != REALSXP ||
        TYPEOF(exception_se) != REALSXP || TYPEOF(exception_eaf) != REALSXP ||
        TYPEOF(exception_flags) != INTSXP ||
        XLENGTH(exception_z) != exception_count ||
        XLENGTH(exception_se) != exception_count ||
        XLENGTH(exception_eaf) != exception_count ||
        XLENGTH(exception_flags) != exception_count) {
      compressor_fail("native decoder exception vectors are malformed");
    }
    const int* exception_rows = INTEGER(exception_row);
    const double* exception_z_values = REAL(exception_z);
    const double* exception_se_values = REAL(exception_se);
    const double* exception_eaf_values = REAL(exception_eaf);
    const int* exception_flag_values = INTEGER(exception_flags);
    for (R_xlen_t i = 0; i < exception_count; ++i) {
      const int row = exception_rows[i];
      const int flags = exception_flag_values[i];
      if (flags <= 0 || (flags & ~7) != 0) {
        compressor_fail("native decoder exception flags are invalid");
      }
      if (row < 0 || static_cast<R_xlen_t>(row) >= n) continue;
      if (flags & 1) z_out[row] = exception_z_values[i];
      if (flags & 2) se_out[row] = exception_se_values[i];
      if (flags & 4) eaf_out[row] = exception_eaf_values[i];
      if (want_beta) beta_out[row] = z_out[row] * se_out[row];
      if (want_p && (flags & 1)) {
        p_out[row] = std::erfc(std::abs(z_out[row]) * kInvSqrtTwo);
      }
    }
  
    int index = 0;
    SET_VECTOR_ELT(output, index++, z);
    if (want_beta) SET_VECTOR_ELT(output, index++, beta);
    SET_VECTOR_ELT(output, index++, se);
    SET_VECTOR_ELT(output, index++, eaf);
    if (want_p) SET_VECTOR_ELT(output, index++, p);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, output_count));
    ++protect_count;
    index = 0;
    SET_STRING_ELT(names, index++, scalar_name("z"));
    if (want_beta) SET_STRING_ELT(names, index++, scalar_name("beta"));
    SET_STRING_ELT(names, index++, scalar_name("standard_error"));
    SET_STRING_ELT(names, index++, scalar_name("effect_allele_frequency"));
    if (want_p) SET_STRING_ELT(names, index++, scalar_name("p_value"));
    Rf_setAttrib(output, R_NamesSymbol, names);
    UNPROTECT(protect_count);
    return output;
  }, "native Pcodec error");
}

// Read the current on-disk native Pcodec streams directly. A former bridge
// below reads an intermediate representation; this reader instead consumes the
// `.cpr` stream files and block index used by current stores. It returns compact
// integer/numeric vectors so R never has to perform one readBin()/decompress/
// unlist cycle per block.
extern "C" SEXP compressor_read_pcodec_native_codes(
    SEXP files,
    SEXP position_blocks,
    SEXP substitution_blocks,
    SEXP z_blocks,
    SEXP eaf_blocks,
    SEXP se_blocks,
    SEXP exception_blocks,
    SEXP row_count,
    SEXP streams,
    SEXP exception_codec,
    SEXP threads) {
#ifdef COMPRESSOR_NATIVE_PCODEC
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(files) != STRSXP || XLENGTH(files) != 6 ||
        TYPEOF(streams) != STRSXP || TYPEOF(exception_codec) != STRSXP ||
        XLENGTH(exception_codec) != 1 || STRING_ELT(exception_codec, 0) == NA_STRING) {
      throw std::runtime_error("malformed arguments to native Pcodec stream reader");
    }
    const double row_count_value = native_scalar_integer(
      row_count, true, "native Pcodec row count");
    if (row_count_value > static_cast<double>(std::numeric_limits<R_xlen_t>::max())) {
      throw std::runtime_error("native Pcodec row count is invalid");
    }
    const R_xlen_t n = static_cast<R_xlen_t>(row_count_value);
    const double thread_number = native_scalar_integer(
      threads, false, "native Pcodec thread count");
    if (thread_number >= static_cast<double>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error("native Pcodec thread count is too large");
    }
    const std::size_t requested_threads = static_cast<std::size_t>(thread_number);
    const bool need_position = requested_has(streams, "position");
    const bool need_substitution = requested_has(streams, "substitution");
    const bool need_z = requested_has(streams, "z");
    const bool need_se = requested_has(streams, "se");
    const bool need_eaf = need_se || requested_has(streams, "eaf");
    const bool need_numeric = need_z || need_se || need_eaf;
    // "exceptions" alone reads just the exception sidecar (one file handle,
    // one call), e.g. for threshold candidate selection.
    const bool need_exceptions = need_numeric || requested_has(streams, "exceptions");

    std::vector<NativePcodecBlock> positions = read_native_blocks(
      position_blocks, static_cast<R_xlen_t>(n), true);
    std::vector<NativePcodecBlock> substitutions = read_native_blocks(
      substitution_blocks, static_cast<R_xlen_t>(n), true);
    std::vector<NativePcodecBlock> z_values = read_native_blocks(
      z_blocks, static_cast<R_xlen_t>(n), false);
    std::vector<NativePcodecBlock> eaf_values = read_native_blocks(
      eaf_blocks, static_cast<R_xlen_t>(n), false);
    std::vector<NativePcodecBlock> se_values = read_native_blocks(
      se_blocks, static_cast<R_xlen_t>(n), false);
    std::vector<NativePcodecExceptionBlock> exceptions =
      read_native_exception_blocks(exception_blocks);
    require_block_row_parity(positions, substitutions,
                             "native key stream");
    require_block_row_parity(z_values, eaf_values,
                             "native value stream");
    require_block_row_parity(z_values, se_values,
                             "native value stream");
    if (need_exceptions && exceptions.size() != z_values.size()) {
      throw std::runtime_error("native Pcodec exception index does not match value blocks");
    }

    int protect_count = 0;
    SEXP position = R_NilValue;
    SEXP substitution = R_NilValue;
    SEXP z = R_NilValue;
    SEXP se = R_NilValue;
    SEXP eaf = R_NilValue;
    if (need_position) {
      position = PROTECT(Rf_allocVector(REALSXP, n));
      ++protect_count;
    }
    if (need_substitution) {
      substitution = PROTECT(Rf_allocVector(INTSXP, n));
      ++protect_count;
    }
    if (need_z) {
      z = PROTECT(Rf_allocVector(INTSXP, n));
      ++protect_count;
    }
    if (need_se) {
      se = PROTECT(Rf_allocVector(INTSXP, n));
      ++protect_count;
    }
    if (need_eaf) {
      eaf = PROTECT(Rf_allocVector(INTSXP, n));
      ++protect_count;
    }

    std::unique_ptr<NativePcodecFile> exception_file;
    if (need_exceptions) exception_file.reset(new NativePcodecFile(
      native_path(files, 5, "exception")));

    if (need_position) {
      std::vector<std::uint32_t> gaps(static_cast<std::size_t>(n));
      native_decompress_blocks_parallel<std::uint32_t>(
        native_path(files, 0, "position"), positions, kPcoTypeU32,
        requested_threads, gaps);
      for (const NativePcodecBlock& block : positions) {
        if (block.values == 0 || gaps[static_cast<std::size_t>(block.row_start)] != 0) {
          throw std::runtime_error("native Pcodec position block does not reset its delta");
        }
        std::uint64_t current = block.first_position;
        for (std::uint64_t i = 0; i < block.values; ++i) {
          const std::uint64_t gap =
            gaps[static_cast<std::size_t>(block.row_start + i)];
          if (current > std::numeric_limits<std::uint32_t>::max() - gap) {
            throw std::runtime_error("native Pcodec position exceeds uint32");
          }
          current += gap;
          REAL(position)[block.row_start + i] = static_cast<double>(current);
        }
      }
    }
    if (need_substitution) {
      std::vector<std::uint8_t> codes(static_cast<std::size_t>(n));
      native_decompress_blocks_parallel<std::uint8_t>(
        native_path(files, 1, "substitution"), substitutions, kPcoTypeU8,
        requested_threads, codes);
      for (std::size_t row = 0; row < codes.size(); ++row) {
        const int code = static_cast<int>(codes[row]);
        if (code > 15 || (code >> 2) == (code & 3)) {
          throw std::runtime_error("native Pcodec substitution code is invalid");
        }
        INTEGER(substitution)[row] = static_cast<int>(codes[row]);
      }
    }

    if (need_z) {
      std::vector<std::uint16_t> codes(static_cast<std::size_t>(n));
      native_decompress_blocks_parallel<std::uint16_t>(
        native_path(files, 2, "z"), z_values, kPcoTypeU16,
        requested_threads, codes);
      for (std::size_t row = 0; row < codes.size(); ++row) {
        // Widest supported semantic Z profile (Z12: 4094 central codes plus
        // missing and exception). The store's own domain is checked against
        // its manifest in pcodec_native_validate_code_domains().
        if (codes[row] > 4095) {
          throw std::runtime_error("native Pcodec Z code is outside its domain");
        }
        INTEGER(z)[row] = static_cast<int>(codes[row]);
      }
    }
    if (need_se) {
      std::vector<std::uint8_t> codes(static_cast<std::size_t>(n));
      native_decompress_blocks_parallel<std::uint8_t>(
        native_path(files, 4, "SE"), se_values, kPcoTypeU8,
        requested_threads, codes);
      for (std::size_t row = 0; row < codes.size(); ++row) {
        if (codes[row] > 255) {
          throw std::runtime_error("native Pcodec SE code is outside its domain");
        }
        INTEGER(se)[row] = static_cast<int>(codes[row]);
      }
    }
    if (need_eaf) {
      std::vector<std::uint8_t> codes(static_cast<std::size_t>(n));
      native_decompress_blocks_parallel<std::uint8_t>(
        native_path(files, 3, "EAF"), eaf_values, kPcoTypeU8,
        requested_threads, codes);
      for (std::size_t row = 0; row < codes.size(); ++row) {
        if (codes[row] > 255) {
          throw std::runtime_error("native Pcodec EAF code is outside its domain");
        }
        INTEGER(eaf)[row] = static_cast<int>(codes[row]);
      }
    }

    std::vector<int> exception_rows;
    std::vector<double> exception_z;
    std::vector<double> exception_log2se;
    std::vector<double> exception_eaf;
    std::vector<int> exception_flags;
    if (need_exceptions) {
      const std::string codec = CHAR(STRING_ELT(exception_codec, 0));
      std::int64_t previous_row = -1;
      for (std::size_t block_id = 0; block_id < exceptions.size(); ++block_id) {
        native_append_exception_block(*exception_file, exceptions[block_id],
          z_values[block_id], static_cast<std::uint64_t>(n), &previous_row, codec,
          exception_rows, exception_z, exception_log2se, exception_eaf,
          exception_flags);
      }
    }

    SEXP exception_output = PROTECT(Rf_allocVector(VECSXP, 5));
    ++protect_count;
    SEXP exception_row = PROTECT(Rf_allocVector(INTSXP, exception_rows.size()));
    ++protect_count;
    SEXP exception_z_output = PROTECT(Rf_allocVector(REALSXP, exception_z.size()));
    ++protect_count;
    SEXP exception_log2se_output = PROTECT(Rf_allocVector(REALSXP, exception_log2se.size()));
    ++protect_count;
    SEXP exception_eaf_output = PROTECT(Rf_allocVector(REALSXP, exception_eaf.size()));
    ++protect_count;
    SEXP exception_flags_output = PROTECT(Rf_allocVector(INTSXP, exception_flags.size()));
    ++protect_count;
    for (std::size_t i = 0; i < exception_rows.size(); ++i) {
      INTEGER(exception_row)[i] = exception_rows[i];
      REAL(exception_z_output)[i] = exception_z[i];
      REAL(exception_log2se_output)[i] = exception_log2se[i];
      REAL(exception_eaf_output)[i] = exception_eaf[i];
      INTEGER(exception_flags_output)[i] = exception_flags[i];
    }
    SET_VECTOR_ELT(exception_output, 0, exception_row);
    SET_VECTOR_ELT(exception_output, 1, exception_z_output);
    SET_VECTOR_ELT(exception_output, 2, exception_log2se_output);
    SET_VECTOR_ELT(exception_output, 3, exception_eaf_output);
    SET_VECTOR_ELT(exception_output, 4, exception_flags_output);
    SEXP exception_names = PROTECT(Rf_allocVector(STRSXP, 5));
    ++protect_count;
    SET_STRING_ELT(exception_names, 0, scalar_name("row"));
    SET_STRING_ELT(exception_names, 1, scalar_name("z"));
    SET_STRING_ELT(exception_names, 2, scalar_name("log2se"));
    SET_STRING_ELT(exception_names, 3, scalar_name("eaf"));
    SET_STRING_ELT(exception_names, 4, scalar_name("flags"));
    Rf_setAttrib(exception_output, R_NamesSymbol, exception_names);

    SEXP output = PROTECT(Rf_allocVector(VECSXP, 6));
    ++protect_count;
    SET_VECTOR_ELT(output, 0, position);
    SET_VECTOR_ELT(output, 1, substitution);
    SET_VECTOR_ELT(output, 2, z);
    SET_VECTOR_ELT(output, 3, se);
    SET_VECTOR_ELT(output, 4, eaf);
    SET_VECTOR_ELT(output, 5, exception_output);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
    ++protect_count;
    SET_STRING_ELT(names, 0, scalar_name("position"));
    SET_STRING_ELT(names, 1, scalar_name("substitution"));
    SET_STRING_ELT(names, 2, scalar_name("z"));
    SET_STRING_ELT(names, 3, scalar_name("se"));
    SET_STRING_ELT(names, 4, scalar_name("eaf"));
    SET_STRING_ELT(names, 5, scalar_name("exceptions"));
    Rf_setAttrib(output, R_NamesSymbol, names);
    UNPROTECT(protect_count);
    return output;
  }, "unknown native Pcodec stream reader error");
#else
  Rf_error("native Pcodec is not available in this build");
#endif
  return R_NilValue;
}

#ifdef COMPRESSOR_NATIVE_PCODEC
namespace {

// Per-worker lazily opened stream files: each worker opens every stream file
// at most once and reuses the handle for all of its blocks.
class SelectWorkerFiles {
 public:
  explicit SelectWorkerFiles(const std::vector<std::string>& paths)
    : paths_(paths), files_(paths.size()) {}
  NativePcodecFile& get(std::size_t i) {
    if (!files_[i]) files_[i].reset(new NativePcodecFile(paths_[i]));
    return *files_[i];
  }
 private:
  const std::vector<std::string>& paths_;
  std::vector<std::unique_ptr<NativePcodecFile>> files_;
};

template <typename Fn>
void select_run_tasks(std::size_t task_count, std::size_t requested_threads,
                      const std::vector<std::string>& paths, Fn fn) {
  if (task_count == 0) return;
  const std::size_t workers_n = std::max<std::size_t>(
    1, std::min<std::size_t>(requested_threads, task_count));
  if (workers_n == 1) {
    SelectWorkerFiles files(paths);
    for (std::size_t t = 0; t < task_count; ++t) fn(t, files);
    return;
  }
  std::atomic<std::size_t> next(0);
  std::mutex error_mutex;
  std::exception_ptr first_error;
  std::vector<std::thread> workers;
  workers.reserve(workers_n);
  for (std::size_t w = 0; w < workers_n; ++w) {
    workers.emplace_back([&]() {
      try {
        SelectWorkerFiles files(paths);
        while (true) {
          const std::size_t t = next.fetch_add(1, std::memory_order_relaxed);
          if (t >= task_count) break;
          fn(t, files);
        }
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!first_error) first_error = std::current_exception();
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  if (first_error) std::rethrow_exception(first_error);
}

struct SelectKeyOut {
  std::vector<int> rows;
  std::vector<double> position;
  std::vector<int> substitution;
};

struct SelectExcOut {
  std::vector<int> index;  // zero-based index into the selection
  std::vector<double> z, log2se, eaf;
  std::vector<int> flags;
};

template <typename T>
std::vector<T> select_decode(NativePcodecFile& file, const NativePcodecBlock& block,
                             unsigned char dtype) {
  return native_decompress_block<T>(file, block, dtype);
}

// First index in [lo, n) with values[index] >= target, searching outward
// from lo in doubling steps (galloping) and finishing with a binary search.
template <typename Get>
std::size_t select_gallop(std::size_t lo, std::size_t n, double target, Get get) {
  if (lo >= n || get(lo) >= target) return lo;
  std::size_t step = 1;
  std::size_t prev = lo;  // get(prev) < target
  std::size_t probe = lo + 1;
  while (probe < n && get(probe) < target) {
    prev = probe;
    step <<= 1;
    probe = lo + step;
  }
  std::size_t first = prev + 1;
  std::size_t last = std::min(probe, n);
  while (first < last) {
    const std::size_t mid = first + (last - first) / 2;
    if (get(mid) < target) first = mid + 1; else last = mid;
  }
  return first;
}

// Rows of one decoded key block whose code (position * 16 + substitution) is
// in the sorted unique selection `sel`, in row order. Key blocks are sorted by
// code, so a galloping merge of the two sorted sequences replaces one binary
// search of the whole selection per row; a block that is not sorted by code
// keeps the per-row search. Duplicate codes in a block are all taken, as the
// per-row search does.
template <typename Take>
void select_key_join(const std::vector<double>& pos,
                     const std::vector<std::uint8_t>& subs,
                     const double* sel, std::size_t sel_len, Take take) {
  const std::size_t nb = pos.size();
  if (!nb || !sel_len) return;
  std::vector<double> code(nb);
  bool sorted = true;
  for (std::size_t i = 0; i < nb; ++i) {
    code[i] = pos[i] * 16 + static_cast<double>(subs[i]);
    if (i && code[i] < code[i - 1]) sorted = false;
  }
  if (!sorted) {
    for (std::size_t i = 0; i < nb; ++i) {
      if (std::binary_search(sel, sel + sel_len, code[i])) take(i);
    }
    return;
  }
  auto block_at = [&](std::size_t k) { return code[k]; };
  auto sel_at = [&](std::size_t k) { return sel[k]; };
  std::size_t j = static_cast<std::size_t>(
    std::lower_bound(sel, sel + sel_len, code[0]) - sel);
  std::size_t i = 0;
  while (i < nb && j < sel_len) {
    if (code[i] == sel[j]) {
      take(i);
      ++i;  // keep j: the next row may repeat this code
    } else if (code[i] < sel[j]) {
      i = select_gallop(i + 1, nb, sel[j], block_at);
    } else {
      j = select_gallop(j + 1, sel_len, code[i], sel_at);
    }
  }
}

SEXP select_int_vector(const std::vector<int>& v) {
  SEXP out = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(v.size())));
  if (!v.empty()) std::memcpy(INTEGER(out), v.data(), v.size() * sizeof(int));
  UNPROTECT(1);
  return out;
}
SEXP select_real_vector(const std::vector<double>& v) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(v.size())));
  if (!v.empty()) std::memcpy(REAL(out), v.data(), v.size() * sizeof(double));
  UNPROTECT(1);
  return out;
}

}  // namespace
#endif

// Selective native reader. One selection (mode 0: sorted unique zero-based
// row ids; mode 1: inclusive global-position window c(lower, upper); mode 2:
// sorted unique numeric keys position * 16 + substitution). Decodes only the
// blocks the selection touches, opening each stream file once per worker, and
// returns raw codes (and exception records) for the selected rows; the exact
// R-side arithmetic then turns codes into values.
extern "C" SEXP compressor_read_pcodec_native_select(
    SEXP files,
    SEXP position_blocks,
    SEXP substitution_blocks,
    SEXP z_blocks,
    SEXP eaf_blocks,
    SEXP se_blocks,
    SEXP exception_blocks,
    SEXP row_count,
    SEXP streams,
    SEXP exception_codec,
    SEXP threads,
    SEXP mode_sexp,
    SEXP selection,
    SEXP need_identity_sexp,
    SEXP delta_sexp) {
#ifdef COMPRESSOR_NATIVE_PCODEC
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(files) != STRSXP || XLENGTH(files) != 6 ||
        TYPEOF(streams) != STRSXP || TYPEOF(exception_codec) != STRSXP ||
        XLENGTH(exception_codec) != 1 || STRING_ELT(exception_codec, 0) == NA_STRING ||
        TYPEOF(selection) != REALSXP ||
        TYPEOF(need_identity_sexp) != LGLSXP || TYPEOF(delta_sexp) != LGLSXP) {
      throw std::runtime_error("malformed arguments to native Pcodec selective reader");
    }
    const double row_count_value = native_scalar_integer(
      row_count, true, "native Pcodec row count");
    const R_xlen_t n = static_cast<R_xlen_t>(row_count_value);
    const double thread_number = native_scalar_integer(
      threads, false, "native Pcodec thread count");
    const std::size_t requested_threads = static_cast<std::size_t>(thread_number);
    const int mode = static_cast<int>(native_scalar_integer(mode_sexp, true, "native select mode"));
    const bool need_identity = LOGICAL(need_identity_sexp)[0] == 1;
    const bool delta = LOGICAL(delta_sexp)[0] == 1;
    const bool need_z = requested_has(streams, "z");
    const bool need_se = requested_has(streams, "se");
    const bool need_eaf = need_se || requested_has(streams, "eaf");
    const bool need_numeric = need_z || need_se || need_eaf;
    if (mode < 0 || mode > 2 || (mode > 0 && !need_identity)) {
      throw std::runtime_error("invalid native select mode");
    }

    std::vector<NativePcodecBlock> positions = read_native_blocks(position_blocks, n, true);
    std::vector<NativePcodecBlock> substitutions = read_native_blocks(substitution_blocks, n, true);
    std::vector<NativePcodecBlock> z_values = read_native_blocks(z_blocks, n, false);
    std::vector<NativePcodecBlock> eaf_values = read_native_blocks(eaf_blocks, n, false);
    std::vector<NativePcodecBlock> se_values = read_native_blocks(se_blocks, n, false);
    std::vector<NativePcodecExceptionBlock> exceptions = read_native_exception_blocks(exception_blocks);
    require_block_row_parity(positions, substitutions, "native key stream");
    require_block_row_parity(z_values, eaf_values, "native value stream");
    require_block_row_parity(z_values, se_values, "native value stream");
    if (need_numeric && exceptions.size() != z_values.size()) {
      throw std::runtime_error("native Pcodec exception index does not match value blocks");
    }
    std::vector<std::string> paths;
    const char* labels[6] = {"position", "substitution", "z", "EAF", "SE", "exception"};
    for (int i = 0; i < 6; ++i) paths.push_back(native_path(files, i, labels[i]));
    const std::string codec = CHAR(STRING_ELT(exception_codec, 0));

    // last_position is carried as column 7 of the position block matrix.
    std::vector<double> last_position(positions.size(), 0.0);
    if (mode == 1 || mode == 2) {
      const SEXP dims = Rf_getAttrib(position_blocks, R_DimSymbol);
      if (INTEGER(dims)[1] < 7) throw std::runtime_error("native select needs last positions");
      for (std::size_t b = 0; b < positions.size(); ++b) {
        last_position[b] = static_cast<double>(matrix_uint64(
          position_blocks, static_cast<R_xlen_t>(b), 6, "native Pcodec block last position"));
      }
    }

    const double* sel = REAL(selection);
    const std::size_t sel_len = static_cast<std::size_t>(XLENGTH(selection));
    std::vector<int> sel_rows;  // mode 0 only
    std::vector<double> target_positions;  // mode 2
    if (mode == 0) {
      sel_rows.reserve(sel_len);
      for (std::size_t i = 0; i < sel_len; ++i) {
        if (!R_FINITE(sel[i]) || sel[i] < 0 || sel[i] >= static_cast<double>(n) ||
            (i > 0 && sel[i] <= sel[i - 1])) {
          throw std::runtime_error("native select rows must be sorted unique valid row ids");
        }
        sel_rows.push_back(static_cast<int>(sel[i]));
      }
    } else if (mode == 1) {
      if (sel_len != 2) throw std::runtime_error("native select window needs two bounds");
    } else {
      for (std::size_t i = 0; i < sel_len; ++i) {
        if (i > 0 && sel[i] <= sel[i - 1]) {
          throw std::runtime_error("native select keys must be sorted unique");
        }
        const double p = std::floor(sel[i] / 16);
        if (target_positions.empty() || target_positions.back() != p) target_positions.push_back(p);
      }
    }

    double source_bytes = 0;
    std::vector<int> selected;           // selected rows (sorted)
    std::vector<double> sel_position;
    std::vector<int> sel_substitution;

    auto block_hits = [&](const std::vector<NativePcodecBlock>& blocks,
                          const std::vector<int>& rows,
                          std::vector<std::size_t>& out) {
      std::size_t cursor = 0;
      for (std::size_t b = 0; b < blocks.size() && cursor < rows.size(); ++b) {
        const std::size_t lo = static_cast<std::size_t>(
          std::lower_bound(rows.begin() + cursor, rows.end(),
                           static_cast<int>(blocks[b].row_start)) - rows.begin());
        const std::size_t hi = static_cast<std::size_t>(
          std::lower_bound(rows.begin() + lo, rows.end(),
                           static_cast<int>(blocks[b].row_stop)) - rows.begin());
        if (hi > lo) out.push_back(b);
        cursor = hi;
      }
    };

    if (need_identity) {
      std::vector<std::size_t> key_candidates;
      if (mode == 0) {
        block_hits(positions, sel_rows, key_candidates);
      } else {
        for (std::size_t b = 0; b < positions.size(); ++b) {
          const double first = static_cast<double>(positions[b].first_position);
          if (mode == 1) {
            if (last_position[b] >= sel[0] && first <= sel[1]) key_candidates.push_back(b);
          } else {
            const auto lo = std::lower_bound(target_positions.begin(), target_positions.end(), first);
            const auto hi = std::upper_bound(target_positions.begin(), target_positions.end(), last_position[b]);
            if (lo < hi) key_candidates.push_back(b);
          }
        }
      }
      std::vector<SelectKeyOut> key_out(key_candidates.size());
      select_run_tasks(key_candidates.size(), requested_threads, paths,
        [&](std::size_t t, SelectWorkerFiles& wf) {
          const std::size_t b = key_candidates[t];
          const NativePcodecBlock& pb = positions[b];
          const NativePcodecBlock& sb = substitutions[b];
          std::vector<std::uint32_t> gaps = select_decode<std::uint32_t>(wf.get(0), pb, kPcoTypeU32);
          std::vector<std::uint8_t> subs = select_decode<std::uint8_t>(wf.get(1), sb, kPcoTypeU8);
          std::vector<double> pos(gaps.size());
          if (delta) {
            double current = static_cast<double>(pb.first_position);
            for (std::size_t i = 0; i < gaps.size(); ++i) {
              current += static_cast<double>(gaps[i]);
              pos[i] = current;
            }
          } else {
            for (std::size_t i = 0; i < gaps.size(); ++i) pos[i] = static_cast<double>(gaps[i]);
          }
          SelectKeyOut& out = key_out[t];
          auto take = [&](std::size_t i) {
            out.rows.push_back(static_cast<int>(pb.row_start + i));
            out.position.push_back(pos[i]);
            out.substitution.push_back(static_cast<int>(subs[i]));
          };
          if (mode == 0) {
            auto it = std::lower_bound(sel_rows.begin(), sel_rows.end(), static_cast<int>(pb.row_start));
            for (; it != sel_rows.end() && static_cast<std::uint64_t>(*it) < pb.row_stop; ++it) {
              take(static_cast<std::size_t>(*it) - static_cast<std::size_t>(pb.row_start));
            }
          } else if (mode == 1) {
            for (std::size_t i = 0; i < pos.size(); ++i) {
              if (pos[i] >= sel[0] && pos[i] <= sel[1]) take(i);
            }
          } else {
            select_key_join(pos, subs, sel, sel_len, take);
          }
        });
      for (std::size_t t = 0; t < key_candidates.size(); ++t) {
        const std::size_t b = key_candidates[t];
        source_bytes += static_cast<double>(positions[b].length) +
          static_cast<double>(substitutions[b].length);
        const SelectKeyOut& out = key_out[t];
        selected.insert(selected.end(), out.rows.begin(), out.rows.end());
        sel_position.insert(sel_position.end(), out.position.begin(), out.position.end());
        sel_substitution.insert(sel_substitution.end(), out.substitution.begin(), out.substitution.end());
      }
    } else {
      selected = sel_rows;
    }

    const std::size_t m = selected.size();
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 12));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 12));
    const char* result_names[12] = {"rows", "position", "substitution", "z", "eaf", "se",
                                    "exc_index", "exc_z", "exc_log2se", "exc_eaf",
                                    "exc_flags", "source_bytes"};
    for (int i = 0; i < 12; ++i) SET_STRING_ELT(names, i, scalar_name(result_names[i]));
    Rf_setAttrib(result, R_NamesSymbol, names);
    SET_VECTOR_ELT(result, 0, select_int_vector(selected));
    if (need_identity) {
      SET_VECTOR_ELT(result, 1, select_real_vector(sel_position));
      SET_VECTOR_ELT(result, 2, select_int_vector(sel_substitution));
    }

    if (m > 0 && need_numeric) {
      SEXP z_out = need_z ? Rf_allocVector(INTSXP, static_cast<R_xlen_t>(m)) : R_NilValue;
      if (need_z) SET_VECTOR_ELT(result, 3, z_out);
      SEXP eaf_out = need_eaf ? Rf_allocVector(INTSXP, static_cast<R_xlen_t>(m)) : R_NilValue;
      if (need_eaf) SET_VECTOR_ELT(result, 4, eaf_out);
      SEXP se_out = need_se ? Rf_allocVector(INTSXP, static_cast<R_xlen_t>(m)) : R_NilValue;
      if (need_se) SET_VECTOR_ELT(result, 5, se_out);
      int* zp = need_z ? INTEGER(z_out) : nullptr;
      int* ep = need_eaf ? INTEGER(eaf_out) : nullptr;
      int* sp = need_se ? INTEGER(se_out) : nullptr;

      std::vector<std::size_t> value_candidates;
      block_hits(z_values, selected, value_candidates);
      std::vector<SelectExcOut> exc_out(value_candidates.size());
      select_run_tasks(value_candidates.size(), requested_threads, paths,
        [&](std::size_t t, SelectWorkerFiles& wf) {
          const std::size_t b = value_candidates[t];
          const NativePcodecBlock& vb = z_values[b];
          const std::size_t lo = static_cast<std::size_t>(std::lower_bound(
            selected.begin(), selected.end(), static_cast<int>(vb.row_start)) - selected.begin());
          const std::size_t hi = static_cast<std::size_t>(std::lower_bound(
            selected.begin() + lo, selected.end(), static_cast<int>(vb.row_stop)) - selected.begin());
          if (need_z) {
            std::vector<std::uint16_t> v = select_decode<std::uint16_t>(wf.get(2), z_values[b], kPcoTypeU16);
            for (std::size_t i = lo; i < hi; ++i) zp[i] = static_cast<int>(v[selected[i] - vb.row_start]);
          }
          if (need_eaf) {
            std::vector<std::uint8_t> v = select_decode<std::uint8_t>(wf.get(3), eaf_values[b], kPcoTypeU8);
            for (std::size_t i = lo; i < hi; ++i) ep[i] = static_cast<int>(v[selected[i] - vb.row_start]);
          }
          if (need_se) {
            std::vector<std::uint8_t> v = select_decode<std::uint8_t>(wf.get(4), se_values[b], kPcoTypeU8);
            for (std::size_t i = lo; i < hi; ++i) sp[i] = static_cast<int>(v[selected[i] - vb.row_start]);
          }
          if (exceptions[b].count > 0) {
            std::vector<int> rows, flags;
            std::vector<double> ez, el, ee;
            native_append_exception_block(wf.get(5), exceptions[b], vb,
              static_cast<std::uint64_t>(n), nullptr, codec, rows, ez, el, ee, flags);
            SelectExcOut& out = exc_out[t];
            for (std::size_t k = 0; k < rows.size(); ++k) {
              auto it = std::lower_bound(selected.begin() + lo, selected.begin() + hi, rows[k]);
              if (it != selected.begin() + hi && *it == rows[k]) {
                out.index.push_back(static_cast<int>(it - selected.begin()));
                out.z.push_back(ez[k]);
                out.log2se.push_back(el[k]);
                out.eaf.push_back(ee[k]);
                out.flags.push_back(flags[k]);
              }
            }
          }
        });
      SelectExcOut all;
      for (std::size_t t = 0; t < value_candidates.size(); ++t) {
        const std::size_t b = value_candidates[t];
        if (need_z) source_bytes += static_cast<double>(z_values[b].length);
        if (need_eaf) source_bytes += static_cast<double>(eaf_values[b].length);
        if (need_se) source_bytes += static_cast<double>(se_values[b].length);
        if (exceptions[b].count > 0) source_bytes += static_cast<double>(exceptions[b].length);
        const SelectExcOut& o = exc_out[t];
        all.index.insert(all.index.end(), o.index.begin(), o.index.end());
        all.z.insert(all.z.end(), o.z.begin(), o.z.end());
        all.log2se.insert(all.log2se.end(), o.log2se.begin(), o.log2se.end());
        all.eaf.insert(all.eaf.end(), o.eaf.begin(), o.eaf.end());
        all.flags.insert(all.flags.end(), o.flags.begin(), o.flags.end());
      }
      SET_VECTOR_ELT(result, 6, select_int_vector(all.index));
      SET_VECTOR_ELT(result, 7, select_real_vector(all.z));
      SET_VECTOR_ELT(result, 8, select_real_vector(all.log2se));
      SET_VECTOR_ELT(result, 9, select_real_vector(all.eaf));
      SET_VECTOR_ELT(result, 10, select_int_vector(all.flags));
    }
    SET_VECTOR_ELT(result, 11, Rf_ScalarReal(source_bytes));
    UNPROTECT(2);
    return result;
  }, "unknown native Pcodec selective reader error");
#else
  Rf_error("native Pcodec is not available in this build");
#endif
  return R_NilValue;
}



// Domain checks for decoded code vectors, in one allocation-free pass per
// stream. Returns 0 when every code is valid, otherwise the number of the
// first failing check (1 Z, 2 SE, 3 EAF, 4 position order/uint32 domain,
// 5 position beyond the genome table, 6 substitution); R maps the number to
// the error message. Any stream may be NULL to skip it.
extern "C" SEXP compressor_validate_native_codes(
    SEXP z, SEXP se, SEXP eaf, SEXP position, SEXP substitution,
    SEXP limits) {
  return compressor_guard([&]() -> SEXP {
    if (TYPEOF(limits) != REALSXP || XLENGTH(limits) != 4) {
      compressor_fail("malformed limits for native Pcodec code validation");
    }
    const double* limit = REAL(limits);
    auto check_int = [](SEXP codes, double max_code) -> bool {
      if (codes == R_NilValue) return true;
      if (TYPEOF(codes) != INTSXP) return false;
      const int* value = INTEGER(codes);
      const R_xlen_t n = XLENGTH(codes);
      bool ok = true;
      for (R_xlen_t i = 0; i < n; ++i) {
        ok &= value[i] >= 0 && static_cast<double>(value[i]) <= max_code;
      }
      return ok;
    };
    int status = 0;
    if (!check_int(z, limit[0])) status = 1;
    else if (!check_int(se, limit[1])) status = 2;
    else if (!check_int(eaf, limit[2])) status = 3;
    if (status == 0 && position != R_NilValue) {
      if (TYPEOF(position) != REALSXP) {
        status = 4;
      } else {
        const double* value = REAL(position);
        const R_xlen_t n = XLENGTH(position);
        double previous = 0.0;
        bool ordered = true;
        bool inside = true;
        for (R_xlen_t i = 0; i < n; ++i) {
          const double current = value[i];
          // NaN fails every comparison, so it is rejected here as well.
          ordered &= current >= 0.0 && current <= 4294967295.0 &&
                     (i == 0 || current >= previous);
          inside &= !(current >= limit[3]);
          previous = current;
        }
        if (!ordered) status = 4;
        else if (!inside) status = 5;
      }
    }
    if (status == 0 && substitution != R_NilValue) {
      if (TYPEOF(substitution) != INTSXP) {
        status = 6;
      } else {
        const int* value = INTEGER(substitution);
        const R_xlen_t n = XLENGTH(substitution);
        bool ok = true;
        for (R_xlen_t i = 0; i < n; ++i) {
          const int code = value[i];
          ok &= code >= 0 && code <= 15 && (code >> 2) != (code & 3);
        }
        if (!ok) status = 6;
      }
    }
    return Rf_ScalarInteger(status);
  }, "native Pcodec error");
}

// Build the identity columns of a full read from global positions and 4-bit
// substitution codes in one native pass. Strings are taken from small level
// tables (chromosome labels and the four bases), so every row stores a shared
// CHARSXP pointer; nothing is formatted or allocated per row. The semantics
// match the R definition in pcodec_native_key_columns(): the chromosome is
// findInterval(position, offsets) clamped to the table, the local position is
// position - offset + 1, REF is code >> 2 and ALT is code & 3.
extern "C" SEXP compressor_pcodec_key_columns(
    SEXP position, SEXP substitution, SEXP offsets, SEXP labels, SEXP wanted) {
  return compressor_guard([&]() -> SEXP {
    if ((TYPEOF(position) != REALSXP && TYPEOF(position) != INTSXP) ||
        TYPEOF(substitution) != INTSXP || TYPEOF(offsets) != REALSXP ||
        TYPEOF(labels) != STRSXP || TYPEOF(wanted) != LGLSXP ||
        XLENGTH(wanted) != 4 || XLENGTH(offsets) < 1 ||
        XLENGTH(offsets) != XLENGTH(labels)) {
      compressor_fail("malformed arguments to native Pcodec key column builder");
    }
    const R_xlen_t n = XLENGTH(position);
    const int* want = LOGICAL(wanted);
    const bool want_chromosome = want[0] == TRUE;
    const bool want_local = want[1] == TRUE;
    const bool want_reference = want[2] == TRUE;
    const bool want_alternate = want[3] == TRUE;
    if ((want_reference || want_alternate) && XLENGTH(substitution) != n) {
      compressor_fail("native Pcodec key columns need one substitution per position");
    }
    const R_xlen_t k = XLENGTH(offsets);
    const double* offset = REAL(offsets);
    for (R_xlen_t i = 0; i < k; ++i) {
      if (!R_FINITE(offset[i]) || (i > 0 && offset[i] < offset[i - 1])) {
        compressor_fail("native Pcodec chromosome offsets must be finite and sorted");
      }
    }
  
    int protect_count = 0;
    SEXP output = PROTECT(Rf_allocVector(VECSXP, 4));
    ++protect_count;
    SEXP chromosome = R_NilValue, local = R_NilValue;
    SEXP reference = R_NilValue, alternate = R_NilValue;
    if (want_chromosome) {
      chromosome = Rf_allocVector(STRSXP, n);
      SET_VECTOR_ELT(output, 0, chromosome);
    }
    if (want_local) {
      local = Rf_allocVector(INTSXP, n);
      SET_VECTOR_ELT(output, 1, local);
    }
    if (want_reference) {
      reference = Rf_allocVector(STRSXP, n);
      SET_VECTOR_ELT(output, 2, reference);
    }
    if (want_alternate) {
      alternate = Rf_allocVector(STRSXP, n);
      SET_VECTOR_ELT(output, 3, alternate);
    }
  
    if (want_chromosome || want_local) {
      const bool is_real = TYPEOF(position) == REALSXP;
      const double* real_position = is_real ? REAL(position) : nullptr;
      const int* int_position = is_real ? nullptr : INTEGER(position);
      int* local_out = want_local ? INTEGER(local) : nullptr;
      R_xlen_t current = 0;
      for (R_xlen_t row = 0; row < n; ++row) {
        double value;
        if (is_real) {
          value = real_position[row];
        } else {
          value = int_position[row] == NA_INTEGER ? NA_REAL
                                                  : static_cast<double>(int_position[row]);
        }
        if (ISNAN(value)) {
          if (want_chromosome) SET_STRING_ELT(chromosome, row, NA_STRING);
          if (want_local) local_out[row] = NA_INTEGER;
          continue;
        }
        // Sorted stores stay inside the current chromosome or move forward by
        // one; anything else falls back to a binary search.
        if (!(offset[current] <= value && (current + 1 == k || value < offset[current + 1]))) {
          if (current + 1 < k && value >= offset[current + 1] &&
              (current + 2 == k || value < offset[current + 2])) {
            ++current;
          } else {
            const double* hit = std::upper_bound(offset, offset + k, value);
            current = hit == offset ? 0 : static_cast<R_xlen_t>(hit - offset) - 1;
          }
        }
        if (want_chromosome) {
          SET_STRING_ELT(chromosome, row, STRING_ELT(labels, current));
        }
        if (want_local) {
          const double bp = value - offset[current] + 1.0;
          local_out[row] = (bp >= static_cast<double>(INT_MIN + 1) &&
                            bp <= static_cast<double>(INT_MAX))
            ? static_cast<int>(bp) : NA_INTEGER;
        }
      }
    }
  
    if (want_reference || want_alternate) {
      SEXP bases[4];
      const char* base_names[4] = {"A", "C", "G", "T"};
      for (int i = 0; i < 4; ++i) {
        bases[i] = PROTECT(Rf_mkChar(base_names[i]));
        ++protect_count;
      }
      const int* code = INTEGER(substitution);
      for (R_xlen_t row = 0; row < n; ++row) {
        // Valid codes are 0..15. Other values mirror the R fallback, which
        // indexes the base table with bitwShiftR()/bitwAnd() on the code.
        const int value = code[row];
        if (want_reference) {
          SET_STRING_ELT(reference, row, value >= 0 && value <= 15
            ? bases[value >> 2] : NA_STRING);
        }
        if (want_alternate) {
          SET_STRING_ELT(alternate, row, value != NA_INTEGER
            ? bases[static_cast<unsigned>(value) & 3u] : NA_STRING);
        }
      }
    }
    UNPROTECT(protect_count);
    return output;
  }, "native Pcodec error");
}

