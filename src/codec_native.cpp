#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <R.h>
#include <Rinternals.h>

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

std::string bridge_file(const std::string& directory, const char* name) {
  if (!directory.empty() && directory.back() == '/') return directory + name;
  return directory + "/" + name;
}

template <typename T>
std::vector<T> read_bridge_vector(const std::string& path,
                                  std::size_t expected = std::numeric_limits<std::size_t>::max()) {
  std::ifstream stream(path, std::ios::binary | std::ios::ate);
  if (!stream) throw std::runtime_error("cannot open Pcodec bridge file: " + path);
  const std::streamoff byte_count = stream.tellg();
  if (byte_count < 0 || byte_count % static_cast<std::streamoff>(sizeof(T)) != 0) {
    throw std::runtime_error("Pcodec bridge file has an invalid byte length: " + path);
  }
  const std::size_t count = static_cast<std::size_t>(byte_count / sizeof(T));
  if (expected != std::numeric_limits<std::size_t>::max() && count != expected) {
    throw std::runtime_error(
      "Pcodec bridge column has " + std::to_string(count) +
      " values; expected " + std::to_string(expected) + ": " + path);
  }
  std::vector<T> values(count);
  stream.seekg(0, std::ios::beg);
  if (byte_count > 0 &&
      !stream.read(reinterpret_cast<char*>(values.data()), byte_count)) {
    throw std::runtime_error("cannot read complete Pcodec bridge file: " + path);
  }
  if constexpr (sizeof(T) > 1) {
    const std::uint16_t endian_probe = 1;
    const bool little_endian =
      *reinterpret_cast<const std::uint8_t*>(&endian_probe) == 1;
    if (!little_endian) {
      for (T& value : values) {
        auto* first = reinterpret_cast<unsigned char*>(&value);
        std::reverse(first, first + sizeof(T));
      }
    }
  }
  return values;
}

bool requested_has(SEXP requested, const char* target) {
  const R_xlen_t n = XLENGTH(requested);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP value = STRING_ELT(requested, i);
    if (value == NA_STRING) {
      throw std::runtime_error("native Pcodec bridge requested columns contain NA");
    }
    if (std::string(CHAR(value)) == target) return true;
  }
  return false;
}

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
    SEXP se_residual_max) {
  if (TYPEOF(z_code) != INTSXP || TYPEOF(se_code) != INTSXP ||
      TYPEOF(eaf_code) != INTSXP) {
    Rf_error("native decoder requires integer code vectors");
  }
  const R_xlen_t n = XLENGTH(z_code);
  if (XLENGTH(se_code) != n || XLENGTH(eaf_code) != n) {
    Rf_error("native decoder code vectors must have equal lengths");
  }

  const int z_bits_value = Rf_asInteger(z_bits);
  const int se_bits_value = Rf_asInteger(se_bits);
  const int eaf_bits_value = Rf_asInteger(eaf_bits);
  const int z_count_value = Rf_asInteger(z_count);
  const int se_count_value = Rf_asInteger(se_count);
  const int eaf_count_value = Rf_asInteger(eaf_count);
  const double se_residual_min_value = Rf_asReal(se_residual_min);
  const double se_residual_max_value = Rf_asReal(se_residual_max);
  const int block_rows_value = std::max(1, Rf_asInteger(block_rows));
  if (z_bits_value <= 0 || se_bits_value <= 0 || eaf_bits_value <= 0 ||
      z_count_value <= 0 || se_count_value <= 0 || eaf_count_value <= 0) {
    Rf_error("native decoder metadata contains invalid code domains");
  }
  if (!R_FINITE(se_residual_min_value) || !R_FINITE(se_residual_max_value) ||
      se_residual_max_value <= se_residual_min_value) {
    Rf_error("native decoder metadata contains an invalid SE residual range");
  }
  if (z_bits_value > 16 || se_bits_value > 16 || eaf_bits_value > 16) {
    Rf_error("native decoder code domains exceed the supported 16-bit limit");
  }

  const std::size_t z_width = static_cast<std::size_t>(1) << z_bits_value;
  const std::size_t se_width = static_cast<std::size_t>(1) << se_bits_value;
  const std::size_t eaf_width = static_cast<std::size_t>(1) << eaf_bits_value;
  if (z_width > 65536 || se_width > 65536 || eaf_width > 65536) {
    Rf_error("native decoder code domains are too large");
  }

  std::vector<double> z_table(z_width, kNaN);
  const double z_step = (Rf_asReal(z_max) - Rf_asReal(z_min)) /
                        static_cast<double>(z_count_value);
  for (int code = 0; code < z_count_value &&
                       static_cast<std::size_t>(code) < z_width; ++code) {
    z_table[static_cast<std::size_t>(code)] = Rf_asReal(z_min) +
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
  std::vector<double> se_table(se_width * eaf_width, kNaN);
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
    Rf_error("native decoder block centres must be numeric");
  }
  // The block centre is constant for every row in a block.  Materialise its
  // multiplicative factor once; calling exp2() in the row loop made the
  // R-facing decoder materially slower than the historical native reader.
  std::vector<double> centre_factors(static_cast<std::size_t>(centre_count));
  for (R_xlen_t i = 0; i < centre_count; ++i) {
    centre_factors[static_cast<std::size_t>(i)] = std::exp2(centre_values[i]);
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
      const R_xlen_t block = row / block_rows_value;
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
    Rf_error("native decoder exception vectors are malformed");
  }
  const int* exception_rows = INTEGER(exception_row);
  const double* exception_z_values = REAL(exception_z);
  const double* exception_se_values = REAL(exception_se);
  const double* exception_eaf_values = REAL(exception_eaf);
  const int* exception_flag_values = INTEGER(exception_flags);
  for (R_xlen_t i = 0; i < exception_count; ++i) {
    const int row = exception_rows[i];
    if (row < 0 || static_cast<R_xlen_t>(row) >= n) continue;
    if (exception_flag_values[i] & 1) z_out[row] = exception_z_values[i];
    if (exception_flag_values[i] & 2) se_out[row] = exception_se_values[i];
    if (exception_flag_values[i] & 4) eaf_out[row] = exception_eaf_values[i];
    if (want_beta) beta_out[row] = z_out[row] * se_out[row];
    if (want_p && (exception_flag_values[i] & 1)) {
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
}

// Load the compact temporary bridge directly into final R vectors. This avoids
// expanding uint16/uint8 codes to 32-bit R integer vectors, avoids readBin()
// copies, and constructs categorical identity columns in the same native pass.
SEXP compressor_read_pcodec_bridge_impl(
    SEXP directory,
    SEXP requested,
    SEXP row_count,
    SEXP centres,
    SEXP z_min,
    SEXP z_max,
    SEXP z_count,
    SEXP se_count,
    SEXP eaf_count,
    SEXP block_rows,
    SEXP chromosome_lengths) {
  if (TYPEOF(directory) != STRSXP || XLENGTH(directory) != 1 ||
      TYPEOF(requested) != STRSXP || TYPEOF(centres) != REALSXP ||
      TYPEOF(chromosome_lengths) != REALSXP) {
    throw std::runtime_error("malformed arguments to native Pcodec bridge reader");
  }
  if (STRING_ELT(directory, 0) == NA_STRING) {
    throw std::runtime_error("native Pcodec bridge directory is NA");
  }
  const std::string root = CHAR(STRING_ELT(directory, 0));
  const double row_count_value = Rf_asReal(row_count);
  if (!R_FINITE(row_count_value) || row_count_value < 0.0 ||
      row_count_value != std::floor(row_count_value) ||
      static_cast<long double>(row_count_value) >
        static_cast<long double>(std::numeric_limits<R_xlen_t>::max())) {
    throw std::runtime_error("Pcodec bridge row count is invalid");
  }
  const R_xlen_t n = static_cast<R_xlen_t>(row_count_value);
  const std::size_t rows = static_cast<std::size_t>(n);
  const int z_count_value = Rf_asInteger(z_count);
  const int se_count_value = Rf_asInteger(se_count);
  const int eaf_count_value = Rf_asInteger(eaf_count);
  const int block_rows_value = std::max(1, Rf_asInteger(block_rows));
  if (z_count_value != 510 || se_count_value != 62 || eaf_count_value != 255 ||
      block_rows_value != 65536 || Rf_asReal(z_min) != -3.5 || Rf_asReal(z_max) != 3.5) {
    throw std::runtime_error("native Pcodec bridge codec constants are invalid");
  }

  const bool want_chromosome = requested_has(requested, "chromosome");
  const bool want_position = requested_has(requested, "base_pair_location");
  const bool want_effect = requested_has(requested, "effect_allele");
  const bool want_other = requested_has(requested, "other_allele");
  const bool want_z = requested_has(requested, "z");
  const bool want_beta = requested_has(requested, "beta");
  const bool want_se = requested_has(requested, "standard_error");
  const bool want_eaf = requested_has(requested, "effect_allele_frequency");
  const bool want_p = requested_has(requested, "p_value");
  const bool need_z = want_z || want_beta || want_p;
  const bool need_se = want_se || want_beta;
  const bool need_eaf = want_eaf || need_se;
  const bool need_numeric = need_z || need_se || need_eaf;
  const bool compact_identity = XLENGTH(chromosome_lengths) > 0;
  if (compact_identity && XLENGTH(chromosome_lengths) != 24) {
    throw std::runtime_error(
      "compact Pcodec identity requires 24 chromosome lengths");
  }

  std::vector<std::uint8_t> chromosome_codes;
  std::vector<std::int32_t> positions;
  std::vector<std::uint8_t> effect_codes;
  std::vector<std::uint8_t> other_codes;
  std::vector<std::uint32_t> global_positions;
  std::vector<std::uint8_t> substitution_codes;
  if (compact_identity &&
      (want_chromosome || want_position || want_effect || want_other)) {
    global_positions = read_bridge_vector<std::uint32_t>(
      bridge_file(root, "global_position_code.bin"), rows);
    substitution_codes = read_bridge_vector<std::uint8_t>(
      bridge_file(root, "substitution_code.bin"), rows);
  } else {
    if (want_chromosome) {
      chromosome_codes = read_bridge_vector<std::uint8_t>(
        bridge_file(root, "chromosome.bin"), rows);
    }
    if (want_position) {
      positions = read_bridge_vector<std::int32_t>(
        bridge_file(root, "base_pair_location.bin"), rows);
    }
    if (want_effect) {
      effect_codes = read_bridge_vector<std::uint8_t>(
        bridge_file(root, "effect_allele.bin"), rows);
    }
    if (want_other) {
      other_codes = read_bridge_vector<std::uint8_t>(
        bridge_file(root, "other_allele.bin"), rows);
    }
  }

  std::vector<std::uint16_t> z_codes;
  std::vector<std::uint8_t> se_codes;
  std::vector<std::uint8_t> eaf_codes;
  std::vector<std::int32_t> exception_rows;
  std::vector<float> exception_z;
  std::vector<double> exception_se;
  std::vector<float> exception_eaf;
  std::vector<std::uint8_t> exception_flags;
  if (need_numeric) {
    if (need_z) {
      z_codes = read_bridge_vector<std::uint16_t>(bridge_file(root, "z_code.bin"), rows);
    }
    if (need_se) {
      se_codes = read_bridge_vector<std::uint8_t>(bridge_file(root, "se_code.bin"), rows);
    }
    if (need_eaf) {
      eaf_codes = read_bridge_vector<std::uint8_t>(bridge_file(root, "eaf_code.bin"), rows);
    }
    exception_rows = read_bridge_vector<std::int32_t>(bridge_file(root, "exception_row.bin"));
    const std::size_t exception_count = exception_rows.size();
    exception_z = read_bridge_vector<float>(bridge_file(root, "exception_z.bin"), exception_count);
    exception_se = read_bridge_vector<double>(bridge_file(root, "exception_se.bin"), exception_count);
    exception_eaf = read_bridge_vector<float>(bridge_file(root, "exception_eaf.bin"), exception_count);
    exception_flags = read_bridge_vector<std::uint8_t>(
      bridge_file(root, "exception_flags.bin"), exception_count);
  }

  int protect_count = 0;
  SEXP chromosome = R_NilValue;
  SEXP position = R_NilValue;
  SEXP effect = R_NilValue;
  SEXP other = R_NilValue;
  SEXP z = R_NilValue;
  SEXP se = R_NilValue;
  SEXP eaf = R_NilValue;
  SEXP beta = R_NilValue;
  SEXP p = R_NilValue;

  if (want_chromosome) { chromosome = PROTECT(Rf_allocVector(STRSXP, n)); ++protect_count; }
  if (want_position) { position = PROTECT(Rf_allocVector(INTSXP, n)); ++protect_count; }
  if (want_effect) { effect = PROTECT(Rf_allocVector(STRSXP, n)); ++protect_count; }
  if (want_other) { other = PROTECT(Rf_allocVector(STRSXP, n)); ++protect_count; }
  if (need_z) { z = PROTECT(Rf_allocVector(REALSXP, n)); ++protect_count; }
  if (need_se) { se = PROTECT(Rf_allocVector(REALSXP, n)); ++protect_count; }
  if (need_eaf) { eaf = PROTECT(Rf_allocVector(REALSXP, n)); ++protect_count; }
  if (want_beta) { beta = PROTECT(Rf_allocVector(REALSXP, n)); ++protect_count; }
  if (want_p) { p = PROTECT(Rf_allocVector(REALSXP, n)); ++protect_count; }

  const char* chromosome_labels[] = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"
  };
  SEXP chromosome_chars[24];
  for (int i = 0; i < 24; ++i) {
    chromosome_chars[i] = scalar_name(chromosome_labels[i]);
  }
  SEXP bases[] = {
    scalar_name("A"), scalar_name("C"), scalar_name("G"), scalar_name("T")
  };
  if (compact_identity &&
      (want_chromosome || want_position || want_effect || want_other)) {
    std::vector<std::uint64_t> offsets(24);
    std::vector<std::uint64_t> lengths(24);
    std::uint64_t total_length = 0;
    const double* supplied_lengths = REAL(chromosome_lengths);
    for (std::size_t chromosome_index = 0; chromosome_index < 24;
         ++chromosome_index) {
      const double supplied = supplied_lengths[chromosome_index];
      if (!R_FINITE(supplied) || supplied < 1.0 ||
          supplied != std::floor(supplied)) {
        throw std::runtime_error(
          "compact Pcodec identity has an invalid chromosome length");
      }
      offsets[chromosome_index] = total_length;
      lengths[chromosome_index] = static_cast<std::uint64_t>(supplied);
      total_length += lengths[chromosome_index];
    }
    std::size_t chromosome_index = 0;
    std::uint64_t previous_global = 0;
    for (R_xlen_t row = 0; row < n; ++row) {
      const std::uint64_t global =
        global_positions[static_cast<std::size_t>(row)];
      if (global >= total_length) {
        throw std::runtime_error(
          "global position is outside the compact Pcodec genome");
      }
      if (row > 0 && global < previous_global) {
        throw std::runtime_error(
          "compact Pcodec global positions are not sorted");
      }
      previous_global = global;
      while (chromosome_index + 1 < offsets.size() &&
             global >= offsets[chromosome_index + 1]) {
        ++chromosome_index;
      }
      const std::uint64_t local = global - offsets[chromosome_index] + 1;
      if (local > lengths[chromosome_index] ||
          local > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(
          "local position is outside the compact Pcodec chromosome");
      }
      const int substitution =
        substitution_codes[static_cast<std::size_t>(row)];
      const int reference_code = substitution >> 2;
      const int alternate_code = substitution & 3;
      if (substitution < 0 || substitution > 15 ||
          reference_code == alternate_code) {
        throw std::runtime_error(
          "invalid substitution code in compact Pcodec identity");
      }
      if (want_chromosome) {
        SET_STRING_ELT(chromosome, row, chromosome_chars[chromosome_index]);
      }
      if (want_position) INTEGER(position)[row] = static_cast<int>(local);
      if (want_effect) SET_STRING_ELT(effect, row, bases[alternate_code]);
      if (want_other) SET_STRING_ELT(other, row, bases[reference_code]);
    }
  } else {
    for (R_xlen_t row = 0; row < n; ++row) {
      const std::size_t index = static_cast<std::size_t>(row);
      if (want_chromosome) {
        const int code = chromosome_codes[index];
        if (code < 1 || code > 24) {
          throw std::runtime_error(
            "invalid chromosome code in Pcodec bridge");
        }
        SET_STRING_ELT(chromosome, row, chromosome_chars[code - 1]);
      }
      if (want_position) {
        const int value = positions[index];
        if (value < 1) {
          throw std::runtime_error("invalid position in Pcodec bridge");
        }
        INTEGER(position)[row] = value;
      }
      if (want_effect) {
        const int code = effect_codes[index];
        if (code < 0 || code > 3) {
          throw std::runtime_error(
            "invalid effect-allele code in Pcodec bridge");
        }
        SET_STRING_ELT(effect, row, bases[code]);
      }
      if (want_other) {
        const int code = other_codes[index];
        if (code < 0 || code > 3) {
          throw std::runtime_error(
            "invalid other-allele code in Pcodec bridge");
        }
        SET_STRING_ELT(other, row, bases[code]);
      }
    }
  }

  if (need_numeric) {
    const std::size_t z_width = static_cast<std::size_t>(z_count_value + 2);
    const std::size_t se_width = static_cast<std::size_t>(se_count_value + 2);
    const std::size_t eaf_width = static_cast<std::size_t>(eaf_count_value + 1);
    std::vector<double> z_table(z_width, kNaN);
    std::vector<double> p_table(z_width, kNaN);
    const double z_step = (Rf_asReal(z_max) - Rf_asReal(z_min)) /
                          static_cast<double>(z_count_value);
    for (int code = 0; code < z_count_value; ++code) {
      const double value = Rf_asReal(z_min) + (static_cast<double>(code) + 0.5) * z_step;
      z_table[static_cast<std::size_t>(code)] = value;
      p_table[static_cast<std::size_t>(code)] = std::erfc(std::abs(value) * kInvSqrtTwo);
    }
    std::vector<double> eaf_table(eaf_width, kNaN);
    for (int code = 0; code <= eaf_count_value; ++code) {
      eaf_table[static_cast<std::size_t>(code)] = std::pow(
        std::sin((static_cast<double>(code) / static_cast<double>(eaf_count_value)) *
                 kPi / 2.0), 2.0);
    }
    std::vector<double> se_table(se_width * eaf_width, kNaN);
    for (int se_value = 0; se_value < se_count_value; ++se_value) {
      const double residual = -1.0 +
        (static_cast<double>(se_value) + 0.5) * (2.0 / se_count_value);
      for (std::size_t eaf_value = 0; eaf_value < eaf_width; ++eaf_value) {
        const double safe = clamp_probability_eaf(eaf_table[eaf_value]);
        const double correction = -0.5 * std::log2(2.0 * safe * (1.0 - safe));
        se_table[static_cast<std::size_t>(se_value) * eaf_width + eaf_value] =
          std::exp2(residual + correction);
      }
    }
    const R_xlen_t centre_count = XLENGTH(centres);
    const double* centre_values = REAL(centres);
    const R_xlen_t expected_centres = (n + block_rows_value - 1) / block_rows_value;
    if (need_se && centre_count != expected_centres) {
      throw std::runtime_error("Pcodec bridge has the wrong number of SE block centres");
    }
    for (R_xlen_t i = 0; i < centre_count; ++i) {
      if (!R_FINITE(centre_values[i])) {
        throw std::runtime_error("Pcodec bridge has a non-finite SE block centre");
      }
    }
    std::vector<double> centre_factors(static_cast<std::size_t>(centre_count));
    for (R_xlen_t i = 0; i < centre_count; ++i) {
      centre_factors[static_cast<std::size_t>(i)] = std::exp2(centre_values[i]);
    }
    double* z_out = need_z ? REAL(z) : nullptr;
    double* se_out = need_se ? REAL(se) : nullptr;
    double* eaf_out = need_eaf ? REAL(eaf) : nullptr;
    double* beta_out = want_beta ? REAL(beta) : nullptr;
    double* p_out = want_p ? REAL(p) : nullptr;
    for (R_xlen_t row = 0; row < n; ++row) {
      const std::size_t index = static_cast<std::size_t>(row);
      const int z_value = need_z ? z_codes[index] : z_count_value;
      const int se_value = need_se ? se_codes[index] : se_count_value;
      const int eaf_value = need_eaf ? eaf_codes[index] : 0;
      if ((need_z && (z_value < 0 || z_value > z_count_value + 1)) ||
          (need_se && (se_value < 0 || se_value > se_count_value + 1)) ||
          (need_eaf && (eaf_value < 0 || eaf_value > eaf_count_value))) {
        throw std::runtime_error("Pcodec bridge contains a code outside its declared domain");
      }
      const double decoded_z = need_z && z_value >= 0 && z_value < z_count_value
        ? z_table[static_cast<std::size_t>(z_value)] : kNaN;
      const double decoded_eaf = need_eaf && eaf_value >= 0 && eaf_value <= eaf_count_value
        ? eaf_table[static_cast<std::size_t>(eaf_value)] : kNaN;
      double decoded_se = kNaN;
      if (need_se && se_value >= 0 && se_value < se_count_value && centre_count > 0) {
        const R_xlen_t block = row / block_rows_value;
        const R_xlen_t centre = std::min<R_xlen_t>(block, centre_count - 1);
        decoded_se = se_table[static_cast<std::size_t>(se_value) * eaf_width +
                              static_cast<std::size_t>(eaf_value)] *
                     centre_factors[static_cast<std::size_t>(centre)];
      }
      if (need_z) z_out[row] = decoded_z;
      if (need_se) se_out[row] = decoded_se;
      if (need_eaf) eaf_out[row] = decoded_eaf;
      if (want_beta) beta_out[row] = decoded_z * decoded_se;
      if (want_p) p_out[row] = z_value >= 0 && z_value < z_count_value
        ? p_table[static_cast<std::size_t>(z_value)] : kNaN;
    }
    for (std::size_t i = 0; i < exception_rows.size(); ++i) {
      const int row = exception_rows[i];
      if (row < 0 || static_cast<R_xlen_t>(row) >= n) {
        throw std::runtime_error("exception row is outside the Pcodec bridge");
      }
      const int flags = exception_flags[i];
      if (flags <= 0 || (flags & ~7) != 0 ||
          (i > 0 && exception_rows[i - 1] >= row)) {
        throw std::runtime_error("Pcodec bridge exception records are invalid");
      }
      if ((flags & 1) && need_z) REAL(z)[row] = static_cast<double>(exception_z[i]);
      if ((flags & 2) && need_se) REAL(se)[row] = exception_se[i];
      if ((flags & 4) && need_eaf) REAL(eaf)[row] = static_cast<double>(exception_eaf[i]);
      if (want_beta) REAL(beta)[row] = REAL(z)[row] * REAL(se)[row];
      if (want_p && (flags & 1)) {
        REAL(p)[row] = std::erfc(std::abs(REAL(z)[row]) * kInvSqrtTwo);
      }
    }
  }

  const R_xlen_t output_count = XLENGTH(requested);
  SEXP output = PROTECT(Rf_allocVector(VECSXP, output_count)); ++protect_count;
  SEXP names = PROTECT(Rf_allocVector(STRSXP, output_count)); ++protect_count;
  for (R_xlen_t i = 0; i < output_count; ++i) {
    const std::string name = CHAR(STRING_ELT(requested, i));
    SEXP value = R_NilValue;
    if (name == "chromosome") value = chromosome;
    else if (name == "base_pair_location") value = position;
    else if (name == "effect_allele") value = effect;
    else if (name == "other_allele") value = other;
    else if (name == "z") value = z;
    else if (name == "beta") value = beta;
    else if (name == "standard_error") value = se;
    else if (name == "effect_allele_frequency") value = eaf;
    else if (name == "p_value") value = p;
    else throw std::runtime_error("unknown requested Pcodec bridge column: " + name);
    SET_VECTOR_ELT(output, i, value);
    SET_STRING_ELT(names, i, STRING_ELT(requested, i));
  }
  Rf_setAttrib(output, R_NamesSymbol, names);
  Rf_setAttrib(output, R_ClassSymbol, Rf_mkString("data.frame"));
  SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2)); ++protect_count;
  INTEGER(row_names)[0] = NA_INTEGER;
  INTEGER(row_names)[1] = -static_cast<int>(n);
  Rf_setAttrib(output, R_RowNamesSymbol, row_names);
  UNPROTECT(protect_count);
  return output;
}

extern "C" SEXP compressor_read_pcodec_bridge(
    SEXP directory,
    SEXP requested,
    SEXP row_count,
    SEXP centres,
    SEXP z_min,
    SEXP z_max,
    SEXP z_count,
    SEXP se_count,
    SEXP eaf_count,
    SEXP block_rows,
    SEXP chromosome_lengths) {
  try {
    return compressor_read_pcodec_bridge_impl(
      directory, requested, row_count, centres, z_min, z_max,
      z_count, se_count, eaf_count, block_rows, chromosome_lengths);
  } catch (const std::exception& exception) {
    Rf_error("%s", exception.what());
  } catch (...) {
    Rf_error("unknown native Pcodec bridge error");
  }
  return R_NilValue;
}
