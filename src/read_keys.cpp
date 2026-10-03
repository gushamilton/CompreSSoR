// Native parser for canonical "chromosome:position:REF:ALT" variant keys.
//
// Accepts only the strict canonical spelling (optional case-insensitive "chr"
// prefix; chromosome 1-22, X, Y or the 23/24 aliases; a decimal position on
// the chromosome; two distinct single A/C/G/T alleles in either case). Any
// other input -- NA, extra or missing fields, whitespace inside the key,
// out-of-range positions, invalid alleles -- returns NULL so the caller falls
// back to the R parser, which keeps the established validation and error
// messages. For accepted input the result equals the R path exactly: sorted
// unique numeric identity codes (global position * 16 + substitution) plus
// the decoded positions and substitution codes.
#include <algorithm>
#include <cstring>
#include <vector>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

namespace {

inline int key_base_code(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

// trimws() whitespace: [ \t\r\n]
inline bool key_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

// Chromosome index into the build table order 1..22, X, Y, or -1.
inline int key_chromosome(const char* p, std::size_t len) {
  if (len == 1) {
    const char c = p[0];
    if (c >= '1' && c <= '9') return c - '1';
    if (c == 'X' || c == 'x') return 22;
    if (c == 'Y' || c == 'y') return 23;
    return -1;
  }
  if (len == 2) {
    if (p[0] < '1' || p[0] > '2' || p[1] < '0' || p[1] > '9') return -1;
    const int value = (p[0] - '0') * 10 + (p[1] - '0');
    if (value <= 22) return value - 1;
    if (value == 23) return 22;
    if (value == 24) return 23;
  }
  return -1;
}

// Strictly parse one key into its identity code; false means "not strict".
inline bool key_parse_one(const char* s, std::size_t len, bool trim,
                          const double* offsets, const double* lengths,
                          double* code) {
  const char* p = s;
  const char* end = s + len;
  if (trim) {
    while (p < end && key_space(*p)) ++p;
    while (end > p && key_space(end[-1])) --end;
  }
  if (end - p >= 3 && (p[0] == 'c' || p[0] == 'C') && (p[1] == 'h' || p[1] == 'H') &&
      (p[2] == 'r' || p[2] == 'R')) {
    p += 3;
  }
  const char* c0 = p;
  while (p < end && *p != ':') ++p;
  if (p == end) return false;
  const int chromosome = key_chromosome(c0, static_cast<std::size_t>(p - c0));
  if (chromosome < 0) return false;
  ++p;
  double position = 0;
  int digits = 0;
  while (p < end && *p >= '0' && *p <= '9') {
    if (++digits > 10) return false;
    position = position * 10 + (*p - '0');
    ++p;
  }
  if (!digits || p == end || *p != ':') return false;
  if (position < 1 || position > lengths[chromosome]) return false;
  ++p;
  // exactly "R:A" remains
  if (end - p != 3 || p[1] != ':') return false;
  const int ref = key_base_code(p[0]);
  const int alt = key_base_code(p[2]);
  if (ref < 0 || alt < 0 || ref == alt) return false;
  *code = (offsets[chromosome] + position - 1) * 16 + static_cast<double>(ref * 4 + alt);
  return true;
}

}  // namespace

extern "C" SEXP compressor_parse_variant_keys(SEXP keys, SEXP offsets_sexp,
                                              SEXP lengths_sexp, SEXP trim_sexp) {
  if (TYPEOF(keys) != STRSXP || TYPEOF(offsets_sexp) != REALSXP ||
      TYPEOF(lengths_sexp) != REALSXP || XLENGTH(offsets_sexp) != 24 ||
      XLENGTH(lengths_sexp) != 24 || TYPEOF(trim_sexp) != LGLSXP ||
      XLENGTH(trim_sexp) != 1) {
    Rf_error("malformed arguments to the native variant key parser");
  }
  const R_xlen_t n = XLENGTH(keys);
  const double* offsets = REAL(offsets_sexp);
  const double* lengths = REAL(lengths_sexp);
  const bool trim = LOGICAL(trim_sexp)[0] == 1;
  std::vector<double> codes(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(keys, i);
    if (s == NA_STRING) return R_NilValue;
    if (!key_parse_one(CHAR(s), static_cast<std::size_t>(LENGTH(s)), trim, offsets,
                       lengths, &codes[static_cast<std::size_t>(i)])) {
      return R_NilValue;
    }
  }
  if (!std::is_sorted(codes.begin(), codes.end())) std::sort(codes.begin(), codes.end());
  codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
  const R_xlen_t m = static_cast<R_xlen_t>(codes.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP position = Rf_allocVector(REALSXP, m);
  SET_VECTOR_ELT(out, 0, position);
  SEXP substitution = Rf_allocVector(INTSXP, m);
  SET_VECTOR_ELT(out, 1, substitution);
  SEXP code_out = Rf_allocVector(REALSXP, m);
  SET_VECTOR_ELT(out, 2, code_out);
  double* pp = REAL(position);
  int* sp = INTEGER(substitution);
  double* cp = REAL(code_out);
  for (R_xlen_t i = 0; i < m; ++i) {
    const double c = codes[static_cast<std::size_t>(i)];
    const double pos = static_cast<double>(static_cast<long long>(c) >> 4);
    pp[i] = pos;
    sp[i] = static_cast<int>(static_cast<long long>(c) & 15LL);
    cp[i] = c;
  }
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("position"));
  SET_STRING_ELT(names, 1, Rf_mkChar("substitution"));
  SET_STRING_ELT(names, 2, Rf_mkChar("codes"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}
