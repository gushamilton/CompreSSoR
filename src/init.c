#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP compressor_decode_native(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP);
extern SEXP compressor_read_pcodec_native_codes(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_read_pcodec_native_select(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_key_columns(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_validate_native_codes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_native_available(void);
extern SEXP compressor_pcodec_compress_u8(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_compress_u16(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_compress_u32(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u8(SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u16(SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u32(SEXP, SEXP);
extern SEXP compressor_zstd_compress(SEXP, SEXP);
extern SEXP compressor_zstd_decompress(SEXP, SEXP);
extern SEXP compressor_unique_strings(SEXP);
extern SEXP compressor_same_strings(SEXP, SEXP);
extern SEXP compressor_column_xxh64(SEXP, SEXP);
extern SEXP compressor_sha256_raw(SEXP);
extern SEXP compressor_sha256_files(SEXP, SEXP);
extern SEXP compressor_files_equal(SEXP, SEXP, SEXP);
extern SEXP compressor_serialize_sha256(SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_compress_blocks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_exception_blocks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_parse_variant_keys(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"compressor_decode_native", (DL_FUNC) &compressor_decode_native, 22},
    {"compressor_read_pcodec_native_codes", (DL_FUNC) &compressor_read_pcodec_native_codes, 11},
    {"compressor_read_pcodec_native_select", (DL_FUNC) &compressor_read_pcodec_native_select, 15},
    {"compressor_pcodec_key_columns", (DL_FUNC) &compressor_pcodec_key_columns, 5},
    {"compressor_validate_native_codes", (DL_FUNC) &compressor_validate_native_codes, 6},
    {"compressor_pcodec_native_available", (DL_FUNC) &compressor_pcodec_native_available, 0},
    {"compressor_pcodec_compress_u8", (DL_FUNC) &compressor_pcodec_compress_u8, 3},
    {"compressor_pcodec_compress_u16", (DL_FUNC) &compressor_pcodec_compress_u16, 3},
    {"compressor_pcodec_compress_u32", (DL_FUNC) &compressor_pcodec_compress_u32, 3},
    {"compressor_pcodec_decompress_u8", (DL_FUNC) &compressor_pcodec_decompress_u8, 2},
    {"compressor_pcodec_decompress_u16", (DL_FUNC) &compressor_pcodec_decompress_u16, 2},
    {"compressor_pcodec_decompress_u32", (DL_FUNC) &compressor_pcodec_decompress_u32, 2},
    {"compressor_zstd_compress", (DL_FUNC) &compressor_zstd_compress, 2},
    {"compressor_zstd_decompress", (DL_FUNC) &compressor_zstd_decompress, 2},
    {"compressor_unique_strings", (DL_FUNC) &compressor_unique_strings, 1},
    {"compressor_same_strings", (DL_FUNC) &compressor_same_strings, 2},
    {"compressor_column_xxh64", (DL_FUNC) &compressor_column_xxh64, 2},
    {"compressor_sha256_raw", (DL_FUNC) &compressor_sha256_raw, 1},
    {"compressor_sha256_files", (DL_FUNC) &compressor_sha256_files, 2},
    {"compressor_files_equal", (DL_FUNC) &compressor_files_equal, 3},
    {"compressor_serialize_sha256", (DL_FUNC) &compressor_serialize_sha256, 4},
    {"compressor_pcodec_compress_blocks", (DL_FUNC) &compressor_pcodec_compress_blocks, 6},
    {"compressor_exception_blocks", (DL_FUNC) &compressor_exception_blocks, 8},
    {"compressor_parse_variant_keys", (DL_FUNC) &compressor_parse_variant_keys, 4},
    {NULL, NULL, 0}
};

void R_init_CompreSSoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
