#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP compressor_decode_native(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_read_pcodec_bridge(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP);
extern SEXP compressor_pcodec_native_available(void);
extern SEXP compressor_pcodec_compress_u8(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_compress_u16(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_compress_u32(SEXP, SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u8(SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u16(SEXP, SEXP);
extern SEXP compressor_pcodec_decompress_u32(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"compressor_decode_native", (DL_FUNC) &compressor_decode_native, 20},
    {"compressor_read_pcodec_bridge", (DL_FUNC) &compressor_read_pcodec_bridge, 11},
    {"compressor_pcodec_native_available", (DL_FUNC) &compressor_pcodec_native_available, 0},
    {"compressor_pcodec_compress_u8", (DL_FUNC) &compressor_pcodec_compress_u8, 3},
    {"compressor_pcodec_compress_u16", (DL_FUNC) &compressor_pcodec_compress_u16, 3},
    {"compressor_pcodec_compress_u32", (DL_FUNC) &compressor_pcodec_compress_u32, 3},
    {"compressor_pcodec_decompress_u8", (DL_FUNC) &compressor_pcodec_decompress_u8, 2},
    {"compressor_pcodec_decompress_u16", (DL_FUNC) &compressor_pcodec_decompress_u16, 2},
    {"compressor_pcodec_decompress_u32", (DL_FUNC) &compressor_pcodec_decompress_u32, 2},
    {NULL, NULL, 0}
};

void R_init_CompreSSoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
