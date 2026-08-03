#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP compressor_decode_native(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compressor_read_pcodec_bridge(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"compressor_decode_native", (DL_FUNC) &compressor_decode_native, 20},
    {"compressor_read_pcodec_bridge", (DL_FUNC) &compressor_read_pcodec_bridge, 10},
    {NULL, NULL, 0}
};

void R_init_CompreSSoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
