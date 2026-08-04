#ifndef OPT_PCODEC_NATIVE_H
#define OPT_PCODEC_NATIVE_H
#include <stddef.h>
extern "C" {
typedef enum CompressorPcoError { COMPRESSOR_PCO_SUCCESS, COMPRESSOR_PCO_INVALID_TYPE, COMPRESSOR_PCO_COMPRESSION_ERROR, COMPRESSOR_PCO_DECOMPRESSION_ERROR } CompressorPcoError;
CompressorPcoError compressor_pco_decompress_into(const void*, size_t, unsigned char, void*, size_t, size_t*);
}
#endif
