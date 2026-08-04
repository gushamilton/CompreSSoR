#ifndef COMPRESSOR_PCODEC_NATIVE_H
#define COMPRESSOR_PCODEC_NATIVE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum CompressorPcoError {
  COMPRESSOR_PCO_SUCCESS,
  COMPRESSOR_PCO_INVALID_TYPE,
  COMPRESSOR_PCO_COMPRESSION_ERROR,
  COMPRESSOR_PCO_DECOMPRESSION_ERROR
} CompressorPcoError;

typedef struct CompressorPcoChunkConfig {
  unsigned int compression_level;
  size_t max_page_n;
} CompressorPcoChunkConfig;

size_t compressor_pco_guarantee_file_size(size_t n, unsigned char dtype);
CompressorPcoError compressor_pco_compress_into(
    const void* nums, size_t n, unsigned char dtype,
    const CompressorPcoChunkConfig* config, void* dst, size_t dst_cap,
    size_t* n_written);
CompressorPcoError compressor_pco_decompress_into(
    const void* compressed, size_t compressed_len, unsigned char dtype,
    void* dst, size_t dst_cap, size_t* n_written);

#ifdef __cplusplus
}
#endif

#endif
