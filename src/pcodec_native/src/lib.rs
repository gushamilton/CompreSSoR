#![allow(clippy::missing_safety_doc)]

use std::ptr;

use libc::{c_uchar, c_int, c_uint, c_void, size_t};
use pco::data_types::{Number, NumberType};
use pco::standalone::guarantee;
use pco::{match_number_enum, ChunkConfig, PagingSpec};

#[repr(C)]
pub enum PcoError {
    PcoSuccess,
    PcoInvalidType,
    PcoCompressionError,
    PcoDecompressionError,
}

#[repr(C)]
pub struct PcoChunkConfig {
    pub compression_level: c_uint,
    pub max_page_n: size_t,
}

impl Default for PcoChunkConfig {
    fn default() -> Self {
        Self { compression_level: 8, max_page_n: 0 }
    }
}

impl PcoChunkConfig {
    fn to_chunk_config(&self) -> ChunkConfig {
        let paging_spec = if self.max_page_n == 0 {
            PagingSpec::default()
        } else {
            PagingSpec::EqualPagesUpTo(self.max_page_n)
        };
        ChunkConfig::default()
            .with_compression_level(self.compression_level as usize)
            .with_paging_spec(paging_spec)
            .with_enable_8_bit(true)
    }
}

fn guarantee_file_size<T: Number>(n: size_t, paging_spec: &PagingSpec) -> size_t {
    guarantee::file_size::<T::L>(n, paging_spec).unwrap_or(0)
}

fn compress_into<T: Number>(
    nums: *const c_void,
    n: size_t,
    config: &ChunkConfig,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if n_written.is_null() {
        return PcoError::PcoCompressionError;
    }
    if n == 0 {
        unsafe { *n_written = 0; }
        return PcoError::PcoSuccess;
    }
    if nums.is_null() || dst.is_null() || n_written.is_null() {
        return PcoError::PcoCompressionError;
    }
    let source = unsafe { std::slice::from_raw_parts(nums as *const T, n) };
    let destination = unsafe { std::slice::from_raw_parts_mut(dst as *mut u8, dst_cap) };
    let original_len = destination.len();
    match pco::standalone::simple_compress_into::<T, _>(source, config, destination) {
        Err(_) => PcoError::PcoCompressionError,
        Ok(remaining) => {
            unsafe { *n_written = original_len - remaining.len(); }
            PcoError::PcoSuccess
        }
    }
}

fn decompress_into<T: Number>(
    compressed: *const c_void,
    compressed_len: size_t,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if compressed.is_null() || n_written.is_null() || (dst.is_null() && dst_cap != 0) {
        return PcoError::PcoDecompressionError;
    }
    let source = unsafe { std::slice::from_raw_parts(compressed as *const u8, compressed_len) };
    match pco::standalone::simple_decompress::<T>(source) {
        Err(_) => PcoError::PcoDecompressionError,
        Ok(values) => {
            if values.len() > dst_cap { return PcoError::PcoDecompressionError; }
            unsafe {
                ptr::copy_nonoverlapping(values.as_ptr(), dst as *mut T, values.len());
                *n_written = values.len();
            }
            PcoError::PcoSuccess
        }
    }
}

#[no_mangle]
pub extern "C" fn compressor_pco_guarantee_file_size(n: size_t, dtype: c_uchar) -> size_t {
    let Some(dtype_enum) = NumberType::from_descriminant(dtype) else { return 0; };
    let paging_spec = PagingSpec::default();
    match_number_enum!(dtype_enum, NumberType<T> => { guarantee_file_size::<T>(n, &paging_spec) })
}

#[no_mangle]
pub unsafe extern "C" fn compressor_pco_compress_into(
    nums: *const c_void,
    n: size_t,
    dtype: c_uchar,
    config: *const PcoChunkConfig,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if n_written.is_null() {
        return PcoError::PcoCompressionError;
    }
    let Some(dtype_enum) = NumberType::from_descriminant(dtype) else { return PcoError::PcoInvalidType; };
    let chunk_config = if config.is_null() {
        PcoChunkConfig::default().to_chunk_config()
    } else {
        (*config).to_chunk_config()
    };
    match_number_enum!(dtype_enum, NumberType<T> => {
        compress_into::<T>(nums, n, &chunk_config, dst, dst_cap, n_written)
    })
}

#[no_mangle]
pub extern "C" fn compressor_pco_decompress_into(
    compressed: *const c_void,
    compressed_len: size_t,
    dtype: c_uchar,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if n_written.is_null() {
        return PcoError::PcoDecompressionError;
    }
    let Some(dtype_enum) = NumberType::from_descriminant(dtype) else { return PcoError::PcoInvalidType; };
    match_number_enum!(dtype_enum, NumberType<T> => {
        decompress_into::<T>(compressed, compressed_len, dst, dst_cap, n_written)
    })
}

#[no_mangle]
pub unsafe extern "C" fn compressor_zstd_compress_into(
    input: *const c_void,
    input_len: size_t,
    level: c_int,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if n_written.is_null() || (input.is_null() && input_len != 0) ||
        (dst.is_null() && dst_cap != 0) {
        return PcoError::PcoCompressionError;
    }
    let source = unsafe { std::slice::from_raw_parts(input as *const u8, input_len) };
    let compressed = match zstd::bulk::compress(source, level) {
        Ok(value) => value,
        Err(_) => return PcoError::PcoCompressionError,
    };
    if compressed.len() > dst_cap {
        return PcoError::PcoCompressionError;
    }
    unsafe {
        ptr::copy_nonoverlapping(compressed.as_ptr(), dst as *mut u8, compressed.len());
        *n_written = compressed.len();
    }
    PcoError::PcoSuccess
}

#[no_mangle]
pub unsafe extern "C" fn compressor_zstd_decompress_into(
    input: *const c_void,
    input_len: size_t,
    dst: *mut c_void,
    dst_cap: size_t,
    n_written: *mut size_t,
) -> PcoError {
    if n_written.is_null() || (input.is_null() && input_len != 0) ||
        (dst.is_null() && dst_cap != 0) {
        return PcoError::PcoDecompressionError;
    }
    let source = unsafe { std::slice::from_raw_parts(input as *const u8, input_len) };
    let decompressed = match zstd::bulk::decompress(source, dst_cap) {
        Ok(value) => value,
        Err(_) => return PcoError::PcoDecompressionError,
    };
    if decompressed.len() > dst_cap {
        return PcoError::PcoDecompressionError;
    }
    unsafe {
        ptr::copy_nonoverlapping(decompressed.as_ptr(), dst as *mut u8, decompressed.len());
        *n_written = decompressed.len();
    }
    PcoError::PcoSuccess
}
