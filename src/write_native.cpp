// Native helpers for the store-writing path that do not depend on Pcodec.
//
// compressor_unique_strings(): one hashing pass that returns the distinct
//   CHARSXPs of a character vector (first-occurrence order) and the 1-based
//   index of every element into them, so an elementwise string transform can
//   be evaluated once per distinct value.
// compressor_serialize_sha256(): the SHA-256 of serialize(object, NULL,
//   version = version) after skipping `skip` leading bytes, computed while
//   streaming the serialization instead of materialising it.  With skip = 14
//   it equals digest::digest(object, algo = "sha256") (binary XDR format).

#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif

#include <array>
#include <atomic>
#include <climits>
#include <cstdio>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

// ---------------------------------------------------------------- SHA-256
struct Sha256 {
  std::uint32_t state[8];
  std::uint64_t bytes;
  unsigned char buffer[64];
  std::size_t used;
};

const std::uint32_t kSha256K[64] = {
  0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
  0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
  0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
  0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
  0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
  0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
  0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
  0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
  0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
  0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
  0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

inline std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

void sha256_init(Sha256* ctx) {
  static const std::uint32_t initial[8] = {
    0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
    0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
  };
  std::memcpy(ctx->state, initial, sizeof(initial));
  ctx->bytes = 0;
  ctx->used = 0;
}

void sha256_block(std::uint32_t* state, const unsigned char* block) {
  std::uint32_t w[64];
  for (int i = 0; i < 16; ++i) {
    w[i] = (static_cast<std::uint32_t>(block[4 * i]) << 24) |
      (static_cast<std::uint32_t>(block[4 * i + 1]) << 16) |
      (static_cast<std::uint32_t>(block[4 * i + 2]) << 8) |
      static_cast<std::uint32_t>(block[4 * i + 3]);
  }
  for (int i = 16; i < 64; ++i) {
    const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
    const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  std::uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
  std::uint32_t e = state[4], f = state[5], g = state[6], h = state[7];
  for (int i = 0; i < 64; ++i) {
    const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
    const std::uint32_t ch = (e & f) ^ (~e & g);
    const std::uint32_t t1 = h + s1 + ch + kSha256K[i] + w[i];
    const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
    const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
    const std::uint32_t t2 = s0 + maj;
    h = g; g = f; f = e; e = d + t1;
    d = c; c = b; b = a; a = t1 + t2;
  }
  state[0] += a; state[1] += b; state[2] += c; state[3] += d;
  state[4] += e; state[5] += f; state[6] += g; state[7] += h;
}

void sha256_update(Sha256* ctx, const unsigned char* data, std::size_t n) {
  ctx->bytes += n;
  if (ctx->used) {
    const std::size_t take = (64 - ctx->used < n) ? 64 - ctx->used : n;
    std::memcpy(ctx->buffer + ctx->used, data, take);
    ctx->used += take;
    data += take;
    n -= take;
    if (ctx->used < 64) return;
    sha256_block(ctx->state, ctx->buffer);
    ctx->used = 0;
  }
  while (n >= 64) {
    sha256_block(ctx->state, data);
    data += 64;
    n -= 64;
  }
  if (n) {
    std::memcpy(ctx->buffer, data, n);
    ctx->used = n;
  }
}

void sha256_hex(Sha256* ctx, char* out /* 65 bytes */) {
  const std::uint64_t bits = ctx->bytes * 8u;
  unsigned char pad[72];
  std::size_t pad_n = (ctx->used < 56) ? 56 - ctx->used : 120 - ctx->used;
  std::memset(pad, 0, sizeof(pad));
  pad[0] = 0x80u;
  for (int i = 0; i < 8; ++i) {
    pad[pad_n + i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
  }
  const std::uint64_t saved = ctx->bytes;
  sha256_update(ctx, pad, pad_n + 8);
  ctx->bytes = saved;
  static const char hex[] = "0123456789abcdef";
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 4; ++j) {
      const unsigned int byte = (ctx->state[i] >> (24 - 8 * j)) & 0xffu;
      out[8 * i + 2 * j] = hex[byte >> 4];
      out[8 * i + 2 * j + 1] = hex[byte & 0xfu];
    }
  }
  out[64] = '\0';
}

// --------------------------------------------------- streaming serialize
// Plain C-style state only: an R error or interrupt during R_Serialize()
// longjmps over this frame, so it must not own anything with a destructor.
struct HashStream {
  Sha256 sha;
  std::uint64_t skip;  // leading serialized bytes still to discard
};

void hash_bytes(HashStream* hs, const unsigned char* data, std::size_t n) {
  if (hs->skip) {
    const std::uint64_t drop = hs->skip < n ? hs->skip : n;
    hs->skip -= drop;
    data += drop;
    n -= static_cast<std::size_t>(drop);
  }
  if (n) sha256_update(&hs->sha, data, n);
}

void hash_out_char(R_outpstream_t stream, int c) {
  const unsigned char byte = static_cast<unsigned char>(c);
  hash_bytes(static_cast<HashStream*>(stream->data), &byte, 1);
}

void hash_out_bytes(R_outpstream_t stream, void* buf, int length) {
  hash_bytes(static_cast<HashStream*>(stream->data),
             static_cast<const unsigned char*>(buf), static_cast<std::size_t>(length));
}

// ------------------------------------------- pipelined streaming serialize
// The main thread runs R_Serialize() and hands 1 MiB chunks to one worker
// thread that hashes them, so serialization and hashing overlap.  The chunk
// queue is bounded, so memory stays a few MiB.  The pipeline lives on the
// heap and is torn down by the R_UnwindProtect() cleanup if serialization
// longjmps (an R error or user interrupt).
struct HashPipeline {
  Sha256 sha;
  std::uint64_t skip = 0;
  std::size_t chunk_bytes = 1u << 20;
  std::size_t max_queued = 8;
  std::vector<unsigned char> current;
  std::deque<std::vector<unsigned char>> queued;
  std::mutex mutex;
  std::condition_variable items;
  std::condition_variable space;
  bool finished = false;
  bool aborted = false;
  bool failed = false;
  std::thread worker;

  void run() {
    while (true) {
      std::vector<unsigned char> chunk;
      {
        std::unique_lock<std::mutex> lock(mutex);
        items.wait(lock, [&] { return aborted || finished || !queued.empty(); });
        if (aborted) return;
        if (queued.empty()) return;  // finished and drained
        chunk.swap(queued.front());
        queued.pop_front();
      }
      space.notify_one();
      sha256_update(&sha, chunk.data(), chunk.size());
    }
  }

  void push_current() {
    std::unique_lock<std::mutex> lock(mutex);
    space.wait(lock, [&] { return aborted || queued.size() < max_queued; });
    queued.emplace_back();
    queued.back().swap(current);
    lock.unlock();
    items.notify_one();
    current.clear();
    current.reserve(chunk_bytes);
  }

  void stop(bool abort_now) {
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (abort_now) aborted = true; else finished = true;
    }
    items.notify_all();
    space.notify_all();
    if (worker.joinable()) worker.join();
  }
};

void pipeline_bytes(HashPipeline* hp, const unsigned char* data, std::size_t n) {
  if (hp->failed) return;
  if (hp->skip) {
    const std::uint64_t drop = hp->skip < n ? hp->skip : n;
    hp->skip -= drop;
    data += drop;
    n -= static_cast<std::size_t>(drop);
  }
  try {
    while (n) {
      const std::size_t room = hp->chunk_bytes - hp->current.size();
      const std::size_t take = room < n ? room : n;
      hp->current.insert(hp->current.end(), data, data + take);
      data += take;
      n -= take;
      if (hp->current.size() == hp->chunk_bytes) hp->push_current();
    }
  } catch (...) {
    hp->failed = true;  // reported after R_Serialize() returns
  }
}

void pipeline_out_char(R_outpstream_t stream, int c) {
  const unsigned char byte = static_cast<unsigned char>(c);
  pipeline_bytes(static_cast<HashPipeline*>(stream->data), &byte, 1);
}

void pipeline_out_bytes(R_outpstream_t stream, void* buf, int length) {
  pipeline_bytes(static_cast<HashPipeline*>(stream->data),
                 static_cast<const unsigned char*>(buf), static_cast<std::size_t>(length));
}

struct PipelineJob {
  SEXP object;
  int version;
  HashPipeline* pipeline;
};

SEXP pipeline_serialize(void* data) {
  PipelineJob* job = static_cast<PipelineJob*>(data);
  struct R_outpstream_st out;
  R_InitOutPStream(&out, static_cast<R_pstream_data_t>(job->pipeline), R_pstream_xdr_format,
                   job->version, pipeline_out_char, pipeline_out_bytes, NULL, R_NilValue);
  R_Serialize(job->object, &out);
  return R_NilValue;
}

void pipeline_cleanup(void* data, Rboolean jump) {
  if (!jump) return;
  HashPipeline* pipeline = static_cast<PipelineJob*>(data)->pipeline;
  pipeline->stop(true);
  delete pipeline;
}

// ------------------------------------------------------- unique strings
inline std::uint64_t pointer_hash(const void* p) {
  std::uint64_t x = reinterpret_cast<std::uintptr_t>(p);
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  return x;
}

}  // namespace

extern "C" SEXP compressor_serialize_sha256(SEXP object, SEXP version, SEXP skip,
                                            SEXP threads) {
  const int serialize_version = Rf_asInteger(version);
  if (serialize_version != 2 && serialize_version != 3) {
    Rf_error("serialization version must be 2 or 3");
  }
  const double skip_bytes = Rf_asReal(skip);
  if (!R_FINITE(skip_bytes) || skip_bytes < 0) Rf_error("skip must be non-negative");
  const int thread_count = Rf_asInteger(threads);
  // The pipeline adds one hashing thread, so it is only used when the caller
  // allows more than one thread.
  HashPipeline* pipeline = (thread_count != NA_INTEGER && thread_count > 1)
    ? new (std::nothrow) HashPipeline() : NULL;
  if (pipeline) {
    sha256_init(&pipeline->sha);
    pipeline->skip = static_cast<std::uint64_t>(skip_bytes);
    bool started = false;
    try {
      pipeline->current.reserve(pipeline->chunk_bytes);
      pipeline->worker = std::thread([pipeline] { pipeline->run(); });
      started = true;
    } catch (...) {
      started = false;
    }
    if (started) {
      PipelineJob job{object, serialize_version, pipeline};
      SEXP cont = PROTECT(R_MakeUnwindCont());
      R_UnwindProtect(pipeline_serialize, &job, pipeline_cleanup, &job, cont);
      UNPROTECT(1);
      bool failed = pipeline->failed;
      if (!failed) {
        try {
          if (!pipeline->current.empty()) pipeline->push_current();
        } catch (...) {
          failed = true;
        }
      }
      pipeline->stop(failed);
      char hex[65];
      if (!failed) sha256_hex(&pipeline->sha, hex);
      delete pipeline;
      if (failed) Rf_error("could not hash the serialized object (out of memory)");
      return Rf_mkString(hex);
    }
    delete pipeline;
  }
  // Single-threaded streaming (threads = 1, or no worker thread available).
  HashStream hs;
  sha256_init(&hs.sha);
  hs.skip = static_cast<std::uint64_t>(skip_bytes);
  struct R_outpstream_st out;
  R_InitOutPStream(&out, static_cast<R_pstream_data_t>(&hs), R_pstream_xdr_format,
                   serialize_version, hash_out_char, hash_out_bytes, NULL, R_NilValue);
  R_Serialize(object, &out);
  char hex[65];
  sha256_hex(&hs.sha, hex);
  return Rf_mkString(hex);
}

extern "C" SEXP compressor_unique_strings(SEXP x) {
  if (TYPEOF(x) != STRSXP) Rf_error("compressor_unique_strings expects a character vector");
  const R_xlen_t n = XLENGTH(x);
  if (n > INT_MAX) Rf_error("character vector is too long");
  SEXP index = PROTECT(Rf_allocVector(INTSXP, n));
  int* idx = INTEGER(index);
  const SEXP* values = STRING_PTR_RO(x);
  std::vector<SEXP> distinct;
  std::size_t capacity = 64;
  std::vector<int> table(capacity, -1);  // slot -> position in `distinct`
  std::size_t mask = capacity - 1;
  SEXP last = NULL;
  int last_id = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const SEXP s = values[i];
    if (s == last) {
      idx[i] = last_id + 1;
      continue;
    }
    std::size_t slot = static_cast<std::size_t>(pointer_hash(s)) & mask;
    int id;
    while (true) {
      const int at = table[slot];
      if (at < 0) {
        id = static_cast<int>(distinct.size());
        distinct.push_back(s);
        table[slot] = id;
        if (distinct.size() * 2 > capacity) {
          capacity *= 2;
          mask = capacity - 1;
          std::vector<int> grown(capacity, -1);
          for (std::size_t k = 0; k < distinct.size(); ++k) {
            std::size_t probe = static_cast<std::size_t>(pointer_hash(distinct[k])) & mask;
            while (grown[probe] >= 0) probe = (probe + 1) & mask;
            grown[probe] = static_cast<int>(k);
          }
          table.swap(grown);
        }
        break;
      }
      if (distinct[static_cast<std::size_t>(at)] == s) {
        id = at;
        break;
      }
      slot = (slot + 1) & mask;
    }
    last = s;
    last_id = id;
    idx[i] = id + 1;
  }
  SEXP unique = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(distinct.size())));
  for (std::size_t k = 0; k < distinct.size(); ++k) {
    SET_STRING_ELT(unique, static_cast<R_xlen_t>(k), distinct[k]);
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, unique);
  SET_VECTOR_ELT(out, 1, index);
  UNPROTECT(3);
  return out;
}

// TRUE when `a` and `b` are character vectors of the same length holding the
// same CHARSXP at every position.  CHARSXPs are cached, so equal pointers mean
// the same bytes and the same encoding mark.
extern "C" SEXP compressor_same_strings(SEXP a, SEXP b) {
  if (TYPEOF(a) != STRSXP || TYPEOF(b) != STRSXP || XLENGTH(a) != XLENGTH(b)) {
    return Rf_ScalarLogical(FALSE);
  }
  const R_xlen_t n = XLENGTH(a);
  const SEXP* pa = STRING_PTR_RO(a);
  const SEXP* pb = STRING_PTR_RO(b);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (pa[i] != pb[i]) return Rf_ScalarLogical(FALSE);
  }
  return Rf_ScalarLogical(TRUE);
}

// ------------------------------------------------- column content hashing
// compressor_column_xxh64(columns, threads): one XXH64 (seed 0) digest per
// column of a list of atomic vectors, computed over a platform-independent
// byte stream (little-endian):
//   logical/integer  int32 per element (NA = INT_MIN)
//   double           IEEE-754 binary64 bits per element, with R's NA as
//                    0x7FF00000000007A2 and every other NaN as
//                    0x7FF8000000000000
//   character        per element uint32 byte length then the bytes of the
//                    CHARSXP; NA is the length 0xFFFFFFFF with no bytes
// Columns are hashed in parallel on up to `threads` threads; each digest is
// independent of the thread count.  Data pointers are taken on the calling
// thread (materialising any ALTREP column) before workers start, and workers
// only read memory.
namespace {

const std::uint64_t kXxhP1 = 11400714785074694791ULL;
const std::uint64_t kXxhP2 = 14029467366897019727ULL;
const std::uint64_t kXxhP3 = 1609587929392839161ULL;
const std::uint64_t kXxhP4 = 9650029242287828579ULL;
const std::uint64_t kXxhP5 = 2870177450012600261ULL;

inline std::uint64_t rotl64(std::uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }

inline std::uint64_t read_le64(const unsigned char* p) {
  std::uint64_t v = 0;
  for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
  return v;
}

inline std::uint32_t read_le32(const unsigned char* p) {
  return static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
    (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
}

inline std::uint64_t xxh_round(std::uint64_t acc, std::uint64_t input) {
  acc += input * kXxhP2;
  acc = rotl64(acc, 31);
  return acc * kXxhP1;
}

inline std::uint64_t xxh_merge(std::uint64_t acc, std::uint64_t val) {
  acc ^= xxh_round(0, val);
  return acc * kXxhP1 + kXxhP4;
}

struct Xxh64 {
  std::uint64_t v1, v2, v3, v4;
  std::uint64_t total;
  unsigned char buffer[32];
  std::size_t used;
};

void xxh64_init(Xxh64* s) {
  s->v1 = kXxhP1 + kXxhP2;
  s->v2 = kXxhP2;
  s->v3 = 0;
  s->v4 = 0 - kXxhP1;
  s->total = 0;
  s->used = 0;
}

inline void xxh64_stripe(Xxh64* s, const unsigned char* p) {
  s->v1 = xxh_round(s->v1, read_le64(p));
  s->v2 = xxh_round(s->v2, read_le64(p + 8));
  s->v3 = xxh_round(s->v3, read_le64(p + 16));
  s->v4 = xxh_round(s->v4, read_le64(p + 24));
}

void xxh64_update(Xxh64* s, const unsigned char* p, std::size_t n) {
  s->total += n;
  if (s->used + n < 32) {
    std::memcpy(s->buffer + s->used, p, n);
    s->used += n;
    return;
  }
  if (s->used) {
    const std::size_t take = 32 - s->used;
    std::memcpy(s->buffer + s->used, p, take);
    xxh64_stripe(s, s->buffer);
    p += take;
    n -= take;
    s->used = 0;
  }
  while (n >= 32) {
    xxh64_stripe(s, p);
    p += 32;
    n -= 32;
  }
  if (n) {
    std::memcpy(s->buffer, p, n);
    s->used = n;
  }
}

std::uint64_t xxh64_digest(const Xxh64* s) {
  std::uint64_t h;
  if (s->total >= 32) {
    h = rotl64(s->v1, 1) + rotl64(s->v2, 7) + rotl64(s->v3, 12) + rotl64(s->v4, 18);
    h = xxh_merge(h, s->v1);
    h = xxh_merge(h, s->v2);
    h = xxh_merge(h, s->v3);
    h = xxh_merge(h, s->v4);
  } else {
    h = s->v3 + kXxhP5;
  }
  h += s->total;
  const unsigned char* p = s->buffer;
  std::size_t n = s->used;
  while (n >= 8) {
    h ^= xxh_round(0, read_le64(p));
    h = rotl64(h, 27) * kXxhP1 + kXxhP4;
    p += 8;
    n -= 8;
  }
  if (n >= 4) {
    h ^= static_cast<std::uint64_t>(read_le32(p)) * kXxhP1;
    h = rotl64(h, 23) * kXxhP2 + kXxhP3;
    p += 4;
    n -= 4;
  }
  while (n) {
    h ^= static_cast<std::uint64_t>(*p) * kXxhP5;
    h = rotl64(h, 11) * kXxhP1;
    ++p;
    --n;
  }
  h ^= h >> 33;
  h *= kXxhP2;
  h ^= h >> 29;
  h *= kXxhP3;
  h ^= h >> 32;
  return h;
}

// Little-endian staging buffer so elements are fed to XXH64 in large runs.
struct LeWriter {
  Xxh64* state;
  unsigned char buf[1 << 16];
  std::size_t used = 0;
  explicit LeWriter(Xxh64* s) : state(s) {}
  inline void room(std::size_t n) {
    if (used + n > sizeof(buf)) flush();
  }
  inline void u32(std::uint32_t v) {
    room(4);
    buf[used] = static_cast<unsigned char>(v);
    buf[used + 1] = static_cast<unsigned char>(v >> 8);
    buf[used + 2] = static_cast<unsigned char>(v >> 16);
    buf[used + 3] = static_cast<unsigned char>(v >> 24);
    used += 4;
  }
  inline void u64(std::uint64_t v) {
    room(8);
    for (int i = 0; i < 8; ++i) buf[used + i] = static_cast<unsigned char>(v >> (8 * i));
    used += 8;
  }
  inline void bytes(const char* p, std::size_t n) {
    if (n > sizeof(buf) / 2) {
      flush();
      xxh64_update(state, reinterpret_cast<const unsigned char*>(p), n);
      return;
    }
    room(n);
    std::memcpy(buf + used, p, n);
    used += n;
  }
  void flush() {
    if (used) xxh64_update(state, buf, used);
    used = 0;
  }
};

struct ColumnJob {
  int kind = 0;  // 1 int32, 2 double, 3 string
  R_xlen_t n = 0;
  const int* ints = NULL;
  const double* reals = NULL;
  const SEXP* strings = NULL;
  SEXP na_string = NULL;
  std::uint64_t digest = 0;
};

void hash_column(ColumnJob* job) {
  Xxh64 state;
  xxh64_init(&state);
  LeWriter out(&state);
  const R_xlen_t n = job->n;
  if (job->kind == 1) {
    for (R_xlen_t i = 0; i < n; ++i) out.u32(static_cast<std::uint32_t>(job->ints[i]));
  } else if (job->kind == 2) {
    const std::uint64_t na_bits = 0x7FF00000000007A2ULL;
    const std::uint64_t nan_bits = 0x7FF8000000000000ULL;
    for (R_xlen_t i = 0; i < n; ++i) {
      const double v = job->reals[i];
      std::uint64_t bits;
      if (v != v) {
        std::uint64_t raw;
        std::memcpy(&raw, &v, sizeof(raw));
        // R_IsNA: a NaN whose low 32-bit word is 1954.
        bits = (static_cast<std::uint32_t>(raw & 0xFFFFFFFFu) == 1954u) ? na_bits : nan_bits;
      } else {
        std::memcpy(&bits, &v, sizeof(bits));
      }
      out.u64(bits);
    }
  } else if (job->kind == 3) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const SEXP s = job->strings[i];
      if (s == job->na_string) {
        out.u32(0xFFFFFFFFu);
        continue;
      }
      const R_len_t len = LENGTH(s);
      out.u32(static_cast<std::uint32_t>(len));
      out.bytes(CHAR(s), static_cast<std::size_t>(len));
    }
  }
  out.flush();
  job->digest = xxh64_digest(&state);
}

}  // namespace

extern "C" SEXP compressor_column_xxh64(SEXP columns, SEXP threads) {
  if (TYPEOF(columns) != VECSXP) Rf_error("columns must be a list");
  const R_xlen_t k = XLENGTH(columns);
  // Validate before any C++ object exists: Rf_error() must not skip a
  // destructor.
  for (R_xlen_t j = 0; j < k; ++j) {
    const int type = TYPEOF(VECTOR_ELT(columns, j));
    if (type != LGLSXP && type != INTSXP && type != REALSXP && type != STRSXP) {
      Rf_error("column %d has an unsupported type for content hashing", (int) j + 1);
    }
  }
  std::vector<ColumnJob> jobs(static_cast<std::size_t>(k));
  for (R_xlen_t j = 0; j < k; ++j) {
    SEXP column = VECTOR_ELT(columns, j);
    ColumnJob& job = jobs[static_cast<std::size_t>(j)];
    job.n = XLENGTH(column);
    switch (TYPEOF(column)) {
      case LGLSXP: job.kind = 1; job.ints = LOGICAL_RO(column); break;
      case INTSXP: job.kind = 1; job.ints = INTEGER_RO(column); break;
      case REALSXP: job.kind = 2; job.reals = REAL_RO(column); break;
      case STRSXP: job.kind = 3; job.strings = STRING_PTR_RO(column); job.na_string = NA_STRING; break;
      default: break;  // rejected above
    }
  }
  int workers = Rf_asInteger(threads);
  if (workers == NA_INTEGER || workers < 1) workers = 1;
  if (static_cast<R_xlen_t>(workers) > k) workers = static_cast<int>(k > 0 ? k : 1);
  bool threaded = false;
  if (workers > 1) {
    std::vector<std::thread> pool;
    std::size_t next = 0;
    std::mutex mutex;
    try {
      for (int w = 0; w < workers; ++w) {
        pool.emplace_back([&] {
          while (true) {
            std::size_t at;
            {
              std::lock_guard<std::mutex> lock(mutex);
              if (next >= jobs.size()) return;
              at = next++;
            }
            hash_column(&jobs[at]);
          }
        });
      }
      threaded = true;
    } catch (...) {
      threaded = false;
    }
    for (std::thread& t : pool) if (t.joinable()) t.join();
    if (!threaded) {
      // Thread creation failed part way: hash every column on this thread.
      for (ColumnJob& job : jobs) job.digest = 0;
    }
  }
  if (!threaded) {
    for (ColumnJob& job : jobs) hash_column(&job);
  }
  SEXP out = PROTECT(Rf_allocVector(STRSXP, k));
  static const char hex[] = "0123456789abcdef";
  for (R_xlen_t j = 0; j < k; ++j) {
    char text[17];
    const std::uint64_t d = jobs[static_cast<std::size_t>(j)].digest;
    for (int i = 0; i < 16; ++i) text[i] = hex[(d >> (60 - 4 * i)) & 0xFu];
    text[16] = '\0';
    SET_STRING_ELT(out, j, Rf_mkChar(text));
  }
  UNPROTECT(1);
  return out;
}

// SHA-256 hex digest of a raw vector (used to combine column digests).
extern "C" SEXP compressor_sha256_raw(SEXP bytes) {
  if (TYPEOF(bytes) != RAWSXP) Rf_error("bytes must be a raw vector");
  Sha256 sha;
  sha256_init(&sha);
  sha256_update(&sha, RAW(bytes), static_cast<std::size_t>(XLENGTH(bytes)));
  char hex[65];
  sha256_hex(&sha, hex);
  return Rf_mkString(hex);
}

// SHA-256 hex digests of whole files, `threads` files at a time.  Equal to
// digest::digest(path, algo = "sha256", file = TRUE).  A file that cannot be
// opened or read yields NA.  The workers use only C stdio and plain buffers,
// never the R API.
namespace {
bool sha256_file(const char* path, char* out /* 65 bytes */) {
  std::FILE* file = std::fopen(path, "rb");
  if (!file) return false;
  Sha256 sha;
  sha256_init(&sha);
  std::vector<unsigned char> buffer(1u << 20);
  bool ok = true;
  while (true) {
    const std::size_t got = std::fread(buffer.data(), 1, buffer.size(), file);
    if (got) sha256_update(&sha, buffer.data(), got);
    if (got < buffer.size()) {
      ok = !std::ferror(file);
      break;
    }
  }
  std::fclose(file);
  if (ok) sha256_hex(&sha, out);
  return ok;
}
}  // namespace

extern "C" SEXP compressor_sha256_files(SEXP paths, SEXP threads) {
  if (TYPEOF(paths) != STRSXP) Rf_error("paths must be a character vector");
  const R_xlen_t k = XLENGTH(paths);
  std::vector<std::string> names(static_cast<std::size_t>(k));
  std::vector<int> present(static_cast<std::size_t>(k), 0);
  for (R_xlen_t i = 0; i < k; ++i) {
    SEXP s = STRING_ELT(paths, i);
    if (s == NA_STRING) continue;
    names[static_cast<std::size_t>(i)] = R_ExpandFileName(Rf_translateChar(s));
    present[static_cast<std::size_t>(i)] = 1;
  }
  std::vector<std::array<char, 65>> digests(static_cast<std::size_t>(k));
  std::vector<int> ok(static_cast<std::size_t>(k), 0);
  int workers = Rf_asInteger(threads);
  if (workers == NA_INTEGER || workers < 1) workers = 1;
  if (static_cast<R_xlen_t>(workers) > k) workers = static_cast<int>(k > 0 ? k : 1);
  std::atomic<std::size_t> next(0);
  auto work = [&]() {
    while (true) {
      const std::size_t at = next.fetch_add(1);
      if (at >= names.size()) return;
      if (present[at]) ok[at] = sha256_file(names[at].c_str(), digests[at].data()) ? 1 : 0;
    }
  };
  bool threaded = false;
  if (workers > 1) {
    std::vector<std::thread> pool;
    try {
      for (int w = 0; w < workers; ++w) pool.emplace_back(work);
      threaded = true;
    } catch (...) {
      threaded = false;
    }
    for (std::thread& t : pool) if (t.joinable()) t.join();
  }
  // Serial path, or the remaining files if thread creation failed part way.
  if (!threaded) work();
  SEXP out = PROTECT(Rf_allocVector(STRSXP, k));
  for (R_xlen_t i = 0; i < k; ++i) {
    const std::size_t at = static_cast<std::size_t>(i);
    SET_STRING_ELT(out, i, ok[at] ? Rf_mkChar(digests[at].data()) : NA_STRING);
  }
  UNPROTECT(1);
  return out;
}
