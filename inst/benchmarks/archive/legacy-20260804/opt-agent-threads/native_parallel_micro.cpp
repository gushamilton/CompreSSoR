#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

#include "pcodec_native.h"

struct Block {
  std::string stream;
  std::string file;
  unsigned char dtype;
  std::size_t offset, length, values, row_start, row_stop;
  std::uint64_t first_position, last_position;
  int block_id;
};

struct Result { std::uint64_t checksum = 0; std::size_t values = 0; };

static std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> x; std::stringstream ss(s); std::string v;
  while (std::getline(ss, v, '\t')) x.push_back(v);
  return x;
}

static std::vector<Block> read_manifest(const std::string& path) {
  std::ifstream in(path); if (!in) throw std::runtime_error("manifest open failed");
  std::string line; std::getline(in, line); std::vector<Block> out;
  while (std::getline(in, line)) {
    if (line.empty()) continue; auto x = split(line); if (x.size() != 11) continue;
    Block b{x[0],x[1],static_cast<unsigned char>(std::stoul(x[2])),
      static_cast<std::size_t>(std::stoull(x[3])),static_cast<std::size_t>(std::stoull(x[4])),
      static_cast<std::size_t>(std::stoull(x[5])),static_cast<std::size_t>(std::stoull(x[6])),
      static_cast<std::size_t>(std::stoull(x[7])),static_cast<std::uint64_t>(std::stoull(x[8])),
      static_cast<std::uint64_t>(std::stoull(x[9])),std::stoi(x[10])};
    out.push_back(std::move(b));
  }
  return out;
}

static std::uint64_t checksum(const void* p, std::size_t n, unsigned char dtype) {
  std::uint64_t s = 0;
  if (dtype == 1) { auto* x = static_cast<const std::uint32_t*>(p); for (std::size_t i=0;i<n;++i) s += x[i]; }
  else if (dtype == 7) { auto* x = static_cast<const std::uint16_t*>(p); for (std::size_t i=0;i<n;++i) s += x[i]; }
  else { auto* x = static_cast<const std::uint8_t*>(p); for (std::size_t i=0;i<n;++i) s += x[i]; }
  return s;
}

static Result decode(const Block& b, bool compact) {
  int fd = ::open(b.file.c_str(), O_RDONLY); if (fd < 0) throw std::runtime_error("data open failed");
  std::vector<std::uint8_t> compressed(b.length);
  std::size_t got = 0;
  while (got < b.length) {
    ssize_t n = ::pread(fd, compressed.data()+got, b.length-got, static_cast<off_t>(b.offset+got));
    if (n <= 0) { ::close(fd); throw std::runtime_error("pread failed"); } got += static_cast<std::size_t>(n);
  }
  ::close(fd);
  std::vector<std::uint8_t> output(b.values * (b.dtype == 1 ? 4 : b.dtype == 7 ? 2 : 1));
  std::size_t written = 0;
  auto status = compressor_pco_decompress_into(compressed.data(), compressed.size(), b.dtype,
                                                output.data(), b.values, &written);
  if (status != COMPRESSOR_PCO_SUCCESS || written != b.values) {
    std::ostringstream msg;
    msg << "decode failed stream=" << b.stream << " block=" << b.block_id
        << " offset=" << b.offset << " length=" << b.length
        << " values=" << b.values << " status=" << static_cast<int>(status)
        << " written=" << written;
    throw std::runtime_error(msg.str());
  }
  Result r; r.values = written;
  // Compact mode still fully decodes, but reports only the compact key streams' work.
  r.checksum = checksum(output.data(), written, b.dtype) + (compact ? 17 : 0);
  return r;
}

static std::vector<Block> select_blocks(const std::vector<Block>& all, const std::string& selection) {
  if (selection == "full") return all;
  std::vector<Block> out;
  for (const auto& b : all) {
    bool keep = false;
    if (selection == "region") {
      keep = b.first_position <= 1000000ULL && b.last_position >= 1ULL;
    } else if (selection == "sparse") {
      // 1,000 uniformly spread row targets; include their containing key/value block.
      std::size_t lo = (b.stream == "position" || b.stream == "substitution") ? 8192 : 65536;
      for (std::size_t row = 0; row < 1124344; row += 1124) {
        if (b.row_start <= row && row < b.row_stop) { keep = true; break; }
      }
    }
    if (keep) out.push_back(b);
  }
  return out;
}

static Result run_tasks(std::vector<std::vector<Block>> tasks, int threads) {
  std::atomic<std::size_t> next{0}; std::mutex err_mu; std::string error;
  std::vector<Result> partial(static_cast<std::size_t>(threads)); std::vector<std::thread> pool;
  for (int t=0; t<threads; ++t) pool.emplace_back([&,t] {
    while (true) {
      auto i = next.fetch_add(1); if (i >= tasks.size()) break;
      try { for (const auto& b : tasks[i]) { auto r = decode(b, b.stream == "position" || b.stream == "substitution"); partial[t].checksum += r.checksum; partial[t].values += r.values; } }
      catch (const std::exception& e) { std::lock_guard<std::mutex> g(err_mu); error = e.what(); break; }
    }
  });
  for (auto& t : pool) t.join(); if (!error.empty()) throw std::runtime_error(error);
  Result r; for (const auto& x : partial) { r.checksum += x.checksum; r.values += x.values; } return r;
}

int main(int argc, char** argv) {
  if (argc < 5) { std::cerr << "usage: micro MANIFEST MODE THREADS SELECTION\n"; return 2; }
  const auto all = read_manifest(argv[1]); std::string mode=argv[2]; int threads=std::stoi(argv[3]);
  auto selected = select_blocks(all, argv[4]);
  if (mode == "compact") {
    selected.erase(std::remove_if(selected.begin(), selected.end(), [](const Block& b) {
      return b.stream != "position" && b.stream != "substitution";
    }), selected.end());
  }
  std::vector<std::vector<Block>> tasks;
  if (mode == "stream") { for (const auto& b : selected) tasks.push_back({b}); }
  else if (mode == "block" || mode == "compact") {
    std::vector<int> ids; for (const auto& b:selected) if (std::find(ids.begin(),ids.end(),b.block_id)==ids.end()) ids.push_back(b.block_id);
    for (int id:ids) { std::vector<Block> x; for (const auto& b:selected) if (b.block_id==id) x.push_back(b); tasks.push_back(std::move(x)); }
  } else { for (const auto& b:selected) tasks.push_back({b}); threads=1; }
  auto start=std::chrono::steady_clock::now(); auto r=run_tasks(std::move(tasks),threads);
  double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::cout << mode << '\t' << threads << '\t' << argv[4] << '\t' << std::fixed << sec << '\t' << r.values << '\t' << r.checksum << '\n';
}
