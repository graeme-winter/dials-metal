#include "parallel.hh"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace mxi {

namespace {

std::size_t hardware() {
  const unsigned n = std::thread::hardware_concurrency();
  return n > 0 ? n : 1;
}

std::atomic<std::size_t> g_threads{0};

} // namespace

void set_parallel_threads(std::size_t threads) { g_threads.store(threads); }

std::size_t parallel_threads() {
  const std::size_t n = g_threads.load();
  return n > 0 ? n : hardware();
}

void for_each_block(std::size_t blocks,
                    const std::function<void(std::size_t)> &body) {
  const std::size_t threads = std::min(parallel_threads(), blocks);
  if (threads <= 1) {
    for (std::size_t b = 0; b < blocks; ++b)
      body(b);
    return;
  }
  std::atomic<std::size_t> next{0};
  std::exception_ptr failure;
  std::mutex failure_lock;
  const auto work = [&] {
    try {
      for (std::size_t b = next.fetch_add(1); b < blocks; b = next.fetch_add(1))
        body(b);
    } catch (...) {
      const std::lock_guard<std::mutex> lock(failure_lock);
      if (!failure)
        failure = std::current_exception();
      next.store(blocks); // the rest stop taking blocks
    }
  };
  std::vector<std::thread> pool;
  for (std::size_t t = 1; t < threads; ++t)
    pool.emplace_back(work);
  work();
  for (std::thread &t : pool)
    t.join();
  if (failure)
    std::rethrow_exception(failure);
}

void for_each_index(std::size_t n,
                    const std::function<void(std::size_t)> &body) {
  const std::size_t blocks = std::min<std::size_t>(n, kParallelBlocks);
  for_each_block(blocks, [&](std::size_t b) {
    for (std::size_t i = block_begin(n, blocks, b);
         i < block_begin(n, blocks, b + 1); ++i)
      body(i);
  });
}

} // namespace mxi
