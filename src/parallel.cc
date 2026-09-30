#include "parallel.hh"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <memory>
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

// Set in a pool worker, so that parallel work started from inside a block runs
// serially rather than waiting on the pool it is itself running in.
thread_local bool t_in_pool = false;

//: Workers made once and kept, asleep between jobs. They were made and joined
//: on every call: the error model's golden-section search calls for_each_block
//: dozens of times a round, and on one core sixteen threads a call cost 0.85 s
//: of the error model's 1.46 -- none of it arithmetic -- and on sixteen cores
//: kept it from getting faster. Which thread runs a block still cannot change
//: a result: the blocks and their order are the caller's.
class Pool {
public:
  explicit Pool(std::size_t workers) {
    for (std::size_t w = 0; w < workers; ++w)
      threads_.emplace_back([this] { work(); });
  }
  Pool(const Pool &) = delete;
  Pool &operator=(const Pool &) = delete;
  ~Pool() {
    {
      const std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    wake_.notify_all();
    for (std::thread &t : threads_)
      t.join();
  }
  std::size_t workers() const { return threads_.size(); }

  void run(std::size_t blocks, const std::function<void(std::size_t)> &body) {
    const std::lock_guard<std::mutex> one_job(job_mutex_); // one job at a time
    {
      const std::lock_guard<std::mutex> lock(mutex_);
      body_ = &body;
      blocks_ = blocks;
      next_.store(0);
      busy_ = threads_.size();
      failure_ = nullptr;
      ++generation_;
    }
    wake_.notify_all();
    take_blocks(); // the caller works too
    std::unique_lock<std::mutex> lock(mutex_);
    done_.wait(lock, [this] { return busy_ == 0; });
    body_ = nullptr;
    if (failure_)
      std::rethrow_exception(failure_);
  }

private:
  void take_blocks() {
    try {
      for (std::size_t b = next_.fetch_add(1); b < blocks_;
           b = next_.fetch_add(1))
        (*body_)(b);
    } catch (...) {
      const std::lock_guard<std::mutex> lock(mutex_);
      if (!failure_)
        failure_ = std::current_exception();
      next_.store(blocks_); // the rest stop taking blocks
    }
  }
  void work() {
    t_in_pool = true;
    std::uint64_t seen = 0;
    for (;;) {
      {
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] { return stop_ || generation_ != seen; });
        if (stop_)
          return;
        seen = generation_;
      }
      take_blocks();
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        --busy_;
      }
      done_.notify_one();
    }
  }

  std::vector<std::thread> threads_;
  std::mutex mutex_, job_mutex_;
  std::condition_variable wake_, done_;
  const std::function<void(std::size_t)> *body_ = nullptr;
  std::size_t blocks_ = 0, busy_ = 0;
  std::atomic<std::size_t> next_{0};
  std::uint64_t generation_ = 0;
  bool stop_ = false;
  std::exception_ptr failure_;
};

std::mutex g_pool_mutex;
std::unique_ptr<Pool> g_pool;

//: The pool for the thread count now asked for, made or remade as it changes.
Pool &pool(std::size_t threads) {
  const std::lock_guard<std::mutex> lock(g_pool_mutex);
  if (!g_pool || g_pool->workers() != threads - 1)
    g_pool = std::make_unique<Pool>(threads - 1);
  return *g_pool;
}

} // namespace

void set_parallel_threads(std::size_t threads) { g_threads.store(threads); }

std::size_t parallel_threads() {
  const std::size_t n = g_threads.load();
  return n > 0 ? n : hardware();
}

void for_each_block(std::size_t blocks,
                    const std::function<void(std::size_t)> &body) {
  const std::size_t threads = std::min(parallel_threads(), blocks);
  if (threads <= 1 || t_in_pool) {
    for (std::size_t b = 0; b < blocks; ++b)
      body(b);
    return;
  }
  pool(parallel_threads()).run(blocks, body);
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
