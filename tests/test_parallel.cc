#include <atomic>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../src/parallel.hh"
#include "check.hh"

namespace mxi {

namespace {

//: A sum over a million terms in fixed blocks, combined in block order.
double blocked_sum() {
  const std::size_t n = 1000000, blocks = kParallelBlocks;
  std::vector<double> part(blocks, 0.0);
  for_each_block(blocks, [&](std::size_t b) {
    double s = 0.0;
    for (std::size_t i = block_begin(n, blocks, b);
         i < block_begin(n, blocks, b + 1); ++i)
      s += 1.0 / (1.0 + static_cast<double>(i));
    part[b] = s;
  });
  double total = 0.0;
  for (double p : part)
    total += p;
  return total;
}

} // namespace

TEST(the_pool_gives_the_same_sum_every_call_and_on_any_count) {
  // Called hundreds of times, as the error model's search calls it: the same
  // bits each time, and the same on one thread as on eight.
  set_parallel_threads(1);
  const double one = blocked_sum();
  set_parallel_threads(8);
  for (int k = 0; k < 300; ++k)
    if (blocked_sum() != one) {
      set_parallel_threads(0);
      check::fail("call " + std::to_string(k) + " differed");
    }
  set_parallel_threads(3); // the pool is made again for a new count
  check::is_true(blocked_sum() == one,
                 "and on three threads, after a change of count");
  set_parallel_threads(0);
}

TEST(parallel_work_started_inside_a_block_runs_serially) {
  // A block that itself calls for_each_block must not wait on the pool it is
  // running in: it runs its inner work serially, and everything is done once.
  set_parallel_threads(4);
  std::atomic<int> count{0};
  for_each_block(8, [&](std::size_t) {
    for_each_block(8, [&](std::size_t) { count.fetch_add(1); });
  });
  set_parallel_threads(0);
  check::equal(static_cast<long long>(count.load()), 64,
               "eight inner blocks in each of eight");
}

TEST(an_exception_in_a_block_reaches_the_caller_and_the_pool_carries_on) {
  set_parallel_threads(4);
  bool caught = false;
  try {
    for_each_block(64, [&](std::size_t b) {
      if (b == 17)
        throw std::runtime_error("block 17");
    });
  } catch (const std::runtime_error &e) {
    caught = std::string(e.what()) == "block 17";
  }
  std::atomic<int> after{0};
  for_each_block(64, [&](std::size_t) { after.fetch_add(1); });
  set_parallel_threads(0);
  check::is_true(caught, "the caller gets the exception");
  check::equal(static_cast<long long>(after.load()), 64,
               "and the pool works after it");
}

} // namespace mxi
