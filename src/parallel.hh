#pragma once

// Work in parallel that gives the same answer on any number of threads: the
// work is cut into a fixed number of blocks, whatever the thread count, each
// block's partial result is kept apart, and they are combined in block order.
// Which thread ran a block cannot change a sum, because no sum crosses a block
// boundary except in that fixed order.

#include <cstddef>
#include <functional>

namespace mxi {

//: How many threads parallel work uses; 0 sets the machine's.
void set_parallel_threads(std::size_t threads);
std::size_t parallel_threads();

//: The fixed number of blocks work is cut into.
constexpr std::size_t kParallelBlocks = 64;

//: [begin, end) of block b of n items in `blocks` blocks.
inline std::size_t block_begin(std::size_t n, std::size_t blocks,
                               std::size_t b) {
  return n * b / blocks;
}

//: body(b) for every block b in [0, blocks), across the threads.
void for_each_block(std::size_t blocks,
                    const std::function<void(std::size_t)> &body);

//: body(i) for every i in [0, n), across the threads: for work whose items are
//: independent, each writing only its own result.
void for_each_index(std::size_t n,
                    const std::function<void(std::size_t)> &body);

} // namespace mxi
