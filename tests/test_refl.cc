// The reflection table reader and writer, at the level of bytes.

#include <fstream>
#include <string>

#include "../src/refl.h"
#include "../src/shoebox.h"
#include "check.h"

using namespace mxi;

TEST(a_table_survives_a_read_and_a_write_byte_for_byte) {
  // The reader was changed from a character-at-a-time stream iterator to a
  // single sized read, and the writer given a reserve. Neither should have
  // touched a byte of what is produced, and the only way to be sure of that is
  // to compare the bytes -- a table that reads back with the right values
  // would pass whatever either change did to the framing.
  Table t;
  t.nrows = 7;
  Column &m = t.int_column("miller_index", "cctbx::miller::index<>", 3);
  for (std::size_t i = 0; i < t.nrows * 3; ++i) m.ints[i] = static_cast<int>(i) - 10;
  Column &x = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  for (std::size_t i = 0; i < t.nrows * 3; ++i) x.reals[i] = 0.5 * static_cast<double>(i);
  Column &f = t.int_column("flags", "std::size_t", 1);
  for (std::size_t i = 0; i < t.nrows; ++i) f.ints[i] = 32;
  Table::Opaque blob;
  blob.type = "Shoebox<>";
  blob.rows = t.nrows;
  blob.bytes = std::string(1234, '\x5a');
  t.set_opaque("shoebox", blob);
  t.identifiers[0] = "0123abcd-0000-1111-2222-333344445555";

  const std::string first = "test_roundtrip_a.refl";
  const std::string second = "test_roundtrip_b.refl";
  write_reflections(first, t);
  write_reflections(second, read_reflections(first));

  const auto slurp = [](const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    in.seekg(0, std::ios::end);
    const std::streamoff n = in.tellg();
    in.seekg(0, std::ios::beg);
    std::string out(static_cast<std::size_t>(n), '\0');
    in.read(out.data(), n);
    return out;
  };
  const std::string a = slurp(first);
  const std::string b = slurp(second);
  check::is_true(!a.empty(), "something was written");
  check::equal(static_cast<long long>(b.size()), static_cast<long long>(a.size()),
               "same size after a round trip");
  check::is_true(a == b, "byte for byte identical after a round trip");
}

TEST(a_column_too_large_for_msgpack_is_refused_rather_than_truncated) {
  // msgpack's largest binary type is bin32, so four gigabytes is the most a
  // column can hold. The length used to wrap silently while the bytes were
  // written anyway: a 4.87 GB table declared its shoebox column as 529638180
  // bytes, which is the true size less exactly 2^32, and everything past that
  // point was unreachable. The file looked fine and was not, and nothing found
  // out until it was walked byte by byte.
  //
  // There is no bin64 to reach for, so this cannot be written differently --
  // the column has to be smaller, and the caller has to be told.
  //
  // The test builds the smallest thing that triggers it rather than a real
  // four gigabyte table, which would not fit anywhere this runs.
  Table table;
  table.nrows = 1;
  Column &c = table.int_column("id", "int", 1);
  c.ints[0] = 0;

  // Under AddressSanitizer a four gigabyte allocation aborts rather than
  // throwing, and the rule itself is tested without one below.
#if defined(__SANITIZE_ADDRESS__)
  check::skip("a four gigabyte allocation under AddressSanitizer");
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
  check::skip("a four gigabyte allocation under AddressSanitizer");
#endif
#endif
  Table::Opaque big;
  big.type = "Shoebox<>";
  big.rows = 1;
  try {
    // 4 GB plus one byte. A machine without the memory cannot run this, and
    // it is SKIPPED: it used to assert true and report "ok" for a check that
    // was never made.
    big.bytes.assign(static_cast<std::size_t>(0xFFFFFFFFull) + 1, '\0');
  } catch (const std::bad_alloc &) {
    check::skip("not enough memory for a four gigabyte column");
  }
  table.set_opaque("shoebox", std::move(big));

  bool refused = false;
  std::string said;
  try {
    write_reflections("/dev/null", table);
  } catch (const ReflError &error) {
    refused = true;
    said = error.what();
  }
  check::is_true(refused, "writing it is refused");
  check::is_true(said.find("4 GB") != std::string::npos,
                 "and the reason names the limit");
  check::is_true(said.find("slice") != std::string::npos,
                 "and says what to do instead");
}

TEST(a_saved_shoebox_follows_dials_mask_convention) {
  // Internally a bad pixel keeps its region bit and loses only Valid, so the
  // profile fit can recover a reflection crossing a module gap. DIALS sets
  // Foreground and Background only on voxels already Valid, so a region bit
  // without Valid never occurs there -- and a table carrying one was rejected
  // as an invalid structure, however far under the size limit it was.
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = 0; box.bbox[1] = 2;
  box.bbox[2] = 0; box.bbox[3] = 2;
  box.bbox[4] = 0; box.bbox[5] = 1;
  box.data.assign(4, 1.0f);
  box.background.assign(4, 0.5f);
  box.mask = {
      static_cast<std::uint8_t>(shoebox_mask::kValid | shoebox_mask::kForeground),
      static_cast<std::uint8_t>(shoebox_mask::kValid | shoebox_mask::kBackground),
      shoebox_mask::kForeground,   // a foreground voxel in a module gap
      shoebox_mask::kBackground,   // and a background one
  };
  to_dials_convention(&box);
  check::equal(static_cast<long long>(box.mask[0]),
               static_cast<long long>(shoebox_mask::kValid | shoebox_mask::kForeground),
               "a measured foreground voxel is left alone");
  check::equal(static_cast<long long>(box.mask[1]),
               static_cast<long long>(shoebox_mask::kValid | shoebox_mask::kBackground),
               "and so is a measured background one");
  check::equal(static_cast<long long>(box.mask[2]), 0,
               "an unmeasured foreground voxel is zero, not Foreground alone");
  check::equal(static_cast<long long>(box.mask[3]), 0,
               "and an unmeasured background voxel likewise");
  for (std::uint8_t m : box.mask) {
    const bool region = (m & (shoebox_mask::kForeground | shoebox_mask::kBackground)) != 0;
    check::is_true(!region || (m & shoebox_mask::kValid) != 0,
                   "no region bit without Valid");
  }
}

TEST(the_msgpack_size_rule_holds_at_its_exact_boundary) {
  // The rule on its own, at the boundary, with nothing allocated -- so it runs
  // everywhere, under the sanitizers included, where the end-to-end test above
  // cannot. msgpack's bin32 length is four bytes: 2^32 - 1 is the largest.
  bool largest_refused = false;
  try {
    check_blob_size(0xFFFFFFFFull);
  } catch (const ReflError &) {
    largest_refused = true;
  }
  check::is_true(!largest_refused, "the largest size bin32 can describe is written");

  bool over_refused = false;
  std::string said;
  try {
    check_blob_size(0xFFFFFFFFull + 1);
  } catch (const ReflError &error) {
    over_refused = true;
    said = error.what();
  }
  check::is_true(over_refused, "one byte more is refused");
  check::is_true(said.find("4 GB") != std::string::npos, "and it says why");
  check::is_true(said.find("slice") != std::string::npos, "and what to do");
}
