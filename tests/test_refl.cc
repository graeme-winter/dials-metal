// The reflection table reader and writer, at the level of bytes.

#include <fstream>
#include <string>

#include "../src/refl.h"
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

  Table::Opaque big;
  big.type = "Shoebox<>";
  big.rows = 1;
  try {
    // 4 GB plus one byte. If the allocation itself fails this machine cannot
    // run the test, which is said rather than reported as a pass.
    big.bytes.assign(static_cast<std::size_t>(0xFFFFFFFFull) + 1, '\0');
  } catch (const std::bad_alloc &) {
    check::is_true(true, "not enough memory to test the limit here");
    return;
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
