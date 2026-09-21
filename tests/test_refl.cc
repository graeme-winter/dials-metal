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
