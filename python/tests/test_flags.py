"""The Python flag table and the C++ one are the same table.

Both are DIALS' Flags enum, dials/array_family/reflection_table.h. The Python
copy had three wrong entries for as long as nothing read it -- integrated_sum
and integrated_prf at bits 11 and 12, which are DIALS' overlapped_bg and
overlapped_fg -- so this holds every name the two share to the same bit.
"""

import pathlib
import re

from mxeq.checks.common import FLAGS

HEADER = pathlib.Path(__file__).resolve().parents[2] / "src" / "refl.hh"


def _cpp_flags():
    text = HEADER.read_text()
    out = {}
    for name, bit in re.findall(r"constexpr std::int64_t k(\w+) = 1 << (\d+);", text):
        # kForegroundIncludesBadPixels -> foreground_includes_bad_pixels
        snake = re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()
        out[snake] = 1 << int(bit)
    return out


def test_every_shared_flag_has_the_same_bit():
    cpp = _cpp_flags()
    shared = set(cpp) & set(FLAGS)
    # Enough overlap for the check to mean something.
    assert len(shared) >= 8, shared
    for name in sorted(shared):
        assert FLAGS[name] == cpp[name], (name, FLAGS[name], cpp[name])


def test_the_flags_dials_scale_depends_on_are_dials_bits():
    # dials.scale's combined intensity needs integrated_sum AND integrated_prf,
    # and a gap-crossing reflection is withheld from it by leaving the first
    # off. These are DIALS' values, from its enum.
    assert FLAGS["integrated_sum"] == 256
    assert FLAGS["integrated_prf"] == 512
    assert FLAGS["foreground_includes_bad_pixels"] == 1 << 14
    assert FLAGS["failed_during_summation"] == 1 << 19
