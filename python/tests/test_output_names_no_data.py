"""What the programs print says nothing about one data set.

Results from the data this was developed on belong in the documents and the
commit messages. In output they mislead: mxi_profile printed insulin's
sigma_b and sigma_m beside whatever data it was given. This scans what the
programs can print -- C++ string literals outside comments, and Python's
print, help and report strings -- for names of the data sets and machines
used in development.
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
NAMES = re.compile(
    r"insulin|ins10|thaumatin|ferritin|zenodo|MacBook|M4 Max|RTX 4060|Ryzen",
    re.IGNORECASE,
)


def cpp_literals():
    for path in sorted(
        list((ROOT / "apps").glob("*.cc")) + list((ROOT / "src").rglob("*.cc"))
    ):
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            code = line.split("//", 1)[0]
            for literal in re.findall(r'"((?:[^"\\]|\\.)*)"', code):
                yield path, number, literal


def python_output():
    pattern = re.compile(r"(print\(|help=|lines\.append\()(.*)")
    for path in sorted((ROOT / "python" / "src").rglob("*.py")):
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if line.lstrip().startswith("#"):
                continue
            found = pattern.search(line)
            if found:
                yield path, number, found.group(2)


def test_no_program_prints_the_name_of_a_data_set_or_a_machine():
    hits = [
        f"{path.relative_to(ROOT)}:{number}: {text.strip()[:80]}"
        for path, number, text in list(cpp_literals()) + list(python_output())
        if NAMES.search(text)
    ]
    assert not hits, "printed text names a data set or machine:\n" + "\n".join(hits)
