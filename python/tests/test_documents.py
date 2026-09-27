"""The documents name only things that exist.

A path, a program or a subcommand copied into prose goes stale the first time
the thing it names moves, and nothing checks it. This is the check. It was
written after a review found `CLAUDE.md` pointing at six paths that no longer
existed, a README calling integration "not started" long after it worked, and a
generator named `tests/make_real_data.py` that had been moved -- unchanged, git
said -- into the Python package when the Python was consolidated.

A path in backticks is resolved against the repository root and against the
naming document's own directory, since `python/README.md` names files relative
to `python/`. Paths under `dials/` are DIALS' own tree. A path mentioned only as
history is written without backticks, so it reads as prose and is not checked.
"""

import pathlib
import re
import shutil

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCUMENTS = [
    "README.md",
    "CLAUDE.md",
    "python/README.md",
    "python/CLAUDE.md",
    "docs/integration.md",
    "docs/integration_history.md",
    "docs/gpu.md",
    "docs/spotfinder.md",
    "docs/scaling_plan.md",
]
PATH = re.compile(
    r"`([A-Za-z0-9_./-]+/[A-Za-z0-9_.-]+\.(?:md|cc|h|hh|py|cmake|txt|json))`"
)
PROGRAM = re.compile(r"`(mxi_[a-z]+)`")


def _programs():
    """Every mxi_ program the build defines, read from the build files."""
    names = set()
    for f in [ROOT / "CMakeLists.txt", *(ROOT / "cmake").glob("*.cmake")]:
        names |= set(re.findall(r"add_executable\(\s*(mxi_[a-z]+)", f.read_text()))
    return names


@pytest.mark.parametrize("document", DOCUMENTS)
def test_every_path_named_exists(document):
    doc = ROOT / document
    text = doc.read_text()
    missing = []
    for m in PATH.finditer(text):
        path = m.group(1)
        if path.startswith("dials/"):
            continue
        if not ((ROOT / path).exists() or (doc.parent / path).exists()):
            missing.append(path)
    assert not missing, f"{document} names paths that do not exist: {missing}"


@pytest.mark.parametrize("document", DOCUMENTS)
def test_every_program_named_is_built(document):
    built = _programs()
    assert built, "found no mxi_ programs in the build files"
    named = set(PROGRAM.findall((ROOT / document).read_text()))
    unknown = named - built
    assert not unknown, f"{document} names programs the build does not make: {unknown}"


def test_every_mxeq_subcommand_named_exists():
    from mxeq.cli import build_parser

    parser = build_parser()
    commands = set()
    for action in parser._subparsers._group_actions:
        commands |= set(action.choices)
    assert commands, "mxeq has no subcommands?"
    for document in DOCUMENTS:
        text = (ROOT / document).read_text()
        for m in re.finditer(r"`mxeq ([a-z]+)", text):
            assert m.group(1) in commands, (document, m.group(1), sorted(commands))
