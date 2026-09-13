"""``mxeq`` -- equivalence checks between two runs of an MX pipeline.

    mxeq check strong      a.refl b.refl
    mxeq check indexed     a.refl b.refl -e a.expt
    mxeq check refined     a.expt b.expt
    mxeq check integrated  a.refl b.refl -e a.expt [--operator ...]
    mxeq check scaled      a.refl b.refl -e a.expt [--operator ...]
    mxeq inspect           file.refl

Exit status is zero whenever the check ran.  It is not a verdict: this version
applies no thresholds, and a non-zero status means the comparison could not be
made, not that the two disagreed.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from . import expt, refl
from .checks import indexed, integrated, refined, scaled, strong

BOUNDARIES = ("strong", "indexed", "refined", "integrated", "scaled")


def _operator(text: str | None) -> np.ndarray | None:
    if not text:
        return None
    parts = [p for p in text.replace(",", " ").split() if p]
    if len(parts) != 9:
        raise SystemExit(f"--operator needs nine integers, got {len(parts)}")
    return np.array([int(p) for p in parts], dtype=np.int64).reshape(3, 3)


def _load_expt(path: str | None) -> expt.ExperimentList | None:
    return None if path is None else expt.load(path)


def _guess(table: refl.ReflectionTable) -> str:
    if "intensity.scale.value" in table or "inverse_scale_factor" in table:
        return "scaled"
    if "intensity.prf.value" in table or "intensity.sum.variance" in table:
        return "integrated"
    if "miller_index" in table:
        return "indexed"
    return "strong"


def main(argv: list[str] | None = None) -> int:
    # Reports are long by design and `mxeq check ... | head` is the normal way
    # to read one. Python's default SIGPIPE handling turns that into a
    # traceback on stderr; restoring the Unix default makes it do nothing,
    # which is what every other command line tool does.
    try:
        import signal

        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except (AttributeError, ValueError):  # not POSIX, or not the main thread
        pass

    parser = argparse.ArgumentParser(
        prog="mxeq", description="Compare two runs of an MX pipeline, stage by stage."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    c = sub.add_parser("check", help="compare two files at a pipeline boundary")
    c.add_argument("boundary", choices=(*BOUNDARIES, "auto"))
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("-e", "--expt", help="experiment list, for resolution and symmetry")
    c.add_argument("--expt-b", help="B's experiment list, if it differs from A's")
    c.add_argument(
        "--operator",
        help="nine integers: the reindexing operator from `mxeq check indexed`",
    )
    c.add_argument("--radius", type=float, default=2.0, help="spatial match radius, px")
    c.add_argument("--bins", type=int, default=10, help="number of resolution shells")
    c.add_argument("--json", action="store_true", help="emit JSON instead of text")

    i = sub.add_parser("inspect", help="describe a file without assuming its layout")
    i.add_argument("path")

    args = parser.parse_args(argv)

    if args.command == "inspect":
        with open(args.path, "rb") as f:
            raw = f.read()
        if raw[:1] in (b"{", b"\n", b" "):
            experiments = expt.load(args.path)
            print(f"experiment list, {len(experiments)} experiments")
            for n, e in enumerate(experiments):
                print(f"  [{n}] identifier {e.identifier!r}")
                if e.crystal:
                    cell = " ".join(f"{v:.4f}" for v in e.crystal.cell)
                    print(f"      cell {cell}")
                    print(f"      hall {e.crystal.hall!r}")
                    print(f"      scan varying: {e.crystal.scan_varying}")
                if e.scan:
                    print(f"      images {e.scan.image_range} osc {e.scan.oscillation}")
            return 0
        try:
            table = refl.loads(raw)
        except refl.ReflFormatError as exc:
            print(f"could not read as a reflection table: {exc}\n", file=sys.stderr)
            print("\n".join(refl.structure(raw)))
            return 2
        print(f"reflection table, {table.nrows} rows")
        print(f"identifiers: {table.identifiers}")
        for name in sorted(table.columns):
            values = table.columns[name]
            shape = "" if values.ndim == 1 else f" x{values.shape[1]}"
            print(f"  {name:<34} {table.types[name]}{shape}")
        for name, type_name in sorted(table.opaque.items()):
            print(f"  {name:<34} {type_name}  (not decoded)")
        return 0

    boundary = args.boundary
    experiments_a = _load_expt(args.expt)
    experiments_b = _load_expt(args.expt_b or args.expt)
    operator = _operator(args.operator)

    if boundary == "refined":
        report = refined.check(expt.load(args.a), expt.load(args.b), args.a, args.b)
    else:
        table_a, table_b = refl.load(args.a), refl.load(args.b)
        if boundary == "auto":
            boundary = _guess(table_a)
            print(
                f"mxeq: guessing boundary {boundary!r} from A's columns",
                file=sys.stderr,
            )
        if boundary == "strong":
            report = strong.check(table_a, table_b, args.a, args.b, radius=args.radius)
        elif boundary == "indexed":
            report = indexed.check(
                table_a, table_b, args.a, args.b, experiments_a, radius=args.radius
            )
        elif boundary == "integrated":
            report = integrated.check(
                table_a,
                table_b,
                args.a,
                args.b,
                experiments_a,
                operator=operator,
                n_bins=args.bins,
            )
        else:
            report = scaled.check(
                table_a,
                table_b,
                args.a,
                args.b,
                experiments_a or experiments_b,
                operator=operator,
                n_bins=args.bins,
            )

    print(report.json() if args.json else report.text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
