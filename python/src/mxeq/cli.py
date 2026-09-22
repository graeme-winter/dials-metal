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


def _sniff(path):
    """Is this a reflection table or an experiment list? None if neither.

    By the first byte, which is enough and does not require reading either
    format. A msgpack reflection table begins 0x93, the header for a
    three-element array; an experiment list is JSON and begins with whitespace
    or a brace.
    """
    try:
        with open(path, "rb") as handle:
            first = handle.read(1)
    except OSError:
        return None
    if not first:
        return None
    if first[0] == 0x93:
        return "refl"
    if first[0] in b"{ \t\r\n":
        return "expt"
    return None


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

    t = sub.add_parser(
        "trend",
        help="where two integrated files disagree, binned against what might explain it",
    )
    t.add_argument("a", help="ours")
    t.add_argument("b", help="theirs, the reference")
    t.add_argument(
        "--value",
        action="append",
        help="a column to compare; repeatable. Defaults to the intensities, "
        "their variances and the background.",
    )
    t.add_argument(
        "--value-b",
        help="the column to compare against, if different. Pass the same file "
        "twice to compare our own summed and fitted intensities.",
    )
    t.add_argument("--bins", type=int, default=10, help="bins per variable")
    t.add_argument(
        "--worst",
        type=int,
        default=0,
        help="also list the N reflections that disagree most, with their "
        "indices and positions",
    )
    t.add_argument(
        "--radius",
        type=float,
        default=0.5,
        help="how many images apart two rows may be and still be the same observation",
    )

    hr = sub.add_parser(
        "html",
        help="an HTML report of the same comparison as trend, with graphs",
    )
    hr.add_argument("a", help="ours")
    hr.add_argument("b", help="theirs, the reference")
    hr.add_argument("-o", "--output", default="comparison.html", help="where to write")
    hr.add_argument(
        "--value",
        action="append",
        help="a column to compare; repeatable. Defaults to the intensities, "
        "their variances and the background.",
    )
    hr.add_argument("--bins", type=int, default=12, help="bins per variable")
    hr.add_argument("--title", default="integration comparison")
    hr.add_argument(
        "--plotly",
        default=None,
        help="where the page should load plotly from; by default its CDN, "
        "which needs a network connection once to draw",
    )

    pl = sub.add_parser(
        "profiles", help="draw the reference profiles mxi_integrate learned"
    )
    pl.add_argument("path", help="what --save-profiles wrote")
    pl.add_argument(
        "-o", "--output", default="profiles.png", help="where to write the sections"
    )
    pl.add_argument(
        "--widths", help="also write a plot of profile width by region and axis"
    )
    pl.add_argument(
        "--block",
        type=int,
        help="draw only this block of the scan; with a divided scan there can "
        "be more profiles than fit on a page",
    )

    dg = sub.add_parser(
        "disagree",
        help="write a reflection table of only the reflections two "
        "integrations disagree about, to look at in the image viewer",
    )
    dg.add_argument("ours", help="our integrated.refl")
    dg.add_argument("theirs", help="the reference integrated.refl")
    dg.add_argument("-o", "--output", default="disagree.refl", help="where to write")
    dg.add_argument(
        "--value", default="intensity.sum.value", help="the column to compare"
    )
    dg.add_argument(
        "--value-b",
        help="the column to compare against, if different; pass the same file "
        "twice to compare our own summed and fitted intensities",
    )
    dg.add_argument(
        "--factor",
        type=float,
        default=2.0,
        help="how far apart counts as disagreeing, either way round (2)",
    )
    dg.add_argument(
        "--absolute",
        type=float,
        default=5.0,
        help="ignore pairs where both are below this, since a ratio between "
        "two numbers near zero means nothing (5)",
    )
    dg.add_argument(
        "--limit", type=int, default=0, help="keep only the N worst; 0 for all"
    )

    i = sub.add_parser("inspect", help="describe a file without assuming its layout")
    i.add_argument("path")

    args = parser.parse_args(argv)

    if args.command == "html":
        from . import htmlreport

        values = args.value or [
            "intensity.sum.value",
            "intensity.prf.value",
            "background.mean",
        ]
        print(
            htmlreport.write(
                refl.load(args.a),
                refl.load(args.b),
                args.output,
                values,
                n_bins=args.bins,
                title=args.title,
                plotly_src=args.plotly or htmlreport.PLOTLY_CDN,
            )
        )
        return 0

    if args.command == "disagree":
        from . import disagree

        table, report = disagree.select(
            refl.load(args.ours),
            refl.load(args.theirs),
            value=args.value,
            value_b=args.value_b,
            factor=args.factor,
            absolute=args.absolute,
            limit=args.limit,
        )
        print(report)
        if table.nrows == 0:
            print("nothing to write")
            return 0
        refl.write(args.output, table)
        print(f"wrote {args.output}")
        print(f"  dials.image_viewer imported.expt {args.output}")
        return 0

    if args.command == "profiles":
        from .plots import profiles as profile_plots

        reference = profile_plots.read_profiles(args.path)
        print(
            f"{len(reference.profiles)} profiles, {reference.side} a side, "
            f"{reference.divisions}x{reference.divisions} regions per panel"
        )
        thin = [r for r, n in enumerate(reference.spots) if n < 10]
        if thin:
            print(f"  {len(thin)} cells saw fewer than ten spots: {thin[:12]}")
        print(
            f"  {reference.blocks} scan blocks, "
            f"{reference.divisions}x{reference.divisions} cells a panel"
        )
        print(profile_plots.draw_profiles(reference, args.output, args.block))
        if args.widths:
            print(profile_plots.draw_profile_widths(reference, args.widths))
        return 0

    if args.command == "trend":
        from . import trends

        values = args.value or [
            "intensity.sum.value",
            "intensity.prf.value",
            "intensity.sum.variance",
            "intensity.prf.variance",
            "background.mean",
        ]
        print(
            trends.compare(
                refl.load(args.a),
                refl.load(args.b),
                values,
                value_b=args.value_b,
                n_bins=args.bins,
                radius=args.radius,
                worst=args.worst,
            )
        )
        return 0

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

    # Which kind of file each boundary wants, checked before anything tries to
    # parse one as the other. `check refined` compares experiment lists and the
    # rest compare reflection tables, and passing the wrong pair used to fail
    # inside a UTF-8 decoder with the offending byte's position -- a message
    # that says nothing about what went wrong or what to do.
    # `check refined` compares models, and models live in .expt. But comparing
    # two refined pipelines usually means asking how their residuals compare,
    # and that lives in the .refl -- so a pair of reflection tables is accepted
    # and routed to the comparison that answers it, rather than refused on a
    # technicality.
    kinds = {_sniff(args.a), _sniff(args.b)}
    if args.boundary == "refined" and kinds == {"refl"}:
        print(
            "mxeq: comparing reflections, since both files are tables; "
            "pass the .expt files to compare the models instead",
            file=sys.stderr,
        )
        args.boundary = "indexed"

    for path in (args.a, args.b):
        kind = _sniff(path)
        wanted = "expt" if args.boundary == "refined" else "refl"
        if kind is not None and kind != wanted:
            print(
                f"mxeq: {path} looks like a{'n' if kind == 'expt' else ''} "
                f"{'experiment list (.expt)' if kind == 'expt' else 'reflection table (.refl)'}, "
                f"and 'check {args.boundary}' compares "
                f"{'experiment lists' if wanted == 'expt' else 'reflection tables'}.",
                file=sys.stderr,
            )
            if wanted == "expt":
                print(
                    "       'check refined' compares the models -- cell, "
                    "orientation, detector -- so it wants the .expt files.",
                    file=sys.stderr,
                )
            else:
                print(
                    f"       'check {args.boundary}' compares reflections, so "
                    "it wants the .refl files; the .expt goes in -e.",
                    file=sys.stderr,
                )
            return 2

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
