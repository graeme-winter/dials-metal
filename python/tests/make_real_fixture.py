"""Cut a small fixture out of a real ``.refl``, preserving its structure.

The reason this exists rather than a call to :func:`mxeq.refl.dumps`: the first
version of the reader was wrong about the file layout in two ways, and the
round-trip test against this package's own writer passed anyway, because a
writer built from the same wrong assumption is consistent with it.  A fixture
whose *structure* comes from DIALS is the only kind that can catch that.

So this slices the payload blobs of a real file in place and re-packs the same
document with the row count reduced.  Nothing in ``mxeq`` is involved in the
writing.  The shoebox column is dropped: it is ninety per cent of the file and
no check looks inside one.

Usage, from a machine that has real output:

    python tests/make_real_fixture.py dials/strong.refl metal/strong.refl \\
        tests/data --rows 32

The B row subset is chosen by matching against A, so the pair keeps the
correspondence the full files have and the fixture's expected answer is the
full files' answer.
"""

from __future__ import annotations

import argparse
import pathlib

import msgpack
import numpy as np

TAG = "dials::af::reflection_table"
DROP = ("shoebox",)

#: Bytes per row for each type, needed to slice a blob without decoding it.
WIDTH = {
    "int": 4,
    "std::size_t": 8,
    "double": 8,
    "bool": 1,
    "vec2<double>": 16,
    "vec3<double>": 24,
    "mat3<double>": 72,
    "int6": 24,
    "cctbx::miller::index<>": 12,
}


def load_document(path: str) -> list:
    """Decode with ``raw=False``, which is what makes the re-pack faithful.

    DIALS writes text as msgpack ``str`` and payloads as ``bin``. Decoding with
    ``raw=False`` keeps that distinction -- text comes back as ``str``, blobs as
    ``bytes`` -- so packing again with ``use_bin_type=True`` reproduces the
    original encoding byte for byte, down to the tag being a ``fixstr`` rather
    than a ``str8``. Decoding with ``raw=True`` collapses both to ``bytes`` and
    the re-packed file then differs from DIALS' in its string framing, which
    defeats the point of cutting the fixture from a real file at all.
    """
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def centroids(document: list) -> np.ndarray:
    data = document[-1]["data"]
    blob = data["xyzobs.px.value"][1][1]
    return np.frombuffer(blob, dtype="<f8").reshape(-1, 3)


def slice_document(document: list, rows: np.ndarray) -> bytes:
    """Re-pack ``document`` keeping only ``rows``, by slicing the raw blobs."""
    body = dict(document[-1])
    data = {}
    for name, (type_name, payload) in body["data"].items():
        if name in DROP:
            continue
        width = WIDTH.get(type_name)
        if width is None:
            raise SystemExit(f"no width known for {type_name!r} in column {name!r}")
        _, blob = payload
        picked = b"".join(bytes(blob[r * width : (r + 1) * width]) for r in rows)
        data[name] = [type_name, [len(rows), picked]]
    body["nrows"] = len(rows)
    body["data"] = data
    return msgpack.packb([TAG, document[1], body], use_bin_type=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a", help="reference .refl, e.g. DIALS' own")
    parser.add_argument("b", help="the .refl to compare against it")
    parser.add_argument("out", help="output directory")
    parser.add_argument("--rows", type=int, default=32)
    args = parser.parse_args()

    doc_a, doc_b = load_document(args.a), load_document(args.b)
    xyz_a, xyz_b = centroids(doc_a), centroids(doc_b)

    from scipy.spatial import cKDTree

    # Spread the sample over the scan rather than taking the first n rows,
    # which would all sit on the first image or two.
    order = np.argsort(xyz_a[:, 2])
    rows_a = order[np.linspace(0, len(order) - 1, args.rows).astype(int)]

    distance, rows_b = cKDTree(xyz_b).query(xyz_a[rows_a], k=1)
    if distance.max() > 2.0:
        raise SystemExit(
            f"a sampled row of A has no partner in B within 2 px "
            f"(worst {distance.max():.3f}); the two files are not of the same data"
        )

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "dials_strong.refl").write_bytes(slice_document(doc_a, rows_a))
    (out / "other_strong.refl").write_bytes(slice_document(doc_b, rows_b))
    print(f"wrote {args.rows} rows to {out}, worst separation {distance.max():.3e} px")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
