"""Putting the rows of two tables alongside each other.

Every comparison in this package is a join followed by a diff, and the join is
the part that goes wrong.  Two boundaries need two different joins:

* Before indexing there is no key at all.  ``strong.refl`` rows are in whatever
  order the grouping produced them, and the two tables need not even have the
  same number of rows.  The join is spatial, and its result has three parts:
  matched pairs, rows of A with no partner, and rows of B with no partner.  The
  unmatched counts are not an error term to be minimised away -- they are the
  measurement.

* After indexing there is a key, but it is only unique up to a reindexing
  operator, and on a long scan not even then.  So the keyed join resolves
  duplicate keys by a tie-break column rather than taking the first row and
  hoping.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree


@dataclass
class Matching:
    """Indices of matched rows, and how many rows went unpartnered."""

    index_a: np.ndarray
    index_b: np.ndarray
    n_a: int
    n_b: int
    #: Separation of each matched pair, in whatever units went in.
    distance: np.ndarray | None = None

    @property
    def n_matched(self) -> int:
        return len(self.index_a)

    @property
    def n_only_a(self) -> int:
        return self.n_a - self.n_matched

    @property
    def n_only_b(self) -> int:
        return self.n_b - self.n_matched

    @property
    def fraction_matched(self) -> float:
        """Matched pairs over the larger of the two tables.

        Deliberately the larger and not the smaller: a table containing one
        reflection that happens to be in the other should not score 1.0.
        """
        denominator = max(self.n_a, self.n_b)
        return self.n_matched / denominator if denominator else 0.0


def match_positions(
    xyz_a: np.ndarray, xyz_b: np.ndarray, radius: float, scale: tuple = (1.0, 1.0, 1.0)
) -> Matching:
    """Mutual nearest neighbours within ``radius``.

    ``scale`` divides each axis before the search, so that a z in images can be
    compared against an x and y in pixels on a common footing.  Matching is
    mutual -- i is matched to j only if j is also i's nearest -- which is what
    keeps a dense cluster in A from claiming the same row of B several times.
    """
    xyz_a = np.atleast_2d(np.asarray(xyz_a, dtype=float))
    xyz_b = np.atleast_2d(np.asarray(xyz_b, dtype=float))
    n_a, n_b = len(xyz_a), len(xyz_b)
    if n_a == 0 or n_b == 0:
        return Matching(
            np.empty(0, dtype=int), np.empty(0, dtype=int), n_a, n_b, np.empty(0)
        )

    s = np.asarray(scale, dtype=float)
    a, b = xyz_a / s, xyz_b / s

    tree_a, tree_b = cKDTree(a), cKDTree(b)
    d_ab, j_of_i = tree_b.query(a, k=1, distance_upper_bound=radius)
    d_ba, i_of_j = tree_a.query(b, k=1, distance_upper_bound=radius)

    # query() returns len(tree) for "nothing inside the bound", which would
    # index out of range if used unguarded.
    ok = np.isfinite(d_ab) & (j_of_i < n_b)
    i = np.nonzero(ok)[0]
    j = j_of_i[i]
    mutual = i_of_j[j] == i

    return Matching(
        index_a=i[mutual],
        index_b=j[mutual],
        n_a=n_a,
        n_b=n_b,
        distance=d_ab[i[mutual]],
    )


def _key_array(columns: list[np.ndarray]) -> np.ndarray:
    """Pack several columns into one structured array usable as a join key."""
    flat = []
    for c in columns:
        c = np.asarray(c)
        if c.ndim == 1:
            flat.append(c.astype(np.int64))
        else:
            flat.extend(c[:, k].astype(np.int64) for k in range(c.shape[1]))
    if not flat:
        raise ValueError("a keyed join needs at least one column")
    return np.ascontiguousarray(np.stack(flat, axis=1))


def match_keys(
    columns_a: list[np.ndarray],
    columns_b: list[np.ndarray],
    tie_break_a: np.ndarray | None = None,
    tie_break_b: np.ndarray | None = None,
) -> tuple[Matching, int]:
    """Join on an integer key, returning the matching and the duplicate count.

    Where a key appears more than once on either side, the rows are paired up
    by closest ``tie_break`` value -- the frame number, in practice, which is
    what distinguishes the two observations of a reflection on a scan that goes
    round more than once.  With no tie-break given, duplicate keys are counted
    and dropped rather than guessed at.
    """
    key_a = _key_array(columns_a)
    key_b = _key_array(columns_b)
    if key_a.shape[1] != key_b.shape[1]:
        raise ValueError("the two key sets have different widths")
    n_a, n_b = len(key_a), len(key_b)
    if n_a == 0 or n_b == 0:
        return (
            Matching(np.empty(0, dtype=int), np.empty(0, dtype=int), n_a, n_b),
            0,
        )

    # One code per distinct key, shared across both sides, so that the join
    # becomes integer work rather than row-by-row comparison.
    _, codes = np.unique(np.vstack([key_a, key_b]), axis=0, return_inverse=True)
    codes = codes.ravel()
    code_a, code_b = codes[:n_a], codes[n_a:]

    def ranked(code: np.ndarray, tie: np.ndarray | None):
        """Sort by key, then by tie-break, and number the rows within each key."""
        tie = np.zeros(len(code)) if tie is None else np.asarray(tie, dtype=float)
        order = np.lexsort((tie, code))
        sorted_code = code[order]
        # Index of the first row of each run, broadcast back over the run, so
        # rank counts from zero inside every key.
        new_run = np.empty(len(order), dtype=bool)
        new_run[0] = True
        np.not_equal(sorted_code[1:], sorted_code[:-1], out=new_run[1:])
        start = np.maximum.accumulate(np.where(new_run, np.arange(len(order)), 0))
        rank = np.arange(len(order)) - start
        return order, rank

    order_a, rank_a = ranked(code_a, tie_break_a)
    order_b, rank_b = ranked(code_b, tie_break_b)

    # (key, rank) is unique on each side by construction, so a plain
    # intersection pairs the first observation of a key with the first, the
    # second with the second, and so on in tie-break order.
    width = max(int(rank_a.max()), int(rank_b.max())) + 1
    joint_a = code_a[order_a].astype(np.int64) * width + rank_a
    joint_b = code_b[order_b].astype(np.int64) * width + rank_b
    _, pos_a, pos_b = np.intersect1d(joint_a, joint_b, return_indices=True)

    counts = np.bincount(code_a, minlength=len(codes))
    counts_b = np.bincount(code_b, minlength=len(codes))
    n_duplicate = int(np.sum((counts > 1) | (counts_b > 1)))

    return (
        Matching(
            index_a=order_a[pos_a],
            index_b=order_b[pos_b],
            n_a=n_a,
            n_b=n_b,
        ),
        n_duplicate,
    )
