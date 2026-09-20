"""Finding the reindexing operator between two indexing solutions.

Two runs of indexing on the same images can produce solutions related by a
change of basis: the same lattice, described differently.  Compare the Miller
indices directly and almost none of them agree, which looks like catastrophe
and is not.  So the operator has to be found *before* the indices are joined,
and it has to be found from the data rather than assumed.

The method uses the fact that both tables describe the same observed spots.
Match the rows spatially first -- positions are unambiguous -- and then, over
the matched pairs, count how many satisfy ``h_a . M == h_b`` for each candidate
``M``.  The right operator is the one that agrees for nearly all of them, and
if none does then the two solutions genuinely disagree.

The operator is looked for in two ways: enumerated over the lattice point
groups, and solved for directly by least squares over the matched pairs. The
second is needed because the first only finds symmetry-equivalent reindexings,
while two different reduced cells of one lattice are related by a general
unimodular matrix that belongs to no point group.

On conventions: whether the operator acts on Miller indices as a row-vector
product, a column-vector product, or with a transpose somewhere is a question
this deliberately does not have to answer.  The candidate pool is the union of
the cubic and hexagonal lattice groups, explicitly closed under transposition
and inversion (see :func:`candidate_operators`), so the operator is in the pool
under any of those conventions and the search finds it.  Only the *reported*
matrix depends on the convention, and it is reported as the matrix that was
applied.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    import gemmi
except ImportError:  # pragma: no cover - gemmi is a declared dependency
    gemmi = None


def _close(generators: list[np.ndarray], limit: int = 200) -> list[np.ndarray]:
    """The multiplicative closure of a set of integer matrices."""
    group = {np.eye(3, dtype=np.int64).tobytes(): np.eye(3, dtype=np.int64)}
    frontier = [np.eye(3, dtype=np.int64)] + [g.astype(np.int64) for g in generators]
    for g in generators:
        group[g.astype(np.int64).tobytes()] = g.astype(np.int64)
    while frontier:
        new = []
        for a in frontier:
            for g in generators:
                p = (a @ g.astype(np.int64)).astype(np.int64)
                key = p.tobytes()
                if key not in group:
                    group[key] = p
                    new.append(p)
        frontier = new
        if len(group) > limit:
            raise RuntimeError("generators do not close to a finite group")
    return list(group.values())


#: 4-fold about c, 3-fold about the body diagonal, and inversion: the cubic
#: lattice group, 48 operators, which contains every signed permutation of the
#: axes and so every reindexing available to a lattice with orthogonal axes.
_CUBIC = [
    np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
    np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
    -np.eye(3, dtype=int),
]

#: 6-fold about c, a 2-fold, and inversion: the hexagonal lattice group, 24
#: operators, needed because a hexagonal or rhombohedral lattice's own
#: symmetry is not a subgroup of the cubic one.
_HEXAGONAL = [
    np.array([[1, -1, 0], [1, 0, 0], [0, 0, 1]]),
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
    -np.eye(3, dtype=int),
]


def candidate_operators() -> list[np.ndarray]:
    """Every integer operator a reindexing could plausibly be.

    The pool is closed under transposition and inversion *explicitly*, and not
    because the groups it is built from happen to be.  The cubic group is: its
    metric is the identity, so ``M^T = M^-1`` and the transpose of a member is
    a member.  The hexagonal group is not, for exactly the same reason in
    reverse -- its metric is not the identity, so transposes fall outside it.
    That asymmetry is what would have made the search convention-dependent for
    hexagonal and rhombohedral lattices only, which is the kind of bug that
    survives a long time because it passes on every test case someone thinks
    to write down.

    The closure is not a group and does not need to be.  This is a list of
    matrices to try.
    """
    seen: dict[bytes, np.ndarray] = {}
    for generators in (_CUBIC, _HEXAGONAL):
        for m in _close(generators):
            seen.setdefault(m.tobytes(), m)

    while True:
        extra: dict[bytes, np.ndarray] = {}
        for m in seen.values():
            for candidate in (
                m.T.astype(np.int64),
                np.rint(np.linalg.inv(m.astype(float))).astype(np.int64),
            ):
                key = candidate.tobytes()
                if key not in seen:
                    extra[key] = candidate
        if not extra:
            break
        seen.update(extra)
    return list(seen.values())


@dataclass
class Reindexing:
    operator: np.ndarray
    fraction: float
    is_identity: bool
    #: Fraction agreeing under the identity, for comparison.  If the best
    #: operator barely beats this, the search has found nothing.
    fraction_identity: float
    n_tested: int
    #: Whether the operator leaves the crystal's own metric invariant.  An
    #: operator that does not is not a lattice symmetry, and a high match
    #: fraction for one is a result that wants looking at rather than using.
    metric_compatible: bool | None = None

    def apply(self, hkl: np.ndarray) -> np.ndarray:
        return (np.asarray(hkl, dtype=np.int64) @ self.operator).astype(np.int64)


def find_operator(
    hkl_a: np.ndarray,
    hkl_b: np.ndarray,
    real_space_a: np.ndarray | None = None,
    metric_tolerance: float = 1e-3,
) -> Reindexing:
    """The operator mapping ``hkl_a`` onto ``hkl_b``, over already-matched rows.

    ``hkl_a`` and ``hkl_b`` must be row-aligned: row *i* of each is the same
    observed spot, as established by a spatial match.
    """
    hkl_a = np.atleast_2d(np.asarray(hkl_a, dtype=np.int64))
    hkl_b = np.atleast_2d(np.asarray(hkl_b, dtype=np.int64))
    if hkl_a.shape != hkl_b.shape:
        raise ValueError("hkl_a and hkl_b must be row-aligned")

    # Unindexed reflections carry (0, 0, 0) and would agree under every
    # operator, so they are excluded rather than allowed to flatten the search.
    keep = np.any(hkl_a != 0, axis=1) & np.any(hkl_b != 0, axis=1)
    a, b = hkl_a[keep], hkl_b[keep]
    n = len(a)

    identity = np.eye(3, dtype=np.int64)
    if n == 0:
        return Reindexing(identity, 0.0, True, 0.0, 0)

    best = identity
    best_score = -1.0
    identity_score = 0.0
    for m in candidate_operators():
        score = float(np.mean(np.all((a @ m) == b, axis=1)))
        if np.array_equal(m, identity):
            identity_score = score
        if score > best_score:
            best, best_score = m, score

    # The enumerated pool only contains lattice *point group* operators, which
    # covers reindexings between symmetry-equivalent settings and nothing else.
    # Two perfectly good reduced cells of the same lattice can be related by a
    # general unimodular matrix that is in no point group -- on a rhombohedral
    # insulin cell the operator between one indexing and another came out as
    # [[1,-1,0],[0,-1,1],[0,-1,0]], which the pool does not contain, and the
    # search reported two per cent agreement and no operator found.
    #
    # So solve for it instead of enumerating: h_b = h_a M is linear in M, and
    # least squares over a few thousand pairs recovers it to a part in 1e15.
    # Rounded, it must be integer and unimodular to be a change of basis at
    # all, and it still has to beat the pool on agreement before it is used.
    if n >= 3:
        solved, *_ = np.linalg.lstsq(a.astype(float), b.astype(float), rcond=None)
        rounded = np.rint(solved)
        integral = np.abs(solved - rounded).max() < 0.05
        unimodular = abs(abs(float(np.linalg.det(rounded))) - 1.0) < 1e-6
        if integral and unimodular:
            candidate = rounded.astype(np.int64)
            score = float(np.mean(np.all((a @ candidate) == b, axis=1)))
            if score > best_score:
                best, best_score = candidate, score

    metric_ok = None
    if real_space_a is not None:
        g = (
            np.asarray(real_space_a, dtype=float)
            @ np.asarray(real_space_a, dtype=float).T
        )
        transformed = best.T.astype(float) @ g @ best.astype(float)
        scale = np.abs(g).max() or 1.0
        metric_ok = bool(np.all(np.abs(transformed - g) / scale < metric_tolerance))

    return Reindexing(
        operator=best,
        fraction=best_score,
        is_identity=bool(np.array_equal(best, identity)),
        fraction_identity=identity_score,
        n_tested=n,
        metric_compatible=metric_ok,
    )


def point_group_operators(hall: str | None) -> list[np.ndarray] | None:
    """The rotation parts of a space group's operations, from its Hall symbol.

    Returned so that a found operator can be labelled as being inside the
    space group -- in which case it is a symmetry equivalence and not a
    reindexing at all -- or outside it.
    """
    if hall is None or gemmi is None:
        return None
    try:
        ops = gemmi.symops_from_hall(hall)
    except (ValueError, RuntimeError):
        return None
    out: dict[bytes, np.ndarray] = {}
    for op in ops.sym_ops:
        m = (np.array(op.rot, dtype=np.int64) // op.DEN).astype(np.int64)
        out.setdefault(m.tobytes(), m)
    return list(out.values())


def space_group_name(hall: str | None) -> str | None:
    """A readable Hermann-Mauguin name for a Hall symbol, if gemmi knows one."""
    if hall is None or gemmi is None:
        return None
    try:
        sg = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(hall))
    except (ValueError, RuntimeError):
        return None
    return sg.hm if sg is not None else None
