# The documents

What each is for. A **reference** says what a part does, how to run it and how
well it works, and is kept true to the code; **notes** record invariants and
lessons for whoever changes the code; a **history** is kept as it was written,
its conclusions sometimes overturned later, with an index of those.

| document | kind | about |
|---|---|---|
| `README.md` | reference | building, the chain, the programs, where it stands |
| `docs/review.md` | guide | for reviewers: what this is, how it differs from DIALS, where things are |
| `docs/outstanding.md` | list | every open task, with its evidence |
| `docs/integration.md` | reference | `mxi_integrate` |
| `docs/symmetry.md` | reference | `mxi_symmetry` |
| `docs/scaling.md` | reference | `mxi_scale` |
| `docs/spots.md` | reference | the spot finder, `mxi_find` |
| `docs/gpu.md` | reference | where integration's, indexing's and refinement's time goes, and what a device would take; the refinement target's port, designed and not written |
| `docs/conventions.md` | reference | conventions checked against DIALS' files, and what is not |
| `CLAUDE.md` | notes | invariants and lessons across the pipeline |
| `docs/spots_notes.md` | notes | the spot finder's, from when it was a repository of its own |
| `docs/spotfinder.md` | history | how the spot finder came into this tree |
| `docs/integration_history.md` | history | integration's development, with its overturned conclusions |
| `python/README.md` | reference | `mxeq`, the comparison tools |
| `python/CLAUDE.md` | notes | `mxeq`'s |

Every path and program a document names is checked by
`python/tests/test_documents.py`, which fails on a stale one.
