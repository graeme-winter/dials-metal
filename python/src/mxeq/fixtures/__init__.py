"""Making test data.

Synthetic sweeps, synthetic NXmx files, and cut-down real ones. These were
scattered across `tests/` and `spotfinder/tests/` as loose scripts; a fixture
that two test suites need belongs to neither of them.

They import h5py, which nothing else here does, so the checker still runs where
it is absent.
"""
