"""Plots of what the pipeline produced.

These read the text files the C++ tools write and draw them. They are here
rather than in `docs/` because a script that is run is code, and code that is
not in the package is code nobody tests, formats or finds.

They import matplotlib, which the checker does not. Nothing else in `mxeq` may
import this subpackage, so the checker keeps working where matplotlib is not
installed.
"""
