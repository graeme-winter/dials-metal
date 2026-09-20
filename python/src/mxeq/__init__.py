"""mxeq -- equivalence checks between two runs of an MX processing pipeline.

The package exists to answer one question at each boundary of the pipeline:
given the same images, do these two implementations produce equivalent output?
Not identical -- equivalence is the target, because identity stops being
achievable somewhere around integration and insisting on it would mean either
reproducing DIALS' floating point exactly or abandoning the check.

Nothing here imports cctbx or DIALS. The files are msgpack and JSON, and they
are read directly, so the checks run anywhere.
"""

from . import expt, match, refl, reindex, report, stats

__all__ = ["expt", "match", "refl", "reindex", "report", "stats"]
__version__ = "0.1.0"
