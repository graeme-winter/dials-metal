"""One module per pipeline boundary.

Each exposes ``check(a, b, ...) -> Report``. They are separate modules and not
one parameterised function because the join, the metrics and the failure modes
are different at every boundary, and sharing a body between them would mean
the differences lived in flags.
"""

from . import indexed, integrated, refined, scaled, strong

__all__ = ["strong", "indexed", "refined", "integrated", "scaled"]
