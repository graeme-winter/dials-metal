"""How a check says what it found.

There is no pass or fail here, and that is version one's main design decision.
A check produces named quantities; a human decides what is acceptable.  Once
the distributions have been seen on real data, thresholds can be written
against the JSON output -- which is why the JSON carries the same numbers as
the text and not a rendered summary of them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from .stats import Summary


@dataclass
class Section:
    title: str
    lines: list[str] = field(default_factory=list)
    values: dict[str, Any] = field(default_factory=dict)

    def note(self, text: str) -> None:
        self.lines.append(text)

    def scalar(self, name: str, value: Any, text: str | None = None) -> None:
        self.values[name] = value
        if text is not None:
            self.lines.append(text)
        elif isinstance(value, float):
            self.lines.append(f"  {name:<28} {value: .4g}")
        else:
            self.lines.append(f"  {name:<28} {value}")

    def summary(self, name: str, summary: Summary, fmt: str = "{: .4g}") -> None:
        self.values[name] = {
            "n": summary.n,
            "mean": summary.mean,
            "median": summary.median,
            "std": summary.std,
            "robust_std": summary.robust_std,
            "min": summary.minimum,
            "max": summary.maximum,
            "percentiles": {str(k): v for k, v in summary.percentiles.items()},
        }
        self.lines.append(summary.line(name, fmt))

    def table(
        self, header: list[str], rows: list[list[str]], name: str | None = None
    ) -> None:
        widths = [
            max(len(header[i]), *(len(r[i]) for r in rows)) if rows else len(header[i])
            for i in range(len(header))
        ]
        self.lines.append(
            "  " + "  ".join(h.rjust(widths[i]) for i, h in enumerate(header))
        )
        for row in rows:
            self.lines.append(
                "  " + "  ".join(c.rjust(widths[i]) for i, c in enumerate(row))
            )
        if name is not None:
            self.values[name] = {"header": header, "rows": rows}


@dataclass
class Report:
    boundary: str
    file_a: str
    file_b: str
    sections: list[Section] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def section(self, title: str) -> Section:
        s = Section(title)
        self.sections.append(s)
        return s

    def warn(self, text: str) -> None:
        self.warnings.append(text)

    def text(self) -> str:
        out = [
            f"mxeq {self.boundary}",
            f"  A  {self.file_a}",
            f"  B  {self.file_b}",
        ]
        for s in self.sections:
            out.append("")
            out.append(s.title)
            out.extend(s.lines)
        if self.warnings:
            out.append("")
            out.append("warnings")
            out.extend(f"  ! {w}" for w in self.warnings)
        out.append("")
        out.append("No thresholds are applied. These are measurements, not a verdict.")
        return "\n".join(out)

    def json(self) -> str:
        return json.dumps(
            {
                "boundary": self.boundary,
                "file_a": self.file_a,
                "file_b": self.file_b,
                "sections": {s.title: s.values for s in self.sections},
                "warnings": self.warnings,
            },
            indent=2,
            default=float,
        )
