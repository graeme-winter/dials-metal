"""mxi_scale --anomalous keeps Friedel mates apart in scaling, and only there.

With a strong anomalous signal, I(+) and I(-) of an acentric reflection
differ; merged into one group, the difference is taken for error -- the error
model inflates every sigma and outlier rejection throws the largest differences
away (tests/test_scale.cc measures it). --anomalous scales them as separate
groups. The merging statistics split each pair themselves, so they must report
the same groups either way. Needs MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL.
"""

import os
import re
import subprocess

import pytest

BINARY = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")


def scale(tmp_path, name, *extra):
    run = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return run.stdout


def total_unique(report):
    found = re.search(r"^Total unique\s+(\d+)", report, flags=re.M)
    assert found, "no merging table"
    return int(found[1])


@pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL",
)
def test_anomalous_splits_the_groups_scaled_and_not_the_statistics(tmp_path):
    plain = scale(tmp_path, "plain")
    apart = scale(tmp_path, "apart", "--anomalous")
    assert "Friedel mates kept apart" not in plain
    split = re.search(
        r"kept apart \(--anomalous\): (\d+) groups scaled, from (\d+)", apart
    )
    assert split, apart
    groups, unique = int(split[1]), int(split[2])
    assert unique < groups < 2 * unique, "acentric mates apart, centric ones together"
    assert total_unique(plain) == total_unique(apart) == unique
