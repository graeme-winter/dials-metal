"""Reading the reference profiles mxi_integrate writes."""

import numpy as np
import pytest

from mxeq.plots import profiles


def write_profiles(path, side=5, regions=2, spots=(100, 200)):
    values = []
    for r in range(regions):
        grid = np.zeros((side, side, side))
        grid[side // 2, side // 2, side // 2] = 1.0 + r
        values.append(grid.ravel() / grid.sum())
    with open(path, "w") as f:
        f.write("# reference profiles\n")
        f.write(f"side {side}\n")
        f.write("sigma_b 0.02\n")
        f.write("sigma_m 0.07\n")
        f.write("half_width 3\n")
        f.write("divisions 1\n")
        f.write("panels 1\n")
        for r in range(regions):
            f.write(f"profile {r} spots {spots[r]}\n")
            for v in values[r]:
                f.write(f"{float(v):.17g}\n")
    return values


def test_reads_what_the_integrator_wrote(tmp_path):
    path = tmp_path / "p.txt"
    written = write_profiles(path)
    reference = profiles.read_profiles(str(path))
    assert reference.side == 5
    assert reference.spots == [100, 200]
    assert len(reference.profiles) == 2
    for r in range(2):
        assert np.allclose(reference.profiles[r], written[r])


def test_the_grid_is_reshaped_the_way_it_was_written(tmp_path):
    # The file is written with e1 fastest and e3 slowest, so a peak planted at
    # the centre must come back at the centre. Getting this backwards would
    # transpose every picture and make a broadening in e3 look like one in e1.
    path = tmp_path / "p.txt"
    write_profiles(path, side=7, regions=1, spots=(50,))
    reference = profiles.read_profiles(str(path))
    grid = reference.grid(0)
    assert grid.shape == (7, 7, 7)
    assert grid[3, 3, 3] == pytest.approx(1.0)
    assert grid.sum() == pytest.approx(1.0)


def test_a_truncated_file_is_refused_rather_than_reshaped(tmp_path):
    # A profile with the wrong number of points would otherwise reshape into
    # something plausible and wrong.
    path = tmp_path / "p.txt"
    write_profiles(path, side=5, regions=1, spots=(10,))
    text = path.read_text().splitlines()
    path.write_text("\n".join(text[:-10]) + "\n")
    with pytest.raises(ValueError, match="disagree"):
        profiles.read_profiles(str(path))
