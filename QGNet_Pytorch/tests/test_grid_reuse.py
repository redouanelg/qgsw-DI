"""The cached-grid fast path must be exactly equivalent to rebuilding it."""
from __future__ import annotations

import importlib
import time

import numpy as np
import pytest
import torch

from conftest import synthetic_case, t64

modgrid = importlib.import_module("QGNet_Pytorch.modgrid")
qgsw = importlib.import_module("QGNet_Pytorch.qgsw")

CASE = synthetic_case(ny=48, nx=52, land=True)
CFG = dict(tint=7200.0, dtout=3600.0, dt=600.0, snu=0.0)


def _args():
    return dict(Hi=t64(CASE["ssh"]), c=t64(CASE["c"]),
                lon=t64(CASE["lon"]), lat=t64(CASE["lat"]))


def test_supplied_grid_gives_identical_output():
    a, _ = qgsw.qgsw(**_args(), **CFG)
    grd = modgrid.grid(t64(CASE["ssh"]), t64(CASE["c"]), CFG["snu"],
                       t64(CASE["lon"]), t64(CASE["lat"]))
    b, _ = qgsw.qgsw(**_args(), **CFG, grd=grd)
    assert torch.equal(torch.nan_to_num(a, nan=-1e30), torch.nan_to_num(b, nan=-1e30))
    assert torch.equal(torch.isnan(a), torch.isnan(b))


def test_supplied_grid_is_reusable_across_calls():
    grd = modgrid.grid(t64(CASE["ssh"]), t64(CASE["c"]), CFG["snu"],
                       t64(CASE["lon"]), t64(CASE["lat"]))
    first, _ = qgsw.qgsw(**_args(), **CFG, grd=grd)
    for _ in range(3):
        again, _ = qgsw.qgsw(**_args(), **CFG, grd=grd)
        assert torch.equal(torch.nan_to_num(first, nan=-1e30),
                           torch.nan_to_num(again, nan=-1e30))


def test_shape_mismatch_is_rejected():
    other = synthetic_case(ny=20, nx=22, land=False)
    grd = modgrid.grid(t64(other["ssh"]), t64(other["c"]), 0.0,
                       t64(other["lon"]), t64(other["lat"]))
    with pytest.raises(ValueError, match="supplied grd has shape"):
        qgsw.qgsw(**_args(), **CFG, grd=grd)


# --------------------------------------------------------------------------
# Every cached quantity on the grid is a pure function of (Hi NaN pattern, c,
# snu, lon, lat). Reusing a grid built from different inputs must fail loudly:
# it silently changes the physics otherwise.
# --------------------------------------------------------------------------
def _grid(ssh=None, c=None, snu=0.0, lon=None, lat=None):
    return modgrid.grid(
        t64(CASE["ssh"] if ssh is None else ssh),
        t64(CASE["c"] if c is None else c),
        snu,
        t64(CASE["lon"] if lon is None else lon),
        t64(CASE["lat"] if lat is None else lat),
    )


@pytest.mark.parametrize("field,build,expect", [
    ("snu",  dict(snu=None),                              "snu"),
    ("snu",  dict(snu=50.0),                              "snu"),
    ("c",    dict(c=np.zeros_like(CASE["ssh"]) + 99.0),   "c"),
    ("lon",  dict(lon=CASE["lon"] + 5.0),                 "lon"),
    ("lat",  dict(lat=CASE["lat"] + 5.0),                 "lat"),
])
def test_mismatched_grid_is_rejected(field, build, expect):
    grd = _grid(**build)
    with pytest.raises(ValueError, match=expect):
        qgsw.qgsw(**_args(), **CFG, grd=grd)


def test_mismatched_land_mask_is_rejected():
    other = CASE["ssh"].copy()
    other[5:9, 5:9] = np.nan          # extra land -> different mask and index tables
    grd = _grid(ssh=other)
    with pytest.raises(ValueError, match="land/NaN pattern"):
        qgsw.qgsw(**_args(), **CFG, grd=grd)


def test_matching_grid_with_different_ssh_values_is_accepted():
    """Only the NaN pattern of Hi matters, not its values: that is the whole
    point of reusing one grid across many SSH snapshots of the same region."""
    other = CASE["ssh"] + 0.37
    grd = _grid(ssh=other)
    a, _ = qgsw.qgsw(**_args(), **CFG, grd=grd)
    b, _ = qgsw.qgsw(**_args(), **CFG)
    assert torch.equal(torch.nan_to_num(a, nan=-1e30), torch.nan_to_num(b, nan=-1e30))


def test_scalar_and_broadcast_c_are_equivalent_for_reuse():
    grd = _grid(c=np.zeros_like(CASE["ssh"]) + 1.5)
    args = dict(_args())
    args["c"] = torch.as_tensor(1.5, dtype=torch.float64)   # scalar form
    qgsw.qgsw(**args, **CFG, grd=grd)   # must not raise


def test_grid_build_is_not_quadratic_in_land_points():
    """The NaN halo must be a vectorised dilation, not a per-NaN Python loop."""
    import numpy as np

    big = synthetic_case(ny=256, nx=256, land=False)
    ssh = big["ssh"].copy()
    ssh[np.random.default_rng(0).random(ssh.shape) < 0.25] = np.nan
    args = (t64(ssh), t64(big["c"]), 0.0, t64(big["lon"]), t64(big["lat"]))
    modgrid.grid(*args)  # warm up
    t0 = time.perf_counter()
    modgrid.grid(*args)
    elapsed = time.perf_counter() - t0
    assert elapsed < 1.0, f"grid build took {elapsed:.2f}s for 16k NaNs"
