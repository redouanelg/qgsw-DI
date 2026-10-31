"""The cached-grid fast path must be exactly equivalent to rebuilding it."""
from __future__ import annotations

import importlib
import time

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
