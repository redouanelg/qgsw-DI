"""Machine-precision parity against the original NumPy qgsw-DI implementation.

NOT bit-level: ``torch.cos``/``torch.sin``/``torch.dot`` do not agree with their
NumPy counterparts in the last place, so a handful of grid cells and the CG line
search differ by 1-2 ULP. Everything stays at machine epsilon and does not
amplify over thousands of steps. ``test_regression.py`` is the bit-level guard.

These tests are skipped when qgsw-DI cannot be found. Point QGSW_DI_PATH at a
checkout to enable them.
"""
from __future__ import annotations

import numpy as np
import pytest
import torch

from conftest import assert_close_with_nans, npy, synthetic_case, t64

# Tolerances. Everything here should land at machine epsilon; the values are
# the measured worst cases with ~4x headroom.
TOL_OP = 1e-13
TOL_RUN = 1e-13
TOL_GRID = 1e-13

CASES = {
    "ssh": None,  # filled from the ssh_case fixture
    "synthetic_land": synthetic_case(land=True),
    "synthetic_noland": synthetic_case(land=False),
}


@pytest.fixture(params=["ssh", "synthetic_land", "synthetic_noland"])
def case(request, ssh_case):
    if request.param == "ssh":
        return ssh_case
    return CASES[request.param]


@pytest.fixture
def grids(case, numpy_model, torch_model):
    n_modgrid = numpy_model[0]
    t_modgrid = torch_model[0]
    gn = n_modgrid.grid(case["ssh"].copy(), case["c"].copy(), 0.0,
                        case["lon"].copy(), case["lat"].copy())
    gt = t_modgrid.grid(t64(case["ssh"]), t64(case["c"]), 0.0,
                        t64(case["lon"]), t64(case["lat"]))
    return gn, gt


# --------------------------------------------------------------------- grid
FLOAT_ATTRS = ["mask", "dx", "dy", "f0", "c", "lon", "lat",
               "mask1d", "H", "c1d", "f01d", "dx1d", "dy1d"]
INT_ATTRS = ["indi", "indj", "vp1", "vp2", "indp",
             "vp2n", "vp2s", "vp2e", "vp2w", "vp2nn", "vp2ss", "vp2ee",
             "vp2ww", "vp2nw", "vp2ne", "vp2se", "vp2sw"]


@pytest.mark.parametrize("attr", FLOAT_ATTRS)
def test_grid_float_attributes(grids, attr):
    gn, gt = grids
    assert_close_with_nans(getattr(gn, attr), npy(getattr(gt, attr)),
                           tol=TOL_GRID, name=attr)


@pytest.mark.parametrize("attr", INT_ATTRS)
def test_grid_index_attributes_are_exact(grids, attr):
    gn, gt = grids
    assert np.array_equal(np.asarray(getattr(gn, attr)), npy(getattr(gt, attr))), attr


def test_grid_scalars(grids):
    gn, gt = grids
    for attr in ("np", "np2", "nx", "ny"):
        assert int(getattr(gn, attr)) == int(getattr(gt, attr)), attr


@pytest.mark.parametrize("attr", ["dx", "dy", "f0", "dx1d", "dy1d", "f01d"])
def test_metric_attributes_are_within_a_few_ulp(grids, attr):
    """The float metrics are *not* bit-identical, only ULP-identical.

    ``dx`` picks up ``torch.cos`` vs ``math.cos`` and ``f0`` picks up
    ``torch.sin`` vs ``numpy.sin``; both round differently in the last place.
    This test pins how large that may get, so a real regression in the metric
    formulas cannot hide behind "it's just floating point".
    """
    gn, gt = grids
    x = np.asarray(getattr(gn, attr), dtype="float64")
    y = npy(getattr(gt, attr))
    d = np.abs(x - y)
    ulp = np.where(x != 0, d / np.spacing(np.abs(x) + 1e-300), 0.0)
    assert ulp.max() <= 4.0, f"{attr} differs by {ulp.max():.1f} ULP"


def test_index_tables_and_mask_are_exactly_equal(grids):
    """The integer bookkeeping and the land mask admit no floating-point excuse."""
    gn, gt = grids
    for attr in ("mask", "indi", "indj", "vp1", "vp2", "indp",
                 "vp2n", "vp2s", "vp2e", "vp2w"):
        assert np.array_equal(np.asarray(getattr(gn, attr)), npy(getattr(gt, attr))), attr


# ----------------------------------------------------------------- operators
def test_h2pv(case, grids, numpy_model, torch_model):
    gn, gt = grids
    qn, = numpy_model[2].h2pv(case["ssh"].copy(), gn)
    qt, = torch_model[2].h2pv(t64(case["ssh"]), gt)
    assert_close_with_nans(qn, npy(qt), tol=TOL_OP, name="h2pv")


def test_h2uv(case, grids, numpy_model, torch_model):
    gn, gt = grids
    un, vn = numpy_model[1].h2uv(case["ssh"].copy(), gn)
    ut, vt = torch_model[1].h2uv(t64(case["ssh"]), gt)
    assert_close_with_nans(un, npy(ut), tol=TOL_OP, name="u")
    assert_close_with_nans(vn, npy(vt), tol=TOL_OP, name="v")


@pytest.mark.parametrize("way", [1.0, -1.0])
def test_qrhs(case, grids, numpy_model, torch_model, way):
    gn, gt = grids
    qn, = numpy_model[2].h2pv(case["ssh"].copy(), gn)
    un, vn = numpy_model[1].h2uv(case["ssh"].copy(), gn)
    rn, = numpy_model[1].qrhs(un.copy(), vn.copy(), qn.copy(), gn, way)
    rt, = torch_model[1].qrhs(t64(un), t64(vn), t64(qn), gt, way)
    assert_close_with_nans(rn, npy(rt), tol=TOL_OP, name=f"qrhs(way={way})")


@pytest.mark.parametrize("nitr", [1, 2, 5, 25])
def test_pv2h(case, grids, numpy_model, torch_model, nitr):
    gn, gt = grids
    qn, = numpy_model[2].h2pv(case["ssh"].copy(), gn)
    hg = np.nan_to_num(case["ssh"], nan=0.0) * 1.01
    hn, = numpy_model[2].pv2h(qn.copy(), hg.copy(), gn, nitr)
    ht, = torch_model[2].pv2h(t64(qn), t64(hg), gt, nitr)
    assert_close_with_nans(hn, npy(ht), tol=TOL_OP, name=f"pv2h(nitr={nitr})")


# --------------------------------------------------------------- integrations
RUN_CONFIGS = [
    pytest.param(dict(tint=5 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0), id="forward5d"),
    pytest.param(dict(tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=None), id="snu_none"),
    pytest.param(dict(tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=50.0), id="diffusive"),
    pytest.param(dict(tint=-2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0), id="backward"),
    pytest.param(dict(tint=86400.0, dtout=21600.0, dt=600.0, snu=0.0), id="subdaily_output"),
]


@pytest.mark.parametrize("cfg", RUN_CONFIGS)
def test_full_integration(case, numpy_model, torch_model, cfg):
    Sn, _ = numpy_model[3].qgsw(Hi=case["ssh"].copy(), c=case["c"].copy(),
                                lon=case["lon"].copy(), lat=case["lat"].copy(), **cfg)
    St, _ = torch_model[3].qgsw(Hi=t64(case["ssh"]), c=t64(case["c"]),
                                lon=t64(case["lon"]), lat=t64(case["lat"]), **cfg)
    assert_close_with_nans(Sn, npy(St), tol=TOL_RUN, name="SSH")


def test_nudged_integration(ssh_case, numpy_model, torch_model):
    cfg = dict(tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0,
               rappel=1.0 / (5 * 86400.0))
    Sn, _ = numpy_model[3].qgsw(Hi=ssh_case["ssh"].copy(), c=ssh_case["c"].copy(),
                                lon=ssh_case["lon"].copy(), lat=ssh_case["lat"].copy(),
                                Hm=ssh_case["sshm"].copy(), **cfg)
    St, _ = torch_model[3].qgsw(Hi=t64(ssh_case["ssh"]), c=t64(ssh_case["c"]),
                                lon=t64(ssh_case["lon"]), lat=t64(ssh_case["lat"]),
                                Hm=t64(ssh_case["sshm"]), **cfg)
    assert_close_with_nans(Sn, npy(St), tol=TOL_RUN, name="SSH nudged")


@pytest.mark.slow
def test_long_run_does_not_amplify_roundoff(ssh_case, numpy_model, torch_model):
    cfg = dict(tint=30 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0)
    Sn, _ = numpy_model[3].qgsw(Hi=ssh_case["ssh"].copy(), c=ssh_case["c"].copy(),
                                lon=ssh_case["lon"].copy(), lat=ssh_case["lat"].copy(), **cfg)
    St, _ = torch_model[3].qgsw(Hi=t64(ssh_case["ssh"]), c=t64(ssh_case["c"]),
                                lon=t64(ssh_case["lon"]), lat=t64(ssh_case["lat"]), **cfg)
    err = assert_close_with_nans(Sn, npy(St), tol=1e-12, name="SSH 30d")
    assert err < 1e-13, f"roundoff amplified to {err:.3e} over 4320 steps"


def test_observation_sampling_matches_where_the_original_works(ssh_case, numpy_model, torch_model):
    """The 2019 original crashed if an obs fell in the first window; compare the rest."""
    rng = np.random.default_rng(1)
    obs = np.stack([
        rng.uniform(ssh_case["lon"].min() + 1, ssh_case["lon"].max() - 1, 40),
        rng.uniform(ssh_case["lat"].min() + 1, ssh_case["lat"].max() - 1, 40),
        rng.uniform(600.0, 2 * 86400.0, 40),   # avoid the t=0 window
    ], axis=1)
    cfg = dict(tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0)
    _, hgn = numpy_model[3].qgsw(Hi=ssh_case["ssh"].copy(), c=ssh_case["c"].copy(),
                                 lon=ssh_case["lon"].copy(), lat=ssh_case["lat"].copy(),
                                 obsspace=obs.copy(), **cfg)
    _, hgt = torch_model[3].qgsw(Hi=t64(ssh_case["ssh"]), c=t64(ssh_case["c"]),
                                 lon=t64(ssh_case["lon"]), lat=t64(ssh_case["lat"]),
                                 obsspace=t64(obs), **cfg)
    assert_close_with_nans(hgn, npy(hgt), tol=1e-12, name="obs")
    assert np.isfinite(npy(hgt)).sum() > 0


def test_observation_in_the_first_window(ssh_case, numpy_model, torch_model):
    """An obs at t=0 used to hit `UnboundLocalError` in the original.

    Unpatched copies of qgsw-DI reference `h` before it is assigned; patched
    copies sample the initial condition, which is what the port does.
    """
    obs = np.array([[290.0, 36.0, 0.0], [291.0, 37.0, 86400.0]])
    cfg = dict(tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0)
    _, hgt = torch_model[3].qgsw(Hi=t64(ssh_case["ssh"]), c=t64(ssh_case["c"]),
                                 lon=t64(ssh_case["lon"]), lat=t64(ssh_case["lat"]),
                                 obsspace=t64(obs), **cfg)
    assert torch.isfinite(hgt).all()
    try:
        _, hgn = numpy_model[3].qgsw(Hi=ssh_case["ssh"].copy(), c=ssh_case["c"].copy(),
                                     lon=ssh_case["lon"].copy(), lat=ssh_case["lat"].copy(),
                                     obsspace=obs.copy(), **cfg)
    except UnboundLocalError:
        pytest.skip("qgsw-DI copy still has the unpatched obsspace bug")
    assert_close_with_nans(hgn, npy(hgt), tol=1e-12, name="obs at t=0")
