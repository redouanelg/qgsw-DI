"""Shared fixtures for the QGNet_Pytorch test suite."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

PKG_ROOT = Path(__file__).resolve().parents[1]
REPO_PARENT = PKG_ROOT.parent
if str(REPO_PARENT) not in sys.path:
    sys.path.insert(0, str(REPO_PARENT))


def _candidate_numpy_dirs():
    env = os.environ.get("QGSW_DI_PATH")
    if env:
        yield Path(env).expanduser().resolve()
    yield REPO_PARENT              # this package lives inside the qgsw-DI repo
    yield REPO_PARENT / "qgsw-DI"  # ... or next to a separate checkout
    yield REPO_PARENT / "qgsw_DI"


@pytest.fixture(scope="session")
def numpy_model():
    """Import the original NumPy package, or skip if it is not available."""
    for d in _candidate_numpy_dirs():
        if (d / "qgsw.py").is_file() and (d / "modgrid.py").is_file():
            break
    else:
        pytest.skip("qgsw-DI not found; set QGSW_DI_PATH to enable cross-checks")

    if not hasattr(np, "NAN"):  # let an unpatched copy run on numpy>=2
        np.NAN = np.nan
    sys.path.insert(0, str(d))
    try:
        for m in ("modgrid", "moddyn", "modelliptic", "qgsw"):
            sys.modules.pop(m, None)
        import moddyn, modelliptic, modgrid, qgsw  # noqa: E402

        mods = (modgrid, moddyn, modelliptic, qgsw)
        # Do not let the bare-name modules shadow anything later on.
        for m in ("modgrid", "moddyn", "modelliptic", "qgsw"):
            sys.modules.pop(m, None)
        return mods
    finally:
        sys.path.remove(str(d))


@pytest.fixture(scope="session")
def torch_model():
    import importlib

    return (
        importlib.import_module("QGNet_Pytorch.modgrid"),
        importlib.import_module("QGNet_Pytorch.moddyn"),
        importlib.import_module("QGNet_Pytorch.modelliptic"),
        importlib.import_module("QGNet_Pytorch.qgsw"),
    )


def find_ssh_nc():
    for p in (PKG_ROOT / "ssh.nc", REPO_PARENT / "ssh.nc"):
        if p.is_file():
            return p
    return None


@pytest.fixture(scope="session")
def ssh_case():
    """The bundled Gulf Stream SSH test case as float64 NumPy arrays."""
    path = find_ssh_nc()
    if path is None:
        pytest.skip("ssh.nc not available")
    netCDF4 = pytest.importorskip("netCDF4")
    with netCDF4.Dataset(str(path)) as nc:
        ssh = np.array(nc.variables["SSH"][0, :, :], dtype="float64")
        sshm = np.array(nc.variables["SSHm"][:, :], dtype="float64")
        lon1d = np.array(nc.variables["lon"][:], dtype="float64")
        lat1d = np.array(nc.variables["lat"][:], dtype="float64")
    ssh = np.where(np.isfinite(ssh), ssh, np.nan)
    sshm = np.where(np.isfinite(sshm), sshm, np.nan)
    lon, lat = np.meshgrid(lon1d, lat1d)
    return dict(ssh=ssh, sshm=sshm, lon=lon, lat=lat, c=np.zeros_like(ssh) + 2.5)


def synthetic_case(ny=40, nx=48, seed=0, land=True):
    """A small smooth synthetic field, optionally with land blobs."""
    rng = np.random.default_rng(seed)
    lon, lat = np.meshgrid(300.0 + 0.1875 * np.arange(nx), 30.0 + 0.15 * np.arange(ny))
    yy, xx = np.mgrid[0:ny, 0:nx]
    ssh = (
        0.4 * np.sin(2 * np.pi * xx / nx) * np.cos(2 * np.pi * yy / ny)
        + 0.2 * np.sin(6 * np.pi * xx / nx + 1.0)
        + 0.05 * rng.standard_normal((ny, nx))
    )
    if land:
        ssh[ny // 2 - 3 : ny // 2 + 3, nx // 3 : nx // 3 + 5] = np.nan
        ssh[-6:, :4] = np.nan
    return dict(ssh=ssh, lon=lon, lat=lat, c=np.zeros_like(ssh) + 1.5)


def t64(x):
    return torch.as_tensor(np.ascontiguousarray(x), dtype=torch.float64)


def npy(x):
    return x.detach().cpu().numpy()


def assert_close_with_nans(a, b, tol=1e-12, name=""):
    """Compare arrays requiring an identical NaN pattern."""
    a = np.asarray(a, dtype="float64")
    b = np.asarray(b, dtype="float64")
    assert a.shape == b.shape, f"{name}: shape {a.shape} != {b.shape}"
    na, nb = np.isnan(a), np.isnan(b)
    assert np.array_equal(na, nb), f"{name}: NaN patterns differ ({int((na ^ nb).sum())} cells)"
    fa, fb = a[~na], b[~nb]
    if fa.size == 0:
        return 0.0
    err = float(np.max(np.abs(fa - fb)))
    scale = max(float(np.max(np.abs(fa))), 1.0)
    assert err <= tol * scale, f"{name}: max|Δ|={err:.3e} exceeds {tol:.1e}*{scale:.3g}"
    return err
