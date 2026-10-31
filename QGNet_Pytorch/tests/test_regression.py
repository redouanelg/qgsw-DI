"""Byte-level regression guard.

``reference_outputs.npz`` was produced by a build verified bit-identical to the
NumPy original's structure (see tests/test_parity.py for the cross-check against
qgsw-DI itself). Any refactor that changes a single bit of the float64 output
fails here, which is what makes performance work on this solver safe.

Regenerate deliberately with::

    python tests/test_regression.py --regenerate
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import torch

from conftest import find_ssh_nc, npy, synthetic_case, t64

REFERENCE = Path(__file__).with_name("reference_outputs.npz")


def _cases(ssh_case):
    ssh, sshm = ssh_case["ssh"], ssh_case["sshm"]
    lon, lat, c = ssh_case["lon"], ssh_case["lat"], ssh_case["c"]
    syn = synthetic_case(land=True)
    return {
        "fwd5d": dict(Hi=ssh, c=c, lon=lon, lat=lat,
                      tint=5 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0),
        "snuNone": dict(Hi=ssh, c=c, lon=lon, lat=lat,
                        tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=None),
        "snu50": dict(Hi=ssh, c=c, lon=lon, lat=lat,
                      tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=50.0),
        "bwd": dict(Hi=ssh, c=c, lon=lon, lat=lat,
                    tint=-2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0),
        "nudge": dict(Hi=ssh, c=c, lon=lon, lat=lat, Hm=sshm,
                      rappel=1.0 / (5 * 86400.0),
                      tint=2 * 86400.0, dtout=86400.0, dt=600.0, snu=0.0),
        "synth": dict(Hi=syn["ssh"], c=syn["c"], lon=syn["lon"], lat=syn["lat"],
                      tint=86400.0, dtout=21600.0, dt=600.0, snu=0.0),
    }


def _run(kw):
    import importlib

    qgsw = importlib.import_module("QGNet_Pytorch.qgsw")
    kw = {k: (t64(v) if isinstance(v, np.ndarray) else v) for k, v in kw.items()}
    S, _ = qgsw.qgsw(**kw)
    return npy(S)


@pytest.mark.parametrize("name", ["fwd5d", "snuNone", "snu50", "bwd", "nudge", "synth"])
def test_output_is_bit_identical_to_reference(ssh_case, name):
    if not REFERENCE.is_file():
        pytest.skip("reference_outputs.npz missing")
    ref = np.load(REFERENCE)[name]
    got = _run(_cases(ssh_case)[name])
    assert got.shape == ref.shape
    assert np.array_equal(np.isnan(got), np.isnan(ref)), f"{name}: NaN pattern changed"
    a, b = got[~np.isnan(got)], ref[~np.isnan(ref)]
    if not np.array_equal(a, b):
        d = np.abs(a - b)
        raise AssertionError(
            f"{name}: output changed. max|Δ|={d.max():.3e}, "
            f"{int((a != b).sum())}/{a.size} values differ. "
            "If this is intentional, regenerate reference_outputs.npz."
        )


if __name__ == "__main__":  # pragma: no cover
    import argparse
    import sys

    ap = argparse.ArgumentParser()
    ap.add_argument("--regenerate", action="store_true")
    args = ap.parse_args()
    if not args.regenerate:
        sys.exit("pass --regenerate to overwrite the reference")

    import netCDF4

    with netCDF4.Dataset(str(find_ssh_nc())) as nc:
        ssh = np.where(np.isfinite(nc.variables["SSH"][0]), nc.variables["SSH"][0], np.nan).astype("f8")
        sshm = np.where(np.isfinite(nc.variables["SSHm"][:]), nc.variables["SSHm"][:], np.nan).astype("f8")
        lon, lat = np.meshgrid(np.asarray(nc.variables["lon"][:], dtype="f8"),
                               np.asarray(nc.variables["lat"][:], dtype="f8"))
    case = dict(ssh=ssh, sshm=sshm, lon=lon, lat=lat, c=np.zeros_like(ssh) + 2.5)
    np.savez_compressed(REFERENCE, **{k: _run(v) for k, v in _cases(case).items()})
    print(f"wrote {REFERENCE}")
