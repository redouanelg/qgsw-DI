"""Batched integration: a leading batch dimension on Hi, and network checkpointing."""
from __future__ import annotations

import importlib

import numpy as np
import pytest
import torch

from conftest import synthetic_case, t64

modgrid = importlib.import_module("QGNet_Pytorch.modgrid")
moddyn = importlib.import_module("QGNet_Pytorch.moddyn")
modelliptic = importlib.import_module("QGNet_Pytorch.modelliptic")
qgsw = importlib.import_module("QGNet_Pytorch.qgsw")

CASE = synthetic_case(ny=24, nx=28, land=True)
B = 3


def _batch():
    """Three different SSH fields sharing CASE's land mask."""
    rng = np.random.default_rng(1)
    base = CASE["ssh"]
    fields = [base, 0.8 * base, base + 0.03 * rng.standard_normal(base.shape)]
    return t64(np.stack([np.where(np.isnan(base), np.nan, f) for f in fields]))


def _grid():
    return modgrid.grid(t64(CASE["ssh"]), t64(CASE["c"]), 0.0, t64(CASE["lon"]), t64(CASE["lat"]))


def _same(a, b, tol=1e-14):
    a, b = a.detach(), b.detach()
    assert torch.equal(torch.isnan(a), torch.isnan(b))
    d = (torch.nan_to_num(a) - torch.nan_to_num(b)).abs().max()
    scale = max(float(torch.nan_to_num(b).abs().max()), 1.0)
    assert float(d) <= tol * scale, f"max|diff| = {float(d):.3e}"


def _run(Hi, **kw):
    cfg = dict(c=t64(CASE["c"]), lon=t64(CASE["lon"]), lat=t64(CASE["lat"]),
               tint=3600.0, dtout=1800.0, dt=600.0, snu=0.0)
    cfg.update(kw)
    return qgsw.qgsw(Hi=Hi, **cfg)


# ---------------------------------------------------------------- operators
def test_operators_act_per_sample():
    g = _grid()
    H = _batch()
    q, = modelliptic.h2pv(H, g)
    u, v = moddyn.h2uv(H, g)
    for i in range(B):
        qi, = modelliptic.h2pv(H[i], g)
        ui, vi = moddyn.h2uv(H[i], g)
        assert torch.equal(q[i], qi) and torch.equal(u[i], ui) and torch.equal(v[i], vi)
    for way in (1.0, -1.0):
        rq, = moddyn.qrhs(u, v, q, g, way)
        for i in range(B):
            assert torch.equal(rq[i], moddyn.qrhs(u[i], v[i], q[i], g, way)[0])


@pytest.mark.parametrize("nitr", [1, 5])
def test_inversion_has_one_line_search_per_sample(nitr):
    g = _grid()
    H = _batch()
    q, = modelliptic.h2pv(H, g)
    guess = H + 0.01
    h, = modelliptic.pv2h(q, guess, g, nitr=nitr)
    for i in range(B):
        _same(h[i], modelliptic.pv2h(q[i], guess[i], g, nitr=nitr)[0])


# --------------------------------------------------------------- propagator
@pytest.mark.parametrize("kw", [
    dict(),
    dict(tint=-3600.0),
    dict(snu=None),
    dict(rappel=1e-5, Hm="2d"),
    dict(rappel=1e-5, Hm="3d"),
])
def test_batched_run_matches_each_sample_run_alone(kw):
    H = _batch()
    kw = dict(kw)
    hm = kw.pop("Hm", None)
    if hm == "2d":
        kw["Hm"] = t64(CASE["ssh"] * 0.9)
    elif hm == "3d":
        kw["Hm"] = H * 0.9
    S, _ = _run(H, **kw)
    assert S.shape == (B, 3, 24, 28)
    for i in range(B):
        kwi = dict(kw)
        if hm == "3d":
            kwi["Hm"] = kw["Hm"][i]
        Si, _ = _run(H[i], **kwi)
        _same(S[i], Si)


class _Conv(torch.nn.Module):
    def __init__(self):
        super().__init__()
        torch.manual_seed(0)
        self.conv = torch.nn.Conv2d(1, 2, 3, padding=1).double()
        torch.nn.init.normal_(self.conv.weight, std=0.02)
        self.seen = []

    def forward(self, x):
        self.seen.append(tuple(x.shape))
        return self.conv(x)


def test_network_sees_the_whole_batch_and_matches_per_sample():
    H = _batch()
    net = _Conv()
    S, _ = _run(H, nn_module=net, nn_cadence_steps=2)
    assert set(net.seen) == {(B, 1, 24, 28)}
    for i in range(B):
        Si, _ = _run(H[i], nn_module=net, nn_cadence_steps=2)
        _same(S[i], Si)


def test_batched_gradients_equal_the_sum_of_per_sample_gradients():
    H = _batch()
    net = _Conv()

    def loss(S):
        last = S[..., -1, :, :]
        return (torch.nan_to_num(last) ** 2).sum()

    net.zero_grad()
    loss(_run(H, nn_module=net)[0]).backward()
    g_batch = net.conv.weight.grad.clone()
    net.zero_grad()
    for i in range(B):
        loss(_run(H[i], nn_module=net)[0]).backward()
    assert torch.allclose(g_batch, net.conv.weight.grad, rtol=1e-10, atol=0.0)


def test_checkpointing_changes_neither_output_nor_gradients():
    H = _batch()
    out = {}
    for ck in (False, True):
        net = _Conv()
        S, _ = _run(H, nn_module=net, nn_checkpoint=ck)
        (torch.nan_to_num(S[:, -1]) ** 2).sum().backward()
        out[ck] = (S.detach(), net.conv.weight.grad.clone())
    assert torch.equal(torch.nan_to_num(out[False][0]), torch.nan_to_num(out[True][0]))
    assert torch.allclose(out[False][1], out[True][1], rtol=1e-12, atol=0.0)


# ------------------------------------------------------------------- guards
def test_batch_with_different_land_masks_raises():
    H = _batch().clone()
    H[1, 5, 5] = float("nan")
    with pytest.raises(ValueError, match="same land"):
        _run(H)


def test_batch_with_obsspace_raises():
    obs = t64(np.array([[CASE["lon"][5, 5], CASE["lat"][5, 5], 0.0]]))
    with pytest.raises(NotImplementedError, match="batched"):
        _run(_batch(), obsspace=obs)


def test_hm_batch_must_match_hi_batch():
    with pytest.raises(ValueError, match="Hm has shape"):
        _run(_batch(), rappel=1e-5, Hm=_batch()[:2])


def test_four_dimensional_hi_is_rejected():
    with pytest.raises(ValueError, match=r"\(B, ny, nx\)"):
        _run(_batch()[None])
