"""The ``nn_channels`` contract and the paper's reference network."""
from __future__ import annotations

import importlib

import numpy as np
import pytest
import torch

from conftest import synthetic_case, t64

modgrid = importlib.import_module("QGNet_Pytorch.modgrid")
qgsw = importlib.import_module("QGNet_Pytorch.qgsw")
models = importlib.import_module("QGNet_Pytorch.models")

CASE = synthetic_case(ny=24, nx=28, land=True)
CFG = dict(tint=1200.0, dtout=1200.0, dt=600.0, snu=0.0)


def _args():
    return dict(Hi=t64(CASE["ssh"]), c=t64(CASE["c"]),
                lon=t64(CASE["lon"]), lat=t64(CASE["lat"]))


class Spy(torch.nn.Module):
    """Records what it was fed and returns a zero correction."""

    def __init__(self):
        super().__init__()
        self.inputs = []

    def forward(self, x):
        self.inputs.append(x.detach().clone())
        return torch.zeros((1, 2, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)


# ------------------------------------------------------------------ defaults
def test_default_is_ssh_only():
    """arXiv:1911.08856 feeds the ConvNet SSH alone; that must be the default."""
    assert qgsw.DEFAULT_NN_CHANNELS == ("ssh",)
    m = Spy()
    qgsw.qgsw(**_args(), **CFG, nn_module=m)
    assert m.inputs[0].shape[1] == 1


def test_available_channels_are_advertised():
    assert qgsw.NEURAL_CHANNELS == ("ssh", "mask", "f0", "dx", "dy")


# ------------------------------------------------------------------ selection
@pytest.mark.parametrize("spec,width", [
    (("ssh",), 1),
    (("ssh", "mask"), 2),
    (("ssh", "mask", "f0"), 3),
    (("ssh", "mask", "f0", "dx", "dy"), 5),
    ("ssh", 1),                       # bare string is accepted
])
def test_channel_count_matches_the_request(spec, width):
    m = Spy()
    qgsw.qgsw(**_args(), **CFG, nn_module=m, nn_channels=spec)
    assert m.inputs[0].shape[1] == width


def test_channels_are_stacked_in_the_requested_order():
    grd = modgrid.grid(t64(CASE["ssh"]), t64(CASE["c"]), 0.0,
                       t64(CASE["lon"]), t64(CASE["lat"]))
    expected = {
        "ssh": torch.nan_to_num(t64(CASE["ssh"]), nan=0.0),
        "mask": grd.mask / 2.0,
    }
    for spec in (("ssh", "mask"), ("mask", "ssh")):
        m = Spy()
        qgsw.qgsw(**_args(), **CFG, nn_module=m, nn_channels=spec)
        got = m.inputs[0][0]
        for k, name in enumerate(spec):
            assert torch.equal(got[k], expected[name]), f"{spec}: channel {k} is not {name}"


def test_requested_channels_have_the_documented_content():
    grd = modgrid.grid(t64(CASE["ssh"]), t64(CASE["c"]), 0.0,
                       t64(CASE["lon"]), t64(CASE["lat"]))
    m = Spy()
    qgsw.qgsw(**_args(), **CFG, nn_module=m,
              nn_channels=("ssh", "mask", "f0", "dx", "dy"))
    x = m.inputs[0][0]
    assert torch.isfinite(x).all(), "CNN input must never contain NaN/Inf"
    # mask channel: land 0, boundary 0.5, interior 1
    assert torch.equal(x[1], grd.mask / 2.0)
    assert set(np.unique(x[1].numpy()).tolist()) <= {0.0, 0.5, 1.0}
    # standardised geometry channels are zero outside the active interior
    active = grd.mask == 2.0
    for k in (2, 3, 4):
        assert torch.equal(x[k][~active], torch.zeros_like(x[k][~active]))


def test_grid_channels_are_constant_in_time():
    """The motivation for excluding them by default: they carry no state."""
    m = Spy()
    qgsw.qgsw(**_args(), tint=3600.0, dtout=3600.0, dt=600.0, snu=0.0,
              nn_module=m, nn_channels=("ssh", "mask", "f0", "dx", "dy"))
    assert len(m.inputs) > 1
    for later in m.inputs[1:]:
        assert torch.equal(later[0, 1:], m.inputs[0][0, 1:]), "grid channels changed in time"
    assert not torch.equal(m.inputs[-1][0, 0], m.inputs[0][0, 0]), "ssh channel should evolve"


# ----------------------------------------------------------------- validation
@pytest.mark.parametrize("bad,msg", [
    (("ssh", "bogus"), "unknown nn_channels"),
    (("nope",), "unknown nn_channels"),
    (("ssh", "ssh"), "duplicate nn_channels"),
    ((), "at least one channel"),
])
def test_bad_channel_specs_are_rejected(bad, msg):
    with pytest.raises(ValueError, match=msg):
        qgsw.qgsw(**_args(), **CFG, nn_module=Spy(), nn_channels=bad)


def test_module_in_channels_mismatch_is_caught_early():
    net = torch.nn.Conv2d(5, 2, 3, padding=1).double()
    with pytest.raises(ValueError, match="in_channels is 5 but nn_channels requests 1"):
        qgsw.qgsw(**_args(), **CFG, nn_module=net, nn_channels=("ssh",))


def test_module_without_in_channels_is_not_rejected():
    qgsw.qgsw(**_args(), **CFG, nn_module=Spy(), nn_channels=("ssh", "mask"))


# --------------------------------------------------------- behavioural parity
def test_channel_choice_does_not_change_a_zero_correction():
    """Selecting channels must only change the CNN input, never the physics."""
    ref, _ = qgsw.qgsw(**_args(), **CFG)
    for spec in (("ssh",), ("ssh", "mask"), ("ssh", "mask", "f0", "dx", "dy")):
        got, _ = qgsw.qgsw(**_args(), **CFG, nn_module=Spy(), nn_channels=spec)
        assert torch.equal(torch.nan_to_num(ref, nan=-1e30),
                           torch.nan_to_num(got, nan=-1e30)), spec


# ------------------------------------------------------- paper reference model
def test_reference_convnet_has_the_published_parameter_count():
    """arXiv:1911.08856 section 3: 'QG-Net has 2545 trainable parameters'."""
    net = models.QGNetConvNet()
    assert sum(p.numel() for p in net.parameters()) == 2545
    assert net.in_channels == 1


def test_reference_convnet_shapes_and_zero_start():
    net = models.QGNetConvNet().double()
    x = torch.randn(1, 1, 20, 24, dtype=torch.float64)
    out = net(x)
    assert out.shape == (1, 2, 20, 24)
    assert torch.equal(out, torch.zeros_like(out)), "scale must start at zero"


def test_reference_convnet_rejects_the_wrong_width():
    net = models.QGNetConvNet().double()
    with pytest.raises(ValueError, match="expected 1 input channel"):
        net(torch.randn(1, 3, 10, 10, dtype=torch.float64))


def test_reference_convnet_runs_in_the_solver_and_trains():
    net = models.QGNetConvNet().double()
    first, _ = qgsw.qgsw(**_args(), **CFG, nn_module=net)
    ref, _ = qgsw.qgsw(**_args(), **CFG)
    assert torch.equal(torch.nan_to_num(first, nan=-1e30),
                       torch.nan_to_num(ref, nan=-1e30)), "zero-init must reproduce pure QG"

    with torch.no_grad():                       # leave the zero-correction manifold
        net.scale.fill_(1.0)
        net.head.weight.normal_(0.0, 0.05)
    out, _ = qgsw.qgsw(**_args(), **CFG, nn_module=net)
    with torch.no_grad():
        d = (out[-1] - ref[-1]).abs()
        assert float(d[torch.isfinite(d)].max()) > 1e-9, "trained net must change the forecast"

    loss = (out[-1][torch.isfinite(out[-1])] ** 2).mean()
    loss.backward()
    grads = [p.grad for p in net.parameters() if p.grad is not None]
    assert grads and all(torch.isfinite(g).all() for g in grads)
    assert max(float(g.abs().max()) for g in grads) > 0.0
