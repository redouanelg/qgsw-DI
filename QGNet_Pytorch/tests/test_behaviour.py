"""Behaviour, robustness and autograd tests that do not need qgsw-DI."""
from __future__ import annotations

import importlib

import numpy as np
import pytest
import torch

from conftest import npy, synthetic_case, t64

modgrid = importlib.import_module("QGNet_Pytorch.modgrid")
moddyn = importlib.import_module("QGNet_Pytorch.moddyn")
modelliptic = importlib.import_module("QGNet_Pytorch.modelliptic")
qgsw = importlib.import_module("QGNet_Pytorch.qgsw")

SMALL = synthetic_case(ny=24, nx=28, land=True)


def _run(case, **kw):
    cfg = dict(tint=3600.0, dtout=3600.0, dt=600.0, snu=0.0)
    cfg.update(kw)
    return qgsw.qgsw(Hi=t64(case["ssh"]), c=t64(case["c"]),
                     lon=t64(case["lon"]), lat=t64(case["lat"]), **cfg)


# ------------------------------------------------------------------ package
def test_submodules_are_importable_from_the_package():
    import QGNet_Pytorch as pkg

    assert pkg.qgsw is qgsw, "`qgsw` must be the submodule, not the function"
    assert pkg.qgsw_fn is qgsw.qgsw
    for name in ("moddyn", "modelliptic", "modgrid"):
        assert getattr(pkg, name) is importlib.import_module(f"QGNet_Pytorch.{name}")


# ------------------------------------------------------------------ helpers
def test_safe_divide_matches_the_original_branch():
    sd = modelliptic._safe_divide
    one = torch.tensor(1.0, dtype=torch.float64)
    assert float(sd(one, torch.tensor(0.0, dtype=torch.float64))) == 1.0
    assert float(sd(one, torch.tensor(2.0, dtype=torch.float64))) == 0.5


def test_safe_divide_propagates_nan_instead_of_hiding_it():
    sd = modelliptic._safe_divide
    nan = torch.tensor(float("nan"), dtype=torch.float64)
    one = torch.tensor(1.0, dtype=torch.float64)
    assert torch.isnan(sd(one, nan)), "a NaN denominator must not be read as zero"
    assert torch.isnan(sd(nan, one))


def test_nan_to_zero_preserves_infinities():
    from QGNet_Pytorch._torch_utils import nan_to_zero

    x = torch.tensor([float("nan"), float("inf"), -float("inf"), 1.0], dtype=torch.float64)
    out = npy(nan_to_zero(x))
    assert out[0] == 0.0
    assert np.isposinf(out[1]) and np.isneginf(out[2])
    assert out[3] == 1.0


def test_mixed_precision_promotes_and_warns():
    case = SMALL
    with pytest.warns(RuntimeWarning, match="mixed floating dtypes"):
        S, _ = qgsw.qgsw(
            Hi=torch.as_tensor(case["ssh"], dtype=torch.float32),
            c=t64(case["c"]), lon=t64(case["lon"]), lat=t64(case["lat"]),
            tint=600.0, dtout=600.0, dt=600.0, snu=0.0,
        )
    assert S.dtype == torch.float64, "highest precision must win"


def test_device_mismatch_is_rejected():
    from QGNet_Pytorch._torch_utils import infer_dtype_device

    if not torch.backends.mps.is_available() and not torch.cuda.is_available():
        pytest.skip("no second device available")
    dev = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cuda")
    a = torch.zeros(3, dtype=torch.float32)
    b = torch.zeros(3, dtype=torch.float32, device=dev)
    with pytest.raises(ValueError, match="same device"):
        infer_dtype_device(a, b)


# --------------------------------------------------------------- validation
@pytest.mark.parametrize("kw,msg", [
    (dict(dt=0.0), "must be positive"),
    (dict(dtout=0.0), "must be positive"),
    (dict(dtout=300.0, dt=600.0), "greater than or equal"),
    (dict(nn_cadence_steps=0), "positive integer"),
])
def test_degenerate_configurations_raise(kw, msg):
    with pytest.raises(ValueError, match=msg):
        _run(SMALL, **kw)


@pytest.mark.parametrize("missing", ["Hi", "c", "lon", "lat", "tint", "dtout", "dt"])
def test_missing_arguments_raise(missing):
    kw = dict(Hi=t64(SMALL["ssh"]), c=t64(SMALL["c"]), lon=t64(SMALL["lon"]),
              lat=t64(SMALL["lat"]), tint=600.0, dtout=600.0, dt=600.0, snu=0.0)
    kw[missing] = None
    with pytest.raises(ValueError, match=f"{missing} is required"):
        qgsw.qgsw(**kw)


def test_rappel_without_Hm_raises():
    with pytest.raises(ValueError, match="Hm is required"):
        _run(SMALL, rappel=1e-6)


def test_output_frames_are_all_written():
    """torch.empty is only safe if every frame gets assigned."""
    for tint, dtout, dt in ((7200.0, 1800.0, 600.0), (86400.0, 21600.0, 600.0),
                            (3600.0, 600.0, 600.0)):
        S, _ = _run(SMALL, tint=tint, dtout=dtout, dt=dt)
        nindex = int(tint / dtout) + 1
        assert S.shape[0] == nindex
        finite = torch.isfinite(t64(SMALL["ssh"]))
        for k in range(nindex):
            frame = S[k][finite]
            assert torch.isfinite(frame).all(), f"frame {k} not written"
            assert float(frame.abs().max()) < 10.0, f"frame {k} looks uninitialised"


# ---------------------------------------------------------------- behaviour
def test_scalar_c_is_broadcast():
    g = modgrid.grid(t64(SMALL["ssh"]), 1.5, 0.0, t64(SMALL["lon"]), t64(SMALL["lat"]))
    assert g.c.shape == t64(SMALL["ssh"]).shape
    assert torch.allclose(g.c, torch.full_like(g.c, 1.5))


def test_snu_none_differs_from_zero():
    """Documented quirk: the tendency masking only runs when snu is not None."""
    a, _ = _run(SMALL, snu=0.0)
    b, _ = _run(SMALL, snu=None)
    d = (a[-1] - b[-1]).abs()
    assert float(d[torch.isfinite(d)].max()) > 0.0


def test_land_stays_nan_and_boundary_is_frozen():
    S, _ = _run(SMALL, tint=7200.0, dtout=3600.0)
    g = modgrid.grid(t64(SMALL["ssh"]), t64(SMALL["c"]), 0.0,
                     t64(SMALL["lon"]), t64(SMALL["lat"]))
    land = g.mask == 0.0
    bnd = g.mask == 1.0
    Hi = t64(SMALL["ssh"])
    for k in range(S.shape[0]):
        assert torch.isnan(S[k][land]).all(), f"frame {k}: land is not NaN"
        assert torch.equal(S[k][bnd], Hi[bnd]), f"frame {k}: boundary drifted"


def test_backward_run_is_not_a_noop():
    f, _ = _run(SMALL, tint=7200.0, dtout=7200.0)
    b, _ = _run(SMALL, tint=-7200.0, dtout=7200.0)
    d = (f[-1] - b[-1]).abs()
    assert float(d[torch.isfinite(d)].max()) > 1e-9


def test_operators_do_not_mutate_their_inputs():
    g = modgrid.grid(t64(SMALL["ssh"]), t64(SMALL["c"]), 0.0,
                     t64(SMALL["lon"]), t64(SMALL["lat"]))
    h = t64(SMALL["ssh"])
    before = h.clone()
    modelliptic.h2pv(h, g)
    moddyn.h2uv(h, g)
    assert torch.equal(torch.nan_to_num(h), torch.nan_to_num(before))

    q, = modelliptic.h2pv(h, g)
    qb = q.clone()
    hg = torch.nan_to_num(h)
    hgb = hg.clone()
    modelliptic.pv2h(q, hg, g, 5)
    assert torch.equal(torch.nan_to_num(q), torch.nan_to_num(qb))
    assert torch.equal(hg, hgb)


def test_border_nans_do_not_expand_the_halo():
    case = synthetic_case(ny=20, nx=22, land=False)
    ssh = case["ssh"].copy()
    ssh[0, 10] = np.nan
    g = modgrid.grid(t64(ssh), t64(case["c"]), 0.0, t64(case["lon"]), t64(case["lat"]))
    # the NaN itself is land, but the 5x5 halo is not applied for border NaNs
    assert float(g.mask[0, 10]) == 0.0
    assert float(g.mask[3, 10]) == 2.0


# ----------------------------------------------------------------- autograd
def test_gradcheck_on_a_tiny_case():
    case = synthetic_case(ny=12, nx=14, land=False)

    def f(h):
        S, _ = qgsw.qgsw(Hi=h, c=t64(case["c"]), lon=t64(case["lon"]),
                         lat=t64(case["lat"]), tint=1200.0, dtout=1200.0,
                         dt=600.0, snu=0.0)
        return S[-1][3:-4, 3:-4]

    h0 = t64(case["ssh"]).requires_grad_(True)
    assert torch.autograd.gradcheck(f, (h0,), eps=1e-6, atol=1e-6, rtol=1e-4)


def test_model_is_only_piecewise_differentiable_at_the_upwind_kink():
    """Documented limitation, not a bug.

    ``qrhs`` splits the face velocity with ``where(scaled < 0, 0, scaled)``,
    which has a kink at ``scaled == 0``. A perfectly symmetric analytic SSH
    field puts face velocities exactly on that kink, so a central finite
    difference straddles it and gradcheck legitimately fails. Real fields never
    hit it. The original NumPy adjoint has the same kink.
    """
    ny, nx = 12, 14
    lon, lat = np.meshgrid(300.0 + 0.1875 * np.arange(nx), 30.0 + 0.15 * np.arange(ny))
    yy, xx = np.mgrid[0:ny, 0:nx]
    ssh = 0.4 * np.sin(2 * np.pi * xx / nx) * np.cos(2 * np.pi * yy / ny) \
        + 0.2 * np.sin(6 * np.pi * xx / nx + 1.0)
    c = np.zeros_like(ssh) + 1.5

    g = modgrid.grid(t64(ssh), t64(c), 0.0, t64(lon), t64(lat))
    with torch.no_grad():
        u, _ = moddyn.h2uv(t64(ssh), g)
    u_face = 0.5 * (u[2:-2, 2:-2] + u[2:-2, 3:-1])
    assert float(u_face.abs().min()) < 1e-15, "field is not degenerate; test is moot"

    def f(h):
        S, _ = qgsw.qgsw(Hi=h, c=t64(c), lon=t64(lon), lat=t64(lat),
                         tint=1200.0, dtout=1200.0, dt=600.0, snu=0.0)
        return S[-1][3:-4, 3:-4]

    with pytest.raises(Exception):
        torch.autograd.gradcheck(f, (t64(ssh).requires_grad_(True),),
                                 eps=1e-6, atol=1e-6, rtol=1e-4)

    # ... while an infinitesimal perturbation off the kink restores agreement.
    ssh_off = ssh + 0.05 * np.random.default_rng(0).standard_normal((ny, nx))
    assert torch.autograd.gradcheck(f, (t64(ssh_off).requires_grad_(True),),
                                    eps=1e-6, atol=1e-6, rtol=1e-4)


def test_gradients_are_finite_with_land():
    h = t64(SMALL["ssh"]).clone().requires_grad_(True)
    S, _ = qgsw.qgsw(Hi=h, c=t64(SMALL["c"]), lon=t64(SMALL["lon"]),
                     lat=t64(SMALL["lat"]), tint=1800.0, dtout=1800.0,
                     dt=600.0, snu=0.0)
    out = S[-1]
    (out[torch.isfinite(out)] ** 2).mean().backward()
    assert torch.isfinite(h.grad).all()
    assert int((h.grad != 0).sum()) > 0


# ------------------------------------------------------------ neural module
class _Zero(torch.nn.Module):
    def forward(self, x):
        return torch.zeros((1, 2, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)


class _Counter(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.calls = 0
        self.conv = torch.nn.Conv2d(5, 2, 1, bias=False).double()
        torch.nn.init.zeros_(self.conv.weight)

    def forward(self, x):
        self.calls += 1
        return self.conv(x)


def test_zero_network_reproduces_pure_qg_exactly():
    a, _ = _run(SMALL, tint=3600.0)
    b, _ = _run(SMALL, tint=3600.0, nn_module=_Zero())
    assert torch.equal(torch.nan_to_num(a), torch.nan_to_num(b))


@pytest.mark.parametrize("cadence,expected", [(1, 6), (2, 3), (3, 2), (6, 1), (10, 1)])
def test_nn_cadence(cadence, expected):
    m = _Counter()
    _run(SMALL, tint=3600.0, nn_module=m, nn_cadence_steps=cadence)
    assert m.calls == expected


def test_nn_input_has_five_channels_and_is_finite():
    seen = {}

    class Spy(torch.nn.Module):
        def forward(self, x):
            seen["shape"] = tuple(x.shape)
            seen["finite"] = bool(torch.isfinite(x).all())
            return torch.zeros((1, 2, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)

    _run(SMALL, tint=600.0, nn_module=Spy())
    ny, nx = SMALL["ssh"].shape
    assert seen["shape"] == (1, 5, ny, nx)
    assert seen["finite"], "CNN input must never contain NaN/Inf"


def test_nn_output_shape_is_validated():
    class Bad(torch.nn.Module):
        def forward(self, x):
            return torch.zeros((1, 3, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)

    with pytest.raises(ValueError, match="must return shape"):
        _run(SMALL, tint=600.0, nn_module=Bad())

    class NotATensor(torch.nn.Module):
        def forward(self, x):
            return [0.0]

    with pytest.raises(TypeError, match="torch.Tensor"):
        _run(SMALL, tint=600.0, nn_module=NotATensor())


def test_nn_nan_outside_the_valid_support_cannot_leak():
    class NanEdges(torch.nn.Module):
        def forward(self, x):
            o = torch.zeros((1, 2, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)
            o[:, :, :2, :] = float("nan")
            o[:, :, -3:, :] = float("inf")
            o[:, :, :, :2] = float("nan")
            o[:, :, :, -3:] = float("inf")
            return o

    S, _ = _run(SMALL, tint=1200.0, nn_module=NanEdges())
    ref, _ = _run(SMALL, tint=1200.0)
    assert torch.equal(torch.nan_to_num(S), torch.nan_to_num(ref))


class _ConstCorrection(torch.nn.Module):
    """Emit a constant (u_corr, v_corr) so each channel can be probed alone."""

    def __init__(self, u_amp, v_amp):
        super().__init__()
        self.u_amp = u_amp
        self.v_amp = v_amp

    def forward(self, x):
        o = torch.zeros((1, 2, x.shape[2], x.shape[3]), dtype=x.dtype, device=x.device)
        o[0, 0] = self.u_amp
        o[0, 1] = self.v_amp
        return o


@pytest.mark.parametrize("u_amp,v_amp,label", [(0.05, 0.0, "u"), (0.0, 0.05, "v")])
def test_each_correction_channel_is_actually_wired_up(u_amp, v_amp, label):
    """Guards against a channel being dropped, zeroed or swapped.

    Without this, multiplying ``u_corr`` by 0 in qgsw.py passes the whole suite:
    the zero-network and gradient tests are both blind to one channel going
    missing.
    """
    ref, _ = _run(SMALL, tint=1800.0, dtout=1800.0)
    got, _ = _run(SMALL, tint=1800.0, dtout=1800.0,
                  nn_module=_ConstCorrection(u_amp, v_amp))
    d = (ref[-1] - got[-1]).abs()
    d = d[torch.isfinite(d)]
    assert float(d.max()) > 1e-9, f"{label}_corr channel has no effect on the solution"


def test_correction_channels_are_not_swapped():
    """u_corr must reach u and v_corr must reach v, not the other way round."""
    a, _ = _run(SMALL, tint=1800.0, dtout=1800.0, nn_module=_ConstCorrection(0.05, 0.0))
    b, _ = _run(SMALL, tint=1800.0, dtout=1800.0, nn_module=_ConstCorrection(0.0, 0.05))
    d = (a[-1] - b[-1]).abs()
    d = d[torch.isfinite(d)]
    assert float(d.max()) > 1e-9, "u and v corrections are indistinguishable"


def test_correction_matches_an_equivalent_hand_applied_velocity_offset():
    """A constant correction must equal adding the same offset to h2uv output
    on the exact discrete support qrhs reads."""
    g = modgrid.grid(t64(SMALL["ssh"]), t64(SMALL["c"]), 0.0,
                     t64(SMALL["lon"]), t64(SMALL["lat"]))
    h = t64(SMALL["ssh"])
    q, = modelliptic.h2pv(h, g)
    u, v = moddyn.h2uv(h, g)

    valid_u = qgsw._neural_u_valid_mask(g)
    valid_v = qgsw._neural_v_valid_mask(g)
    u_man = u + torch.where(valid_u, torch.full_like(u, 0.05), torch.zeros_like(u))
    v_man = v + torch.where(valid_v, torch.full_like(v, -0.02), torch.zeros_like(v))
    rq_man, = moddyn.qrhs(u_man, v_man, q, g, 1.0)

    u_c, v_c = qgsw._compute_neural_velocity_correction(
        qgsw._neural_input(h, qgsw._static_neural_channels(g)),
        g, _ConstCorrection(0.05, -0.02))
    u_api, v_api = qgsw._apply_neural_velocity_correction(u, v, u_c, v_c, valid_u, valid_v)
    rq_api, = moddyn.qrhs(u_api, v_api, q, g, 1.0)

    assert torch.equal(rq_man, rq_api)


def test_nn_gradients_flow():
    net = torch.nn.Conv2d(5, 2, 3, padding=1).double()
    torch.nn.init.zeros_(net.bias)
    S, _ = _run(SMALL, tint=1800.0, dtout=1800.0, nn_module=net)
    out = S[-1]
    (out[torch.isfinite(out)] ** 2).mean().backward()
    assert torch.isfinite(net.weight.grad).all()
    assert float(net.weight.grad.abs().max()) > 0.0
