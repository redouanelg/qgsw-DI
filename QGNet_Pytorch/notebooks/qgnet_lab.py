"""Experiment code for ``SSHforecastingQGSW_PyTorch_multistep.ipynb``.

Data loading, baselines, correction models, training and the physical
diagnostics live here so the notebook can stay a readable narrative and the
long training run can also be launched from the command line::

    python qgnet_lab.py train resunet_amp              # writes checkpoints/<name>.pt/.json
    python qgnet_lab.py train resunet_amp --device mps --batch 16 --checkpoint

Every model in this file plugs into ``qgsw.qgsw(..., nn_module=...)``: it maps
a batch of SSH fields ``(B, 1, ny, nx)`` to velocity corrections
``(B, 2, ny, nx)`` on the C-grid faces where ``moddyn.h2uv`` puts ``u`` (west
face) and ``v`` (south face).

To use another dataset, pass its path, variable names and box to
:func:`load_patch`, and redo the ``(c, snu)`` tuning: ``C_TUNED`` and
``SNU_TUNED`` are only valid for the NATL60 OSMOSIS box.
"""

from __future__ import annotations

import contextlib
import functools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

NOTEBOOK_DIR = Path(__file__).resolve().parent
TORCH_ROOT = NOTEBOOK_DIR.parent
PROJECT_ROOT = TORCH_ROOT.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from QGNet_Pytorch import moddyn, modgrid  # noqa: E402
from QGNet_Pytorch.qgsw import qgsw  # noqa: E402

DATA_FILE = TORCH_ROOT / "NATL60OSMO-CJM165_y2012m06d14-y2013m10d01.1d_SSHdegrad.nc"
CHECKPOINT_DIR = TORCH_ROOT / "checkpoints"

# Experiment constants, shared by the notebook and the CLI.
PATCH = 128
TRAIN_DAYS = 365
MAX_LEAD = 5
C_DEFAULT = 1.5      # m/s, first baroclinic phase speed used by the original runs
C_TUNED = 2.0        # m/s  } best (c, snu) of a grid search on training starts,
SNU_TUNED = 30.0     # m2/s } 1-3 day RMSE; the hybrids are trained on top of it
DT = 1200            # s; same 1-day RMSE as dt=600 to 4 digits (checked in the notebook)
CADENCE_S = 3600     # the correction is re-evaluated every model hour
DTYPE = torch.float32
G = 9.81


# ----------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------
def _pick(ds, name, candidates, what):
    if name is not None:
        return name
    for cand in candidates:
        if cand in ds.variables:
            return cand
    raise KeyError(f"no {what} variable among {candidates}; pass its name explicitly")


def load_patch(patch=PATCH, *, path=DATA_FILE, ssh_var="degraded_sossheig", lon_var=None,
               lat_var=None, time_var=None, fill_value=99999.0, origin=None):
    """Crop a daily SSH archive to a ``(ny, nx)`` box.

    Args:
        patch: box size, an int for a square or ``(ny, nx)``.
        path: NetCDF file holding ``ssh_var`` as ``(time, y, x)`` in metres.
        lon_var, lat_var, time_var: variable names; by default the usual
            NEMO / CF names are tried (``nav_lon``/``lon``/``longitude``, ...).
            Longitude and latitude may be 1-D (regular grid) or 2-D.
        fill_value: value marking land, turned into NaN (``None`` to skip).
            NaNs already in the file are kept: NaN means land to the model.
        origin: ``(y0, x0)`` of the box's first row/column; default centred.
    """
    import xarray as xr

    py, px = (patch, patch) if np.isscalar(patch) else patch
    with xr.open_dataset(path) as ds:
        lon_var = _pick(ds, lon_var, ("nav_lon", "lon", "longitude", "x"), "longitude")
        lat_var = _pick(ds, lat_var, ("nav_lat", "lat", "latitude", "y"), "latitude")
        ssh_all = ds[ssh_var]
        if ssh_all.ndim != 3:
            raise ValueError(f"{ssh_var} must be (time, y, x), got dims {ssh_all.dims}")
        ny, nx = ssh_all.shape[1:]
        lon = ds[lon_var].values.astype(np.float64)
        lat = ds[lat_var].values.astype(np.float64)
        if lon.ndim == 1 and lat.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        y0, x0 = origin if origin is not None else ((ny - py) // 2, (nx - px) // 2)
        if not (0 <= y0 <= ny - py and 0 <= x0 <= nx - px):
            raise ValueError(f"a {py}x{px} box at {(y0, x0)} does not fit in {(ny, nx)}")
        sl = (slice(y0, y0 + py), slice(x0, x0 + px))
        ssh = ssh_all.values[:, sl[0], sl[1]].astype(np.float64)
        try:
            tname = _pick(ds, time_var, ("time_counter", "time", "t"), "time")
            time_ = ds[tname].values
        except KeyError:
            time_ = np.arange(ssh.shape[0])
    if fill_value is not None:
        ssh[ssh == fill_value] = np.nan
    return dict(ssh=ssh, lon=lon[sl], lat=lat[sl], time=time_)


def default_device():
    """CUDA if present, then Apple MPS, else CPU. All are float32."""
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def splits(nt, train_days=TRAIN_DAYS, max_lead=MAX_LEAD):
    """Start-day indices. Every start has ``max_lead`` target days inside its split,
    and validation and test are separated by ``max_lead`` days."""
    train = np.arange(0, train_days - max_lead)
    rest = np.arange(train_days, nt - max_lead)
    n_val = len(rest) // 2
    val = rest[: n_val - max_lead]
    test = rest[n_val:]
    return dict(train=train, val=val, test=test)


class Case:
    """Tensors and a cached grid for one domain (the land mask is constant)."""

    def __init__(self, data, c=C_DEFAULT, snu=0.0, dtype=DTYPE, device="cpu"):
        self.dtype = dtype
        self.device = torch.device(device)
        t = lambda a: torch.as_tensor(a, dtype=dtype, device=self.device)
        self.ssh = t(data["ssh"])
        self.lon = t(data["lon"])
        self.lat = t(data["lat"])
        self.time = data["time"]
        self.set_physics(c, snu)

    def set_physics(self, c, snu):
        self.c = torch.full_like(self.lon, float(c))
        self.c_value = float(c)
        self.snu = snu
        with torch.no_grad():
            self.grd = modgrid.grid(self.ssh[0], self.c, snu, self.lon, self.lat)
        self.active = self.grd.mask == 2.0

    def forecast(self, start, days, nn_module=None, dt=DT, grad=False, checkpoint=False):
        """Forecast from day ``start`` with the t=0 frame dropped.

        ``start`` is an int, giving ``(days, ny, nx)``, or a sequence of start
        days integrated together as one batch, giving ``(B, days, ny, nx)``.
        """
        ctx = torch.enable_grad() if grad else torch.no_grad()
        batched = not np.isscalar(start)
        Hi = self.ssh[torch.as_tensor(np.asarray(start, dtype=np.int64))] if batched else self.ssh[int(start)]
        with ctx:
            out, _ = qgsw(
                Hi=Hi, c=self.c, lon=self.lon, lat=self.lat,
                tint=int(days) * 86400, dtout=86400, dt=dt, snu=self.snu,
                nn_module=nn_module, nn_cadence_steps=max(1, CADENCE_S // dt),
                grd=self.grd, nn_checkpoint=checkpoint,
            )
        return out[..., 1:, :, :]

    def target(self, start, days):
        if not np.isscalar(start):
            return torch.stack([self.target(s, days) for s in start])
        return self.ssh[int(start) + 1: int(start) + 1 + int(days)]

    def mse(self, pred, target):
        d = (pred - target)[..., self.active]
        return (d ** 2).mean(dim=-1)


# ----------------------------------------------------------------------------
# geostrophic diagnostics on cell centres (numpy, float64)
# ----------------------------------------------------------------------------
def _ddx(a, dx):
    return np.gradient(a, axis=-1) / dx


def _ddy(a, dy):
    return np.gradient(a, axis=-2) / dy


def face_to_centre(du, dv):
    """Average face values to the cell-centre velocity that ``qrhs`` advects with:
    ``0.5*(u[i,j] + u[i,j+1])`` and ``0.5*(v[i,j] + v[i+1,j])``."""
    uc = np.full_like(du, np.nan)
    vc = np.full_like(dv, np.nan)
    uc[..., :, :-1] = 0.5 * (du[..., :, :-1] + du[..., :, 1:])
    vc[..., :-1, :] = 0.5 * (dv[..., :-1, :] + dv[..., 1:, :])
    return uc, vc


def geostrophic_centre(h, grd):
    """Geostrophic velocity at centres, consistent with ``moddyn.h2uv``."""
    with torch.no_grad():
        u, v = moddyn.h2uv(torch.as_tensor(h, dtype=grd.f0.dtype, device=grd.f0.device), grd)
    return face_to_centre(u.cpu().double().numpy(), v.cpu().double().numpy())


def diagnostics(h, grd):
    """Geostrophic velocity, vorticity and the candidate correction terms.

    ``ageo`` is the first-order ageostrophic velocity from the advective
    acceleration of the geostrophic flow, ``u_a = (1/f) k x (u_g . grad) u_g``.
    Its along-stream component is the gradient-wind correction: it slows the
    flow round cyclones and speeds it up round anticyclones.
    """
    f = grd.f0.cpu().double().numpy()
    dx = grd.dx.cpu().double().numpy()
    dy = grd.dy.cpu().double().numpy()
    ug, vg = geostrophic_centre(h, grd)
    zeta = _ddx(vg, dx) - _ddy(ug, dy)
    adv_u = ug * _ddx(ug, dx) + vg * _ddy(ug, dy)
    adv_v = ug * _ddx(vg, dx) + vg * _ddy(vg, dy)
    lap_u = _ddx(_ddx(ug, dx), dx) + _ddy(_ddy(ug, dy), dy)
    lap_v = _ddx(_ddx(vg, dx), dx) + _ddy(_ddy(vg, dy), dy)
    return dict(
        f=f, ug=ug, vg=vg, zeta=zeta,
        ageo_u=-adv_v / f, ageo_v=adv_u / f,
        lap_u=lap_u, lap_v=lap_v,
    )


# ----------------------------------------------------------------------------
# correction models
# ----------------------------------------------------------------------------
class GradientWindCorrection(torch.nn.Module):
    """``gamma * (1/f) k x (u_g . grad) u_g`` evaluated on the C-grid faces.

    ``gamma = 1`` is the first-order ageostrophic (gradient-wind) correction
    with no free parameter; ``learn=True`` fits the single scalar through the
    solver like any other network.
    """

    in_channels = 1

    def __init__(self, grd, gamma=1.0, learn=False):
        super().__init__()
        self.grd = grd
        g = torch.tensor(float(gamma), dtype=grd.f0.dtype)
        self.gamma = torch.nn.Parameter(g) if learn else g

    def forward(self, x):
        grd = self.grd
        h = x[:, 0]
        u, v = moddyn.h2uv(h, grd)
        # centre velocities, as qrhs sees them
        uc = 0.5 * (u[..., :, :-1] + u[..., :, 1:])
        vc = 0.5 * (v[..., :-1, :] + v[..., 1:, :])
        uc = torch.nn.functional.pad(uc, (0, 1))
        vc = torch.nn.functional.pad(vc, (0, 0, 0, 1))
        dx, dy, f = grd.dx, grd.dy, grd.f0

        def ddx(a):
            d = torch.zeros_like(a)
            d[..., :, 1:-1] = (a[..., :, 2:] - a[..., :, :-2]) / (2 * dx[:, 1:-1])
            return d

        def ddy(a):
            d = torch.zeros_like(a)
            d[..., 1:-1, :] = (a[..., 2:, :] - a[..., :-2, :]) / (2 * dy[1:-1, :])
            return d

        au = -(uc * ddx(vc) + vc * ddy(vc)) / f
        av = (uc * ddx(uc) + vc * ddy(uc)) / f
        # centre -> faces (west face of u, south face of v)
        du = torch.zeros_like(au)
        dv = torch.zeros_like(av)
        du[..., :, 1:] = 0.5 * (au[..., :, 1:] + au[..., :, :-1])
        dv[..., 1:, :] = 0.5 * (av[..., 1:, :] + av[..., :-1, :])
        return self.gamma.to(h.dtype) * torch.stack([du, dv], dim=1)


class UniformDrift(torch.nn.Module):
    """A spatially uniform velocity ``(U, V)`` added to the advecting flow: 2 parameters.

    It translates the PV field at a constant velocity on top of the QG
    dynamics, i.e. it can only change how fast and in which direction
    eddies propagate.
    """

    in_channels = 1

    def __init__(self):
        super().__init__()
        self.uv = torch.nn.Parameter(torch.zeros(2))

    def forward(self, x):
        return self.uv.to(x.dtype).view(1, 2, 1, 1).expand(x.shape[0], 2, *x.shape[-2:])


class SplitCorrection(torch.nn.Module):
    """Keep only the uniform part (``part="drift"``) or only the spatially
    varying part (``part="pattern"``) of another model's correction.

    The uniform part is the mean over the interior, per component, recomputed
    at every call, so it may still change from one day to the next.
    """

    in_channels = 1

    def __init__(self, model, part, margin=8):
        super().__init__()
        self.model = model
        self.part = part
        self.margin = margin

    def forward(self, x):
        o = self.model(x)
        m = self.margin
        mean = o[..., m:-m, m:-m].mean(dim=(-2, -1), keepdim=True)
        return mean.expand_as(o) if self.part == "drift" else o - mean


def _groups(ch, max_groups=8):
    g = min(max_groups, ch)
    while ch % g:
        g -= 1
    return g


class ResBlock(torch.nn.Module):
    def __init__(self, cin, cout, norm):
        super().__init__()
        Norm = (lambda c: torch.nn.GroupNorm(_groups(c), c)) if norm else (lambda c: torch.nn.Identity())
        self.proj = torch.nn.Conv2d(cin, cout, 1) if cin != cout else torch.nn.Identity()
        self.body = torch.nn.Sequential(
            torch.nn.Conv2d(cin, cout, 3, padding=1), Norm(cout), torch.nn.GELU(),
            torch.nn.Conv2d(cout, cout, 3, padding=1), Norm(cout),
        )
        self.act = torch.nn.GELU()

    def forward(self, x):
        return self.act(self.body(x) + self.proj(x))


class ResUNet(torch.nn.Module):
    """The notebook's residual U-Net, with two switchable design choices.

    ``amplitude_aware=False`` is the original design: SSH standardised per
    image and GroupNorm in every block. Both remove the overall amplitude of
    the input, so the output cannot scale with the strength of the flow, while
    every ageostrophic correction does (the gradient-wind term is quadratic in
    SSH). ``amplitude_aware=True`` removes only the spatial mean of SSH (the
    flow does not see it), divides by a fixed scale, and drops GroupNorm.
    """

    in_channels = 1

    def __init__(self, base=16, depth=2, amplitude_aware=True, ssh_scale=0.1):
        super().__init__()
        self.amplitude_aware = bool(amplitude_aware)
        self.ssh_scale = float(ssh_scale)
        norm = not self.amplitude_aware
        self.depth = depth
        self.enc = torch.nn.ModuleList()
        c_prev = 1
        for lvl in range(depth + 1):
            self.enc.append(ResBlock(c_prev, base * 2 ** lvl, norm))
            c_prev = base * 2 ** lvl
        self.dec = torch.nn.ModuleList(
            ResBlock(base * 2 ** lvl + base * 2 ** (lvl + 1), base * 2 ** lvl, norm)
            for lvl in reversed(range(depth))
        )
        self.head = torch.nn.Conv2d(base, 2, 1)
        torch.nn.init.zeros_(self.head.weight)
        torch.nn.init.zeros_(self.head.bias)

    def normalise(self, x):
        mean = x.mean(dim=(-2, -1), keepdim=True)
        if self.amplitude_aware:
            return (x - mean) / self.ssh_scale
        std = ((x - mean) ** 2).mean(dim=(-2, -1), keepdim=True).add(1e-6).sqrt()
        return (x - mean) / std

    def forward(self, x):
        x = self.normalise(x)
        skips = []
        f = self.enc[0](x)
        skips.append(f)
        for lvl in range(1, self.depth):
            f = self.enc[lvl](torch.nn.functional.avg_pool2d(f, 2))
            skips.append(f)
        d = self.enc[self.depth](torch.nn.functional.avg_pool2d(f, 2))
        for i, blk in enumerate(self.dec):
            s = skips[self.depth - 1 - i]
            d = torch.nn.functional.interpolate(d, size=s.shape[-2:], mode="bilinear", align_corners=False)
            d = blk(torch.cat([d, s], dim=1))
        return self.head(d)


MODEL_SPECS = {
    # name: (factory(case) -> module, lr)
    "resunet_orig": (lambda case: ResUNet(amplitude_aware=False), 1e-3),
    "resunet_amp": (lambda case: ResUNet(amplitude_aware=True), 1e-3),
    "gradwind_gamma": (lambda case: GradientWindCorrection(case.grd, gamma=0.0, learn=True), 2e-2),
    "uniform_drift": (lambda case: UniformDrift(), 2e-3),
}


def build_model(name, case):
    factory, _ = MODEL_SPECS[name]
    return factory(case).to(dtype=case.dtype, device=case.device)


# ----------------------------------------------------------------------------
# training
# ----------------------------------------------------------------------------
def curriculum(epoch, epochs, max_k):
    """Rollout length for this epoch: 1 day, then 2, ... up to ``max_k``."""
    return 1 + min(max_k - 1, epoch * max_k // epochs)


def train(name, case, idx, *, epochs=12, pairs_per_epoch=120, batch=8, micro_batch=None,
          checkpoint=False, max_k=3, lr=None, lr_drop=0.5, seed=0, log=print,
          out_dir=CHECKPOINT_DIR):
    """Train a correction model end-to-end through the solver.

    Each epoch draws ``pairs_per_epoch`` training start days, rolls the hybrid
    model out for ``k`` days (curriculum 1 -> ``max_k``) and minimises the
    mean daily SSH MSE over interior water. The checkpoint with the best
    1-day validation MSE is kept.

    ``batch`` start days make one optimizer step. They are integrated
    ``micro_batch`` at a time as one batched solver call (default: all of
    them at once) and their gradients accumulated, so ``micro_batch`` trades
    memory for speed without changing the update. ``checkpoint=True``
    recomputes the network's activations in the backward pass, cutting
    memory about 4x for ~15 % more time.
    """
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = build_model(name, case)
    lr = MODEL_SPECS[name][1] if lr is None else lr
    params = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-4 if len(params) > 1 else 0.0)
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt = out_dir / f"{name}.pt"
    hist = dict(name=name, epochs=[], config=dict(
        epochs=epochs, pairs_per_epoch=pairs_per_epoch, batch=batch, micro_batch=micro_batch,
        checkpoint=checkpoint, device=str(case.device), max_k=max_k, lr=lr,
        lr_drop=lr_drop, seed=seed, dt=DT, cadence_s=CADENCE_S, c=case.c_value, snu=case.snu,
        n_params=sum(p.numel() for p in params)))
    val_starts = idx["val"]
    best = math.inf
    prev_k = None
    t0 = time.perf_counter()
    for ep in range(epochs):
        k = curriculum(ep, epochs, max_k)
        if prev_k is not None and k != prev_k:
            for gp in opt.param_groups:
                gp["lr"] *= lr_drop
        prev_k = k
        model.train()
        starts = rng.choice(idx["train"], size=min(pairs_per_epoch, len(idx["train"])), replace=False)
        losses = []
        for b in range(0, len(starts), batch):
            opt.zero_grad(set_to_none=True)
            chunk = starts[b:b + batch]
            mb = micro_batch or len(chunk)
            for m in range(0, len(chunk), mb):
                sub = chunk[m:m + mb]
                pred = case.forecast(sub, k, model, grad=True, checkpoint=checkpoint)
                per_sample = case.mse(pred, case.target(sub, k)).mean(dim=-1)   # (len(sub),)
                (per_sample.sum() / len(chunk)).backward()
                losses.extend(per_sample.detach().cpu().tolist())
            opt.step()
        model.eval()
        val = float(np.mean(np.concatenate([
            case.mse(case.forecast(val_starts[i:i + 16], 1, model), case.target(val_starts[i:i + 16], 1))
            .cpu().numpy()[:, 0] for i in range(0, len(val_starts), 16)])))
        rec = dict(epoch=ep + 1, k=k, lr=opt.param_groups[0]["lr"], train_mse=float(np.mean(losses)),
                   val_mse_1d=val, seconds=time.perf_counter() - t0)
        if name == "gradwind_gamma":
            rec["gamma"] = float(model.gamma.detach())
        if name == "uniform_drift":
            rec["uv"] = model.uv.detach().tolist()
        if val < best:
            best = val
            torch.save(model.state_dict(), ckpt)
            rec["saved"] = True
        hist["epochs"].append(rec)
        log(json.dumps(rec))
        (out_dir / f"{name}.json").write_text(json.dumps(hist, indent=1))
    model.load_state_dict(torch.load(ckpt, map_location=case.device))
    model.eval()
    return model, hist


def load_trained(name, case):
    model = build_model(name, case)
    model.load_state_dict(torch.load(CHECKPOINT_DIR / f"{name}.pt", map_location=case.device))
    model.eval()
    hist = json.loads((CHECKPOINT_DIR / f"{name}.json").read_text())
    return model, hist


# ----------------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------------
EVAL_CHUNK = 16   # start days integrated together when evaluating


def _chunks(starts, size=EVAL_CHUNK):
    starts = np.asarray(starts, dtype=np.int64)
    for i in range(0, len(starts), size):
        yield starts[i:i + size]


def rmse_by_lead(case, starts, model=None, days=MAX_LEAD, persistence=False):
    """``(len(starts), days)`` RMSE over interior water, metres."""
    rows = []
    for sub in _chunks(starts):
        tgt = case.target(sub, days)
        if persistence:
            pred = case.ssh[torch.as_tensor(sub)][:, None].expand_as(tgt)
        else:
            pred = case.forecast(sub, days, model)
        rows.append(case.mse(pred, tgt).sqrt().cpu().double().numpy())
    return np.concatenate(rows)


def forecast_all(case, starts, model=None, days=MAX_LEAD, persistence=False):
    """``(len(starts), days, ny, nx)`` float32 forecasts; ``persistence`` repeats day 0."""
    out = []
    for sub in _chunks(starts):
        if persistence:
            f = case.ssh[torch.as_tensor(sub)][:, None].expand(-1, days, -1, -1)
        else:
            f = case.forecast(sub, days, model)
        out.append(f.float().cpu().numpy())
    return np.concatenate(out)


def targets_all(case, starts, days=MAX_LEAD):
    return np.stack([case.target(s, days).float().cpu().numpy() for s in starts])


@contextlib.contextmanager
def inversion_iterations(nitr):
    """Temporarily run ``pv2h`` with ``nitr`` CG iterations instead of the default 1."""
    from QGNet_Pytorch import modelliptic

    original = modelliptic.pv2h
    modelliptic.pv2h = functools.partial(original, nitr=nitr)
    try:
        yield
    finally:
        modelliptic.pv2h = original


def drift_landscape(case, starts, us, vs, days=3):
    """Mean RMSE ``(len(us), len(vs), days)`` of QG plus a fixed uniform drift."""
    m = UniformDrift().to(dtype=case.dtype, device=case.device)
    out = np.empty((len(us), len(vs), days))
    for i, u in enumerate(us):
        for j, v in enumerate(vs):
            m.uv.data = torch.tensor([u, v], dtype=case.dtype, device=case.device)
            out[i, j] = rmse_by_lead(case, starts, m, days=days).mean(0)
    return out


def bootstrap_ci(x, n=2000, seed=0, q=(2.5, 97.5)):
    """Percentile bootstrap CI of the mean along axis 0."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    means = np.stack([x[rng.integers(0, len(x), len(x))].mean(0) for _ in range(n)])
    return np.percentile(means, q, axis=0)


def correction_fields(model, case, h):
    """The model's face correction for SSH ``h``, averaged to centres (numpy)."""
    with torch.no_grad():
        o = model(torch.as_tensor(h, dtype=case.dtype, device=case.device).reshape(1, 1, *h.shape))[0]
    o = o.detach().cpu()
    return face_to_centre(o[0].double().numpy(), o[1].double().numpy())


def radial_spectrum(field, dx_m):
    """Isotropic power spectrum of a 2-D field (Hann-windowed, mean removed)."""
    a = np.nan_to_num(field - np.nanmean(field))
    ny, nx = a.shape
    w = np.outer(np.hanning(ny), np.hanning(nx))
    F = np.abs(np.fft.fftshift(np.fft.fft2(a * w))) ** 2
    ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dx_m))
    kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx_m))
    K = np.hypot(*np.meshgrid(kx, ky))
    kbins = np.linspace(0, K.max() / math.sqrt(2), min(ny, nx) // 2)
    which = np.digitize(K.ravel(), kbins)
    spec = np.array([F.ravel()[which == i].mean() if np.any(which == i) else np.nan
                     for i in range(1, len(kbins))])
    return 0.5 * (kbins[1:] + kbins[:-1]), spec


def pattern_shift(a, b, margin=12):
    """Displacement (cells, ``(dj, di)``) that best maps field ``a`` onto ``b``.

    Phase correlation of the windowed interior anomalies, with a parabolic
    sub-cell refinement of the peak. Positive ``dj`` is eastward, ``di``
    northward. Used to measure how fast the SSH pattern as a whole moves.
    """
    sl = (slice(margin, -margin), slice(margin, -margin))
    a = np.nan_to_num(a[sl] - np.nanmean(a[sl]))
    b = np.nan_to_num(b[sl] - np.nanmean(b[sl]))
    w = np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    Fa, Fb = np.fft.fft2(a * w), np.fft.fft2(b * w)
    R = Fb * np.conj(Fa)
    r = np.real(np.fft.ifft2(R / (np.abs(R) + 1e-12)))
    i0, j0 = np.unravel_index(np.argmax(r), r.shape)

    def refine(m, p, c0):
        n = len(m)
        lo, hi = m[(p - 1) % n], m[(p + 1) % n]
        den = lo - 2 * c0 + hi
        return p + (0.5 * (lo - hi) / den if den != 0 else 0.0)

    di = refine(r[:, j0], i0, r[i0, j0])
    dj = refine(r[i0, :], j0, r[i0, j0])
    ny, nx = r.shape
    di = di - ny if di > ny / 2 else di
    dj = dj - nx if dj > nx / 2 else dj
    return dj, di


def interior_window(shape, margin=8):
    """Boolean mask excluding ``margin`` cells next to the open boundary."""
    m = np.zeros(shape, dtype=bool)
    m[margin:-margin, margin:-margin] = True
    return m


def distance_to_border(shape):
    ny, nx = shape
    j, i = np.meshgrid(np.arange(nx), np.arange(ny))
    return np.minimum.reduce([i, j, ny - 1 - i, nx - 1 - j])


def correction_samples(model, case, starts, margin=8):
    """Pool the learned correction and the candidate physical terms over ``starts``.

    Returns flat arrays on interior cells (``margin`` away from the border):
    the correction ``du, dv``, the geostrophic flow ``ug, vg``, ``zeta/f`` and
    the candidate terms (first-order ageostrophic ``ageo``, Laplacian ``lap``).
    """
    keep = interior_window((case.grd.ny, case.grd.nx), margin)
    cols = {k: [] for k in ("du", "dv", "ug", "vg", "ro", "ageo_u", "ageo_v", "lap_u", "lap_v")}
    for s in starts:
        h = case.ssh[int(s)].cpu().double().numpy()
        d = diagnostics(h, case.grd)
        du, dv = correction_fields(model, case, h)
        vals = dict(du=du, dv=dv, ug=d["ug"], vg=d["vg"], ro=d["zeta"] / d["f"],
                    ageo_u=d["ageo_u"], ageo_v=d["ageo_v"], lap_u=d["lap_u"], lap_v=d["lap_v"])
        for k, v in vals.items():
            cols[k].append(v[keep])
    return {k: np.concatenate(v) for k, v in cols.items()}


def regress_correction(samp):
    """Decompose the learned correction: uniform drift first, then candidate terms.

    ``drift_share`` is the fraction of the correction's mean square carried by
    its spatially uniform part ``(<du>, <dv>)``. The remainder is then fitted
    on each candidate vector field with one scalar coefficient for both
    components (respecting rotation invariance); ``r2`` is the share of the
    remainder's variance each term explains, alone and all four together.
    """
    du, dv = samp["du"], samp["dv"]
    U, V = float(du.mean()), float(dv.mean())
    total = float(np.sum(du ** 2 + dv ** 2))
    y = np.concatenate([du - U, dv - V])
    sst = float(np.sum(y ** 2))
    out = dict(drift=(U, V), drift_share=1 - sst / total, terms={})
    terms = {
        "u_g (rescaling)": (samp["ug"], samp["vg"]),
        "k x u_g (cross-stream)": (-samp["vg"], samp["ug"]),
        "lap u_g (viscosity)": (samp["lap_u"], samp["lap_v"]),
        "(1/f) k x (u_g.grad)u_g (gradient wind)": (samp["ageo_u"], samp["ageo_v"]),
    }
    X = {}
    for k, (xu, xv) in terms.items():
        X[k] = np.concatenate([xu - xu.mean(), xv - xv.mean()])
    for k, x in X.items():
        a = float(x @ y / (x @ x))
        out["terms"][k] = dict(coef=a, r2=1 - float(np.sum((y - a * x) ** 2)) / sst)
    A = np.stack(list(X.values()), axis=1)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    out["terms"]["all four"] = dict(coef=None, r2=1 - float(np.sum((y - A @ coef) ** 2)) / sst)
    return out


def along_stream_fraction(cu, cv, ug, vg):
    """Component of a correction along the geostrophic flow, as a fraction of its speed."""
    return (cu * ug + cv * vg) / (ug ** 2 + vg ** 2)


def binned(x, y, edges, weights=None):
    """Weighted mean of ``y`` in bins of ``x``."""
    which = np.digitize(x, edges)
    w = np.ones_like(y) if weights is None else weights
    out = np.full(len(edges) - 1, np.nan)
    for i in range(1, len(edges)):
        sel = which == i
        if sel.sum() > 50:
            out[i - 1] = np.sum(w[sel] * y[sel]) / np.sum(w[sel])
    return out


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Train a QG correction model on the NATL60 OSMOSIS box.")
    ap.add_argument("command", choices=["train"])
    ap.add_argument("model", choices=sorted(MODEL_SPECS))
    ap.add_argument("threads", nargs="?", type=int, default=4, help="CPU threads (default 4)")
    ap.add_argument("--device", default="cpu", help="cpu, mps, cuda or auto (default cpu)")
    ap.add_argument("--batch", type=int, default=8, help="start days per optimizer step")
    ap.add_argument("--micro-batch", type=int, default=None, help="start days per solver call")
    ap.add_argument("--checkpoint", action="store_true", help="recompute network activations")
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--pairs-per-epoch", type=int, default=120)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    device = default_device() if a.device == "auto" else a.device
    data = load_patch()
    case = Case(data, c=C_TUNED, snu=SNU_TUNED, device=device)
    idx = splits(len(data["ssh"]))
    train(a.model, case, idx, epochs=a.epochs, pairs_per_epoch=a.pairs_per_epoch, batch=a.batch,
          micro_batch=a.micro_batch, checkpoint=a.checkpoint, log=lambda r: print(r, flush=True))
