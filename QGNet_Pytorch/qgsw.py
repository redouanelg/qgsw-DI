"""Forward QG shallow-water propagator implemented in PyTorch."""

from __future__ import annotations

import torch
import torch.utils.checkpoint
from scipy.interpolate import griddata

try:
    from . import moddyn, modelliptic, modgrid
    from ._torch_utils import as_tensor, infer_dtype_device, maybe_tensor_scalar, scalar_float, sign_scalar
except ImportError:  # pragma: no cover
    import moddyn
    import modelliptic
    import modgrid
    from _torch_utils import as_tensor, infer_dtype_device, maybe_tensor_scalar, scalar_float, sign_scalar


def _required(name, value):
    if value is None:
        raise ValueError(f"{name} is required")


def _observation_indices(obsspace, way, time_abs, dt):
    return torch.where((way * obsspace[:, 2] >= time_abs - dt / 2.0) & (way * obsspace[:, 2] < time_abs + dt / 2.0))[0]


def _fill_observations(out, iobs, h, lon, lat, obsspace):
    if iobs.numel() > 0:
        obs_values = griddata(
            (lon.detach().cpu().numpy().ravel(), lat.detach().cpu().numpy().ravel()),
            h.detach().cpu().numpy().ravel(),
            (
                obsspace[iobs, 0].detach().cpu().numpy().squeeze(),
                obsspace[iobs, 1].detach().cpu().numpy().squeeze(),
            ),
        )
        out[iobs] = torch.as_tensor(obs_values, dtype=out.dtype, device=out.device)


#: Channels a correction network may be fed, see ``nn_channels``.
#:
#: ``"ssh"`` is the state and is rebuilt every step. The rest are pure
#: functions of the grid: they are constant in time, so they are built once per
#: run -- and on a single fixed domain they are constant across samples too,
#: which is why they are NOT in the default (see the README).
DYNAMIC_NEURAL_CHANNELS = ("ssh",)
STATIC_NEURAL_CHANNELS = ("mask", "f0", "dx", "dy")
NEURAL_CHANNELS = DYNAMIC_NEURAL_CHANNELS + STATIC_NEURAL_CHANNELS

#: Default network input: SSH alone, as in arXiv:1911.08856.
DEFAULT_NN_CHANNELS = ("ssh",)


def _normalized_grid_channel(values, active, eps=1e-12):
    """Standardise ``values`` over the active cells, zero elsewhere.

    Written without any host-device synchronisation so it stays cheap on GPU.
    """
    values = torch.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    valid = (active & torch.isfinite(values)).to(values.dtype)
    count = valid.sum().clamp_min(1.0)
    mean = (values * valid).sum() / count
    var = (((values - mean) * valid) ** 2).sum() / count
    std = var.sqrt().clamp_min(eps)
    normalized = (values - mean) / std
    return torch.where(active, normalized, torch.zeros_like(normalized))


def normalize_channel_spec(nn_channels):
    """Validate and canonicalise an ``nn_channels`` argument to a tuple."""
    if nn_channels is None:
        nn_channels = DEFAULT_NN_CHANNELS
    if isinstance(nn_channels, str):
        nn_channels = (nn_channels,)
    spec = tuple(nn_channels)
    if not spec:
        raise ValueError("nn_channels must request at least one channel")
    unknown = [c for c in spec if c not in NEURAL_CHANNELS]
    if unknown:
        raise ValueError(
            f"unknown nn_channels {unknown}; available channels are {list(NEURAL_CHANNELS)}"
        )
    duplicated = sorted({c for c in spec if spec.count(c) > 1})
    if duplicated:
        raise ValueError(f"duplicate nn_channels: {duplicated}")
    return spec


def _build_static_channels(grd, spec):
    """Materialise the grid-only channels named in ``spec``, once per run."""
    active = grd.mask == 2.0
    sources = {"f0": grd.f0, "dx": grd.dx, "dy": grd.dy}
    built = {}
    for name in spec:
        if name == "mask":
            built[name] = grd.mask / 2.0
        elif name in sources:
            built[name] = _normalized_grid_channel(sources[name], active)
    return built


def _neural_input(h, spec, static):
    """Stack the requested channels in the order the caller listed them.

    ``h`` is ``(ny, nx)`` or ``(B, ny, nx)``; the result is ``(B, C, ny, nx)``
    with ``B = 1`` for an unbatched run.
    """
    hb = h if h.ndim == 3 else h.unsqueeze(0)
    channels = [
        torch.nan_to_num(hb, nan=0.0) if name == "ssh" else static[name].expand_as(hb)
        for name in spec
    ]
    return torch.stack(channels, dim=1)


def _check_nn_in_channels(nn_module, expected):
    """Fail early and clearly if the module obviously wants a different width.

    Many modules (``nn.Conv2d``, and anything that copies the convention)
    expose ``in_channels``. Without this the mismatch surfaces as an opaque
    convolution shape error from inside the time loop.
    """
    declared = getattr(nn_module, "in_channels", None)
    if isinstance(declared, int) and declared != expected:
        raise ValueError(
            f"nn_module.in_channels is {declared} but nn_channels requests {expected} "
            f"channel(s). Adjust the module, or pass matching nn_channels."
        )


def _compute_neural_velocity_correction(nn_input, grd, nn_module, batched=False, checkpoint=False):
    """Evaluate ``nn_module`` and split its two channels into (u_corr, v_corr).

    Returns ``(B, ny, nx)`` fields when ``batched``, ``(ny, nx)`` otherwise.
    ``checkpoint`` recomputes the module's activations in the backward pass
    instead of storing them.
    """
    if checkpoint:
        correction = torch.utils.checkpoint.checkpoint(nn_module, nn_input, use_reentrant=False)
    else:
        correction = nn_module(nn_input)
    if not isinstance(correction, torch.Tensor):
        raise TypeError("nn_module must return a torch.Tensor")
    expected_shape = (nn_input.shape[0], 2, grd.ny, grd.nx)
    if tuple(correction.shape) != expected_shape:
        raise ValueError(
            f"nn_module must return shape {expected_shape}, got {tuple(correction.shape)}"
        )
    if batched:
        return correction[:, 0], correction[:, 1]
    return correction[0, 0], correction[0, 1]


def _neural_u_valid_mask(grd):
    """The ``u`` faces that ``qrhs`` reads for active (``mask == 2``) cells.

    ``qrhs`` advects cell ``(i, j)`` with ``0.5*(u[i,j] + u[i,j+1])``, so a face
    is read when the cell on either side of it is active. Masking on the face's
    own cell alone would drop the east face of the easternmost active cells and
    halve their correction.
    """
    valid_mask = torch.zeros((grd.ny, grd.nx), dtype=torch.bool, device=grd.mask.device)
    if grd.ny > 4 and grd.nx > 4:
        active = grd.mask == 2.0
        valid_mask[:, 1:] = active[:, 1:] | active[:, :-1]
    return valid_mask


def _neural_v_valid_mask(grd):
    """The ``v`` faces that ``qrhs`` reads for active cells: ``v[i,j]``, ``v[i+1,j]``."""
    valid_mask = torch.zeros((grd.ny, grd.nx), dtype=torch.bool, device=grd.mask.device)
    if grd.ny > 4 and grd.nx > 4:
        active = grd.mask == 2.0
        valid_mask[1:, :] = active[1:, :] | active[:-1, :]
    return valid_mask


def _apply_neural_velocity_correction(u, v, u_corr, v_corr, valid_u, valid_v):
    """Add the network correction on the valid support only.

    ``torch.where`` is used rather than a multiplicative mask so that a NaN or
    Inf produced by the network outside the support cannot leak into ``u``/``v``
    (``nan * 0.0`` is still ``nan``).
    """
    zero_u = torch.zeros((), dtype=u.dtype, device=u.device)
    zero_v = torch.zeros((), dtype=v.dtype, device=v.device)
    u = u + torch.where(valid_u, u_corr.to(dtype=u.dtype, device=u.device), zero_u)
    v = v + torch.where(valid_v, v_corr.to(dtype=v.dtype, device=v.device), zero_v)
    return u, v


def qgsw(
    Hi=None,
    c=None,
    lon=None,
    lat=None,
    tint=None,
    dtout=None,
    dt=None,
    obsspace=None,
    Hm=None,
    rappel=None,
    snu=None,
    nn_module=None,
    nn_cadence_steps=1,
    nn_channels=DEFAULT_NN_CHANNELS,
    grd=None,
    nn_checkpoint=False,
):
    """Run the forward QG shallow-water model.

    Explicit forward Euler on potential vorticity followed by one elliptic
    inversion per step::

        u, v   = h2uv(h^n)
        q^n+1  = q^n + dt * qrhs(u, v, q^n)          [- dt*rappel*(q^n - Qm)]
        h^n+1  = pv2h(q^n+1, hguess=2*h^n - h^n-1)

    Args:
        Hi: initial SSH, metres, ``(ny, nx)`` or a batch ``(B, ny, nx)`` of
            independent initial states integrated together. NaN marks land;
            every sample of a batch must share the same land mask.
        c: first baroclinic phase speed, m/s. 2-D or scalar.
        lon, lat: 2-D degree grids.
        tint: integration time in seconds. A negative value integrates
            backwards (the upwind sense is flipped via ``way = sign(tint)``).
        dtout: output period in seconds. Should be an integer multiple of
            ``dt``; otherwise output frames are not evenly spaced in time.
        dt: time step in seconds.
        obsspace: optional ``(N, 3)`` array of ``[lon, lat, time]``.
        Hm: mean SSH used as the nudging target. Required when ``rappel`` is
            set. ``(ny, nx)``, or ``(B, ny, nx)`` for one target per sample.
            It may be NaN on land only; NaN over water raises.
        rappel: Newtonian nudging rate in s^-1 applied to PV.
        snu: PV diffusivity in m^2/s. Pass ``0.`` rather than ``None`` -- see
            the warning in :func:`moddyn.qrhs`.
        nn_module: optional correction network, see below.
        nn_cadence_steps: how often the network is evaluated.
        nn_checkpoint: recompute the network's activations in the backward
            pass instead of storing them (``torch.utils.checkpoint``). Memory
            then grows with the number of solver steps, not of network
            activations. Every call is evaluated twice, so BatchNorm running
            statistics are updated twice per call in ``train()`` mode.
        nn_channels: sequence of channel names fed to ``nn_module``, stacked in
            the order given. Defaults to ``("ssh",)``. See
            :data:`NEURAL_CHANNELS` for the full set.
        grd: optional prebuilt :class:`modgrid.grid`. Reuse it across calls
            that share the same geometry, land mask, ``c`` and ``snu`` -- grid
            construction is otherwise repeated on every call, which dominates
            the cost of short rollouts during training. The grid is validated
            against ``Hi``/``c``/``snu``/``lon``/``lat`` and a mismatch raises,
            because every cached quantity on it depends on those inputs.
            Do not reuse a grid built from ``c``/``lat``/``lon``/``snu``
            tensors that require grad: its cached coefficients carry that
            graph, so the second ``backward()`` fails with "Trying to backward
            through the graph a second time". Rebuild it every iteration when
            learning those parameters.

    Returns:
        ``(SSH, hg)`` where ``SSH`` has shape ``(int(|tint|/dtout)+1, ny, nx)``,
        or ``(B, int(|tint|/dtout)+1, ny, nx)`` for a batched ``Hi``, and
        ``hg`` holds the SSH sampled at ``obsspace`` (or ``None``; not
        supported for a batch).

    Batched runs give each sample its own PV-inversion line search, so a
    sample agrees with the same sample run alone to round-off (~1e-16); the
    network sees the whole batch at once, ``(B, C, ny, nx)``, and must return
    ``(B, 2, ny, nx)``.

    Inputs can be NumPy arrays or torch tensors. Outputs are torch tensors on
    the inferred input device. **Use float64 for parity with the original
    solver**; float32 drifts by ~1e-7 per simulated day, and Apple MPS supports
    float32 only.

    When ``nn_module`` is provided it receives ``(1, len(nn_channels), ny, nx)``
    and must return ``(1, 2, ny, nx)``; the two output channels correct
    geostrophic ``u`` and ``v`` before the PV RHS is evaluated, on exactly the
    faces ``qrhs`` reads for active (``mask == 2``) cells.

    ``nn_channels`` defaults to ``("ssh",)``, matching the QG-Net of
    arXiv:1911.08856, whose ConvNet "takes SSH as input". The other available
    channels are ``"mask"`` (``grd.mask / 2``: land 0, boundary 0.5, interior
    1) and ``"f0"``, ``"dx"``, ``"dy"`` (standardised over active interior
    cells, zero elsewhere), e.g.::

        nn_channels=("ssh", "mask")                      # land-aware
        nn_channels=("ssh", "mask", "f0", "dx", "dy")    # multi-region

    Be deliberate about the extra channels: they are pure functions of the
    grid, so on a single fixed domain they are constant across every sample and
    act as a positional encoding, which breaks the translation equivariance of
    the correction and makes it non-transferable to another region.

    ``nn_cadence_steps`` controls how often the CNN is evaluated. A value of 1
    evaluates it every internal model step; larger values evaluate it every N
    steps and hold the latest correction in between. The held correction keeps
    its autograd graph, so gradients still flow to every step it covers.

    Note:
        The returned ``SSH`` buffer is allocated with ``torch.empty`` (as the
        original used ``numpy.empty``). Every frame is written for any
        configuration where ``dtout`` is a multiple of ``dt``.
    """
    for name, value in (("Hi", Hi), ("c", c), ("lon", lon), ("lat", lat), ("tint", tint), ("dtout", dtout), ("dt", dt)):
        _required(name, value)

    dtype, device = infer_dtype_device(Hi, c, lon, lat, Hm)
    Hi = as_tensor(Hi, dtype=dtype, device=device)
    if Hi.ndim not in (2, 3):
        raise ValueError(f"Hi must be (ny, nx) or (B, ny, nx), got shape {tuple(Hi.shape)}")
    batched = Hi.ndim == 3
    if batched:
        # One grid serves the whole batch, and the grid's land mask comes from
        # the NaN pattern of SSH, so every sample must share it.
        if not bool((torch.isnan(Hi) == torch.isnan(Hi[0])).all()):
            raise ValueError("every sample of a batched Hi must have the same land (NaN) pattern")
        if obsspace is not None:
            raise NotImplementedError("obsspace sampling is not supported for a batched Hi")
    Hi2d = Hi[0] if batched else Hi
    c = as_tensor(c, dtype=dtype, device=device)
    lon = as_tensor(lon, dtype=dtype, device=device)
    lat = as_tensor(lat, dtype=dtype, device=device)
    obsspace = as_tensor(obsspace, dtype=dtype, device=device) if obsspace is not None else None
    Hm = as_tensor(Hm, dtype=dtype, device=device) if Hm is not None else None
    rappel = maybe_tensor_scalar(rappel, Hi)

    way = sign_scalar(tint)
    tint_abs = abs(scalar_float(tint))
    dtout = scalar_float(dtout)
    dt = scalar_float(dt)
    if dt <= 0.0 or dtout <= 0.0:
        raise ValueError("dt and dtout must be positive")
    nn_cadence_steps = int(nn_cadence_steps)
    if nn_cadence_steps < 1:
        raise ValueError("nn_cadence_steps must be a positive integer")

    if grd is None:
        grd = modgrid.grid(Hi2d, c, snu, lon, lat)
    else:
        # Every cached quantity on the grid depends on these arguments, so a
        # mismatched cache silently changes the physics rather than failing.
        grd.check_matches(Hi2d, c, snu, lon, lat)
    time_abs = 0.0
    index_time = 0
    if obsspace is not None:
        hg = torch.full((obsspace.shape[0],), torch.nan, dtype=dtype, device=device)
        _fill_observations(hg, _observation_indices(obsspace, way, time_abs, dt), Hi, lon, lat, obsspace)
    else:
        hg = None

    nindex_time = int(tint_abs / dtout) + 1
    SSH = torch.empty(Hi.shape[:-2] + (nindex_time, grd.ny, grd.nx), dtype=dtype, device=device)
    SSH[..., index_time, :, :] = Hi

    nstep = int(tint_abs / dt)
    stepout = int(dtout / dt)
    if stepout <= 0:
        raise ValueError("dtout must be greater than or equal to dt")

    if rappel is not None:
        if Hm is None:
            raise ValueError("Hm is required when rappel is not None")
        if tuple(Hm.shape) not in ((grd.ny, grd.nx), tuple(Hi.shape)):
            raise ValueError(f"Hm has shape {tuple(Hm.shape)} but Hi has {tuple(Hi.shape)}")
        # A NaN in Hm over water reaches the advection stencil through Qm and
        # turns the whole field NaN within a few steps, without any error.
        if bool(torch.isnan(Hm)[..., grd.mask >= 1.0].any()):
            raise ValueError(
                "Hm is NaN on cells where Hi is water; its land mask must not "
                "extend beyond that of Hi"
            )
        Qm, = modelliptic.h2pv(Hm, grd)

    h = Hi.clone()
    q, = modelliptic.h2pv(h, grd)
    hb = h.clone()
    u_corr = None
    v_corr = None
    if nn_module is not None:
        # The grid channels and the valid-support masks never change during a
        # run; building them once keeps the step loop allocation-free.
        nn_spec = normalize_channel_spec(nn_channels)
        _check_nn_in_channels(nn_module, len(nn_spec))
        static_channels = _build_static_channels(grd, nn_spec)
        valid_u = _neural_u_valid_mask(grd)
        valid_v = _neural_v_valid_mask(grd)

    for step in range(nstep):
        time_abs = (step + 1) * dt
        if (step + 1) % stepout == 0:
            index_time += 1

        hguess = 2.0 * h - hb
        hb = h.clone()
        qb = q.clone()

        u, v = moddyn.h2uv(h, grd)
        if nn_module is not None:
            if step % nn_cadence_steps == 0:
                u_corr, v_corr = _compute_neural_velocity_correction(
                    _neural_input(h, nn_spec, static_channels), grd, nn_module,
                    batched=batched, checkpoint=nn_checkpoint,
                )
            u, v = _apply_neural_velocity_correction(u, v, u_corr, v_corr, valid_u, valid_v)
        rq, = moddyn.qrhs(u, v, qb, grd, way)
        if rappel is not None:
            q = qb + dt * (rq - rappel * (qb - Qm))
        else:
            q = qb + dt * rq
        h, = modelliptic.pv2h(q, hguess, grd)

        if (step + 1) % stepout == 0 and index_time < nindex_time:
            SSH[..., index_time, :, :] = h

        if obsspace is not None:
            _fill_observations(hg, _observation_indices(obsspace, way, time_abs, dt), h, lon, lat, obsspace)

    return SSH, hg
