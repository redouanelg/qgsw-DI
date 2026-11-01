"""Dynamical QGSW operators implemented with PyTorch tensors."""

from __future__ import annotations

import torch

try:
    from ._torch_utils import as_tensor, nan_to_zero
except ImportError:  # pragma: no cover
    from _torch_utils import as_tensor, nan_to_zero


def h2uv(h, grd):
    """Convert SSH to geostrophic velocity components.

    ``u = -(g/f0) dh/dy`` and ``v = +(g/f0) dh/dx``, discretised exactly as the
    original ``moddyn.h2uv``: 4-point corner averages differenced across one
    cell, which reduces to

        u[i,j] = -(g/f0) * 0.25*(h[i+1,j-1]+h[i+1,j]-h[i-1,j-1]-h[i-1,j]) / dy
        v[i,j] = +(g/f0) * 0.25*(h[i,j+1]+h[i-1,j+1]-h[i,j-1]-h[i-1,j-1]) / dx

    so ``u`` lives on the western face and ``v`` on the southern face of the
    cell (Arakawa C-grid staggering). The first/last row and column are left at
    zero. NaNs coming from land neighbours are zeroed rather than masked, which
    is also what the original does.
    """
    h = as_tensor(h, like=grd.f0)
    u = torch.zeros((grd.ny, grd.nx), dtype=h.dtype, device=h.device)
    v = torch.zeros((grd.ny, grd.nx), dtype=h.dtype, device=h.device)
    g = grd.g

    # grd.g_over_f0 == g / f0, precomputed. IEEE negation is exact, so
    # (-g)/f0 and -(g/f0) are the same bit pattern.
    u[1:-1, 1:] = (
        -grd.g_over_f0[1:-1, 1:]
        * (
            0.25 * (h[1:-1, 1:] + h[1:-1, :-1] + h[2:, :-1] + h[2:, 1:])
            - 0.25 * (h[1:-1, 1:] + h[:-2, 1:] + h[:-2, :-1] + h[1:-1, :-1])
        )
        / grd.dy[1:-1, 1:]
    )
    v[1:, 1:-1] = (
        grd.g_over_f0[1:, 1:-1]
        * (
            0.25 * (h[1:, 1:-1] + h[1:, 2:] + h[:-1, 2:] + h[:-1, 1:-1])
            - 0.25 * (h[1:, 1:-1] + h[:-1, 1:-1] + h[:-1, :-2] + h[1:, :-2])
        )
        / grd.dx[1:, 1:-1]
    )
    return nan_to_zero(u), nan_to_zero(v)


def _split_flux(a, way, minus_strict=True):
    """Upwind split of a face-averaged velocity, following the original exactly.

    ``uminus`` is zeroed where ``> 0`` whereas ``vminus`` is zeroed where
    ``>= 0``. The asymmetry is inherited from the NumPy code and is numerically
    inert, but it is preserved so the two implementations agree bit for bit.

    Returns ``(plus, minus, scaled)`` where ``scaled = way * 0.5 * a`` is the
    unsplit value, reused by the beta term instead of being recomputed.

    Note:
        Reusing ``scaled`` for the beta term, and folding ``way * 0.5`` into it
        rather than into the beta coefficient, is bit-exact ONLY because
        ``way`` is +-1 and ``0.5`` is a power of two: both scalings are then
        exact in IEEE754 and the reassociation cannot change any rounding.
        ``way`` comes from ``_torch_utils.sign_scalar``, which can only return
        +-1.0 or 0.0. An arbitrary ``way`` would break bit parity with the
        NumPy model (measured: ``way=3`` differs on ~18% of cells).
    """
    scaled = way * 0.5 * a
    zero = torch.zeros_like(scaled)
    plus = torch.where(scaled < 0.0, zero, scaled)
    minus = torch.where(scaled > 0.0 if minus_strict else scaled >= 0.0, zero, scaled)
    return plus, minus, scaled


def qrhs(u, v, q, grd, way):
    """Compute the PV tendency ``dq/dt``.

    Third-order upwind-biased (QUICK-type) advection of ``q`` by the C-grid
    velocities, plus an explicit ``-beta*v`` term (``beta`` is obtained by
    centred differencing of ``f0``, because ``q`` excludes planetary vorticity),
    plus optional harmonic diffusion ``snu * laplacian(q)``.

    Warning:
        ``snu=0.`` and ``snu=None`` are NOT equivalent. The masking
        ``rq[mask <= 1] = 0`` lives inside the ``snu is not None`` branch, so
        with ``snu=None`` the tendency is left non-zero over the coastal halo.
        This reproduces the original ``moddyn.qrhs``; pass ``snu=0.`` for the
        published configuration.
    """
    u = as_tensor(u, like=grd.f0)
    v = as_tensor(v, like=grd.f0)
    q = as_tensor(q, like=grd.f0)
    rq = torch.zeros((grd.ny, grd.nx), dtype=q.dtype, device=q.device)
    if grd.ny <= 4 or grd.nx <= 4:
        return (rq,)

    uplus, uminus, _ = _split_flux(u[2:-2, 2:-2] + u[2:-2, 3:-1], way, minus_strict=True)
    vplus, vminus, vscaled = _split_flux(v[2:-2, 2:-2] + v[3:-1, 2:-2], way, minus_strict=False)

    # grd.six_dx_int == 6*dx, grd.six_dy_int == 6*dy, grd.beta_int ==
    # (f0[i+1]-f0[i-1])/(2*dy); all precomputed by modgrid.grid.
    six_dx = grd.six_dx_int
    six_dy = grd.six_dy_int
    qc = q[2:-2, 2:-2]

    interior = (
        rq[2:-2, 2:-2]  # zeros; kept so the operation order matches the original
        - uplus / six_dx * (2.0 * q[2:-2, 3:-1] + 3.0 * qc - 6.0 * q[2:-2, 1:-3] + q[2:-2, :-4])
        + uminus / six_dx * (q[2:-2, 4:] - 6.0 * q[2:-2, 3:-1] + 3.0 * qc + 2.0 * q[2:-2, 1:-3])
        - vplus / six_dy * (2.0 * q[3:-1, 2:-2] + 3.0 * qc - 6.0 * q[1:-3, 2:-2] + q[:-4, 2:-2])
        + vminus / six_dy * (q[4:, 2:-2] - 6.0 * q[3:-1, 2:-2] + 3.0 * qc + 2.0 * q[1:-3, 2:-2])
    )
    # vscaled is exactly `way * 0.5 * (v[2:-2,2:-2] + v[3:-1,2:-2])`.
    interior = interior - grd.beta_int * vscaled

    if grd.snu is not None:
        interior = (
            interior
            + grd.snu_over_dx2_int * (q[2:-2, 3:-1] + q[2:-2, 1:-3] - 2.0 * qc)
            + grd.snu_over_dy2_int * (q[3:-1, 2:-2] + q[1:-3, 2:-2] - 2.0 * qc)
        )
        rq[2:-2, 2:-2] = interior
        rq[grd.mask_inactive] = 0.0
    else:
        rq[2:-2, 2:-2] = interior
    return (rq,)
