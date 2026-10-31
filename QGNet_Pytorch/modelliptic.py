"""Elliptic SSH/PV transforms implemented with PyTorch tensors."""

from __future__ import annotations

import torch

try:
    from ._torch_utils import as_tensor
except ImportError:  # pragma: no cover
    from _torch_utils import as_tensor


def _safe_divide(numerator, denominator, fallback=1.0):
    """Branchless divide that falls back to ``fallback`` when ``denominator`` is zero.

    Mirrors the original NumPy ``if val2 == 0.: s = 1. else: s = val1 / val2``
    without the ``float(...detach().cpu())`` host-device sync in the CG loop.

    A NaN denominator must NOT be treated as a zero: in the original that case
    propagates NaN into ``h`` and makes a diverging solve immediately visible.
    ``torch.abs(nan) > 0`` is False, so the zero test is written as ``== 0`` and
    the NaN/Inf case is routed through the ordinary division.
    """
    is_zero = denominator == 0
    denom_safe = torch.where(is_zero, torch.ones_like(denominator), denominator)
    fallback_value = torch.full_like(denominator, float(fallback))
    return torch.where(is_zero, fallback_value, numerator / denom_safe)


def h2pv(h, grd):
    """Convert SSH to potential vorticity.

    ``q = (g/f0) * laplacian(h) - (g*f0/c**2) * h`` on the interior, the pure
    stretching term alone on the boundary ring (``mask == 1``), and zero on land
    (``mask == 0``).

    Note:
        Planetary vorticity ``f0`` is deliberately NOT included in ``q``; the
        beta effect enters as a separate explicit term in ``moddyn.qrhs``.
        The Laplacian uses locally frozen ``dx[i,j]``/``dy[i,j]`` with no
        metric-term correction, so the operator is only approximately
        self-adjoint on a non-uniform grid.
    """
    h = as_tensor(h, like=grd.f0)
    g = grd.g
    # grd.stretching == -g * f0 / c**2, precomputed; same op order as the
    # original so the result is bit-identical.
    q = grd.stretching * h

    q[1:-1, 1:-1] = (
        g
        / grd.f0[1:-1, 1:-1]
        * (
            (h[2:, 1:-1] + h[:-2, 1:-1] - 2.0 * h[1:-1, 1:-1]) / grd.dy[1:-1, 1:-1] ** 2
            + (h[1:-1, 2:] + h[1:-1, :-2] - 2.0 * h[1:-1, 1:-1]) / grd.dx[1:-1, 1:-1] ** 2
        )
        + grd.stretching[1:-1, 1:-1] * h[1:-1, 1:-1]
    )
    boundary = grd.mask_boundary
    q[boundary] = grd.stretching[boundary] * h[boundary]
    q[grd.mask_land] = 0.0
    return (q,)


def pv2h(q, hg, grd, nitr=1):
    """Invert potential vorticity to SSH with conjugate-gradient iterations.

    Args:
        q: potential vorticity field.
        hg: SSH guess. It also carries the Dirichlet boundary values.
        grd: :class:`modgrid.grid`.
        nitr: TOTAL number of CG steps. The default of 1 performs a single
            optimally-scaled steepest-descent step; the solve is deliberately
            under-converged and relies on the time-extrapolated guess.

    Returns:
        A 1-tuple holding the SSH field. Land points (``mask == 0``) are NaN.

    Note:
        Boundary rows (``mask == 1``) are turned into identity rows whose
        right-hand side is taken from ``hg``, i.e. ``h = hg`` there. Combined
        with ``hguess = 2*h - hb`` in :func:`qgsw.qgsw`, the boundary SSH stays
        frozen at the initial condition for the whole run.
        The operator is not symmetric, so "CG" is nominal; with ``nitr=1`` that
        is irrelevant.
    """
    q = as_tensor(q, like=grd.f0)
    hg = as_tensor(hg, like=q)

    x = hg[grd.indi, grd.indj].clone()
    q1d = q[grd.indi, grd.indj]

    # aaa/bbb depend only on the grid and are cached by modgrid.grid.
    aaa = grd.aaa
    bbb = grd.bbb
    ccc = q1d.clone()

    if grd.vp1.numel() > 0:
        ccc[grd.vp1] = x[grd.vp1]

    avec, = compute_avec(x, aaa, bbb, grd)
    gg = avec - ccc
    p = -gg

    for _ in range(max(int(nitr) - 1, 0)):
        avec, = compute_avec(p, aaa, bbb, grd)
        tmp = torch.dot(p, avec)
        s = _safe_divide(-torch.dot(p, gg), tmp)

        a1 = torch.dot(gg, gg)
        x = x + s * p
        avec, = compute_avec(x, aaa, bbb, grd)
        gg = avec - ccc
        a2 = torch.dot(gg, gg)
        beta = _safe_divide(a2, a1)
        p = -gg + beta * p

    avec, = compute_avec(p, aaa, bbb, grd)
    val1 = -torch.dot(p, gg)
    val2 = torch.dot(p, avec)
    s = _safe_divide(val1, val2)
    x = x + s * p

    h = torch.full_like(hg, torch.nan)
    h[grd.indi, grd.indj] = x
    return (h,)


def compute_avec(vec, aaa, bbb, grd):
    """Apply the Helmholtz operator on the compressed 1-D vector.

    ``aaa``/``bbb`` are accepted for signature compatibility with the original;
    when they are the grid's own cached coefficients the pre-gathered
    ``grd.*_vp2`` views are used to skip four index_select calls per step.
    """
    vec = as_tensor(vec, like=aaa)
    avec = torch.empty((grd.np,), dtype=vec.dtype, device=vec.device)
    if grd.vp2.numel() > 0:
        if aaa is grd.aaa and bbb is grd.bbb:
            aaa_i, bbb_i = grd.aaa_vp2, grd.bbb_vp2
            dx2_i, dy2_i = grd.dx1d2_vp2, grd.dy1d2_vp2
        else:
            aaa_i, bbb_i = aaa[grd.vp2], bbb[grd.vp2]
            dx2_i, dy2_i = grd.dx1d[grd.vp2] ** 2, grd.dy1d[grd.vp2] ** 2
        vec_c = vec[grd.vp2]
        avec[grd.vp2] = (
            aaa_i
            * (
                (vec[grd.vp2e] + vec[grd.vp2w] - 2.0 * vec_c) / dx2_i
                + (vec[grd.vp2n] + vec[grd.vp2s] - 2.0 * vec_c) / dy2_i
            )
            + bbb_i * vec_c
        )
    if grd.vp1.numel() > 0:
        avec[grd.vp1] = vec[grd.vp1]
    return (avec,)
