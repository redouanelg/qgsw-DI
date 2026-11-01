"""Grid construction for the PyTorch QGSW model."""

from __future__ import annotations

import math

import torch

try:
    from ._torch_utils import as_tensor, infer_dtype_device
except ImportError:  # pragma: no cover - supports running from this folder directly
    from _torch_utils import as_tensor, infer_dtype_device


class grid:
    """Torch equivalent of the original ``modgrid.grid`` class.

    The public attributes intentionally keep the original names so the dynamics
    and elliptic inversion can share the same structure as the NumPy code.

    Attributes:
        mask: ``2`` active interior, ``1`` Dirichlet boundary, ``0`` land.
            The border ring is ASYMMETRIC -- 2 cells at the low edges, 3 at the
            high edges -- because the advection stencil is 2 cells wide and the
            original was written that way. Every interior NaN additionally
            flags a 5x5 halo as ``1``; NaNs sitting exactly on the outermost
            row/column do not expand (the original loop ran over
            ``1..ny-2`` / ``1..nx-2``).
        dx, dy: forward differences in metres, edges replicated. ``dy`` carries
            the same ``cos(lat)`` factor as ``dx``; that is how the original
            computes it and it is reproduced here for bit parity.
        f0: ``2*Omega*sin(lat)`` with ``Omega = 2*pi/86164`` (sidereal day).
            Division by ``f0`` is singular at the equator.
        snu: pass ``0.`` rather than ``None`` to activate the tendency masking
            in :func:`moddyn.qrhs` (see the warning there).
        vp2nn/ss/ee/ww/nw/ne/se/sw: unused by the solver, kept so the attribute
            set matches the NumPy ``grid`` object.
    """

    def __init__(self, SSH, c, snu, lon, lat, dtype=None, device=None):
        dtype, device = infer_dtype_device(SSH, c, lon, lat, dtype=dtype, device=device)
        SSH = as_tensor(SSH, dtype=dtype, device=device)
        lon = as_tensor(lon, dtype=dtype, device=device)
        lat = as_tensor(lat, dtype=dtype, device=device)
        c = as_tensor(c, dtype=dtype, device=device)
        if c.ndim == 0:
            c = torch.zeros_like(SSH) + c
        if SSH.ndim != 2:
            raise ValueError("SSH must be a 2D tensor or array")
        if c.shape != SSH.shape or lon.shape != SSH.shape or lat.shape != SSH.shape:
            raise ValueError("SSH, c, lon and lat must have matching 2D shapes")

        ny, nx = SSH.shape
        mask = torch.full((ny, nx), 2.0, dtype=dtype, device=device)
        mask[:2, :] = 1.0
        mask[:, :2] = 1.0
        mask[-3:, :] = 1.0
        mask[:, -3:] = 1.0

        dx = torch.zeros((ny, nx), dtype=dtype, device=device)
        dy = torch.zeros((ny, nx), dtype=dtype, device=device)
        if ny > 2 and nx > 2:
            lat_mid = lat[1:-1, 1:-1]
            dlonx = lon[1:-1, 2:] - lon[1:-1, 1:-1]
            dlatx = lat[1:-1, 2:] - lat[1:-1, 1:-1]
            dlony = lon[2:, 1:-1] - lon[1:-1, 1:-1]
            dlaty = lat[2:, 1:-1] - lat[1:-1, 1:-1]
            coslat = torch.cos(lat_mid * math.pi / 180.0)
            dx[1:-1, 1:-1] = torch.sqrt((dlonx * 111000.0 * coslat) ** 2 + (dlatx * 111000.0) ** 2)
            dy[1:-1, 1:-1] = torch.sqrt((dlony * 111000.0 * coslat) ** 2 + (dlaty * 111000.0) ** 2)

        # Flag a 5x5 halo around every NaN that the original double loop would
        # have visited, i.e. excluding NaNs sitting on the outermost row or
        # column. A stride-1 max pool with padding 2 is exactly that dilation,
        # clipped to the domain, and avoids the original's per-NaN host sync.
        nan_ssh = torch.isnan(SSH)
        if ny > 2 and nx > 2:
            interior_nan = torch.zeros((ny, nx), dtype=dtype, device=device)
            interior_nan[1:-1, 1:-1] = nan_ssh[1:-1, 1:-1].to(dtype)
            if bool(interior_nan.any()):
                halo = torch.nn.functional.max_pool2d(
                    interior_nan.unsqueeze(0).unsqueeze(0), kernel_size=5, stride=1, padding=2
                )[0, 0]
                mask = torch.where(halo > 0, torch.ones_like(mask), mask)

        if ny > 1:
            dx[0, :] = dx[1, :]
            dx[-1, :] = dx[-2, :]
            dy[0, :] = dy[1, :]
            dy[-1, :] = dy[-2, :]
        if nx > 1:
            dx[:, 0] = dx[:, 1]
            dx[:, -1] = dx[:, -2]
            dy[:, 0] = dy[:, 1]
            dy[:, -1] = dy[:, -2]

        mask[nan_ssh] = 0.0
        f0 = 4.0 * math.pi / 86164.0 * torch.sin(lat * math.pi / 180.0)

        valid = torch.where(mask >= 1.0)
        interior = torch.where(mask[valid] == 2.0)[0]
        boundary = torch.where(mask[valid] == 1.0)[0]
        npoints = int(valid[0].numel())

        self.mask1d = mask[valid]
        self.H = SSH[valid]
        self.c1d = c[valid]
        self.f01d = f0[valid]
        self.dx1d = dx[valid]
        self.dy1d = dy[valid]
        self.indi = valid[0].long()
        self.indj = valid[1].long()
        self.vp1 = boundary.long()
        self.vp2 = interior.long()
        self.indp = torch.zeros((ny, nx), dtype=torch.long, device=device)
        self.indp[valid] = torch.arange(npoints, dtype=torch.long, device=device)

        i2 = self.indi[self.vp2]
        j2 = self.indj[self.vp2]
        self.vp2n = self.indp[i2 + 1, j2]
        self.vp2nn = self.indp[i2 + 2, j2]
        self.vp2s = self.indp[i2 - 1, j2]
        self.vp2ss = self.indp[i2 - 2, j2]
        self.vp2e = self.indp[i2, j2 + 1]
        self.vp2ee = self.indp[i2, j2 + 2]
        self.vp2w = self.indp[i2, j2 - 1]
        self.vp2ww = self.indp[i2, j2 - 2]
        self.vp2nw = self.indp[i2 + 1, j2 - 1]
        self.vp2ne = self.indp[i2 + 1, j2 + 1]
        self.vp2se = self.indp[i2 - 1, j2 + 1]
        self.vp2sw = self.indp[i2 - 1, j2 - 1]

        self.mask = mask
        self.f0 = f0
        self.dx = dx
        self.dy = dy
        self.nan_ssh = nan_ssh
        self.np = npoints
        self.np2 = int(self.vp2.numel())
        self.nx = nx
        self.ny = ny
        self.c = c
        self.snu = None if snu is None else as_tensor(snu, dtype=dtype, device=device)
        self.lon = lon
        self.lat = lat
        self.g = torch.as_tensor(9.81, dtype=dtype, device=device)

        # ---- derived constants, cached so the time loop stays allocation-light.
        # These are pure functions of the grid; recomputing them every step (as
        # the original does) costs ~15% of the runtime and nothing else.
        self.mask_boundary = mask == 1.0
        self.mask_land = mask == 0.0
        self.mask_inactive = mask <= 1.0
        self.stretching = -self.g * f0 / (c**2)  # q = ... + stretching * h

        # Elliptic operator coefficients used by modelliptic.compute_avec.
        aaa = self.g / self.f01d
        bbb = -self.g * self.f01d / self.c1d**2
        if self.vp1.numel() > 0:
            aaa[self.vp1] = 0.0
            bbb[self.vp1] = 1.0
        self.aaa = aaa
        self.bbb = bbb
        # NOTE: cache the squares, not their reciprocals. ``a / b`` and
        # ``a * (1/b)`` are not the same in floating point and the whole point
        # of this port is bit parity with the NumPy model.
        self.dx1d2 = self.dx1d**2
        self.dy1d2 = self.dy1d**2
        self.aaa_vp2 = aaa[self.vp2]
        self.bbb_vp2 = bbb[self.vp2]
        self.dx1d2_vp2 = self.dx1d2[self.vp2]
        self.dy1d2_vp2 = self.dy1d2[self.vp2]

        # Constants of the PV tendency stencil (moddyn.qrhs). All of these are
        # grid-only; factoring them out of the time loop is bit-neutral because
        # the surviving operation order is unchanged.
        self.g_over_f0 = self.g / f0
        if ny > 4 and nx > 4:
            dx_int = dx[2:-2, 2:-2]
            dy_int = dy[2:-2, 2:-2]
            self.six_dx_int = 6.0 * dx_int
            self.six_dy_int = 6.0 * dy_int
            # (f0[i+1] - f0[i-1]) / (2*dy): the beta term coefficient.
            self.beta_int = (f0[3:-1, 2:-2] - f0[1:-3, 2:-2]) / (2.0 * dy_int)
            if self.snu is not None:
                self.snu_over_dx2_int = self.snu / (dx_int**2)
                self.snu_over_dy2_int = self.snu / (dy_int**2)
            else:
                self.snu_over_dx2_int = None
                self.snu_over_dy2_int = None
        else:
            self.six_dx_int = None
            self.six_dy_int = None
            self.beta_int = None
            self.snu_over_dx2_int = None
            self.snu_over_dy2_int = None

    def check_matches(self, SSH, c, snu, lon, lat):
        """Raise if this grid was not built from equivalent inputs.

        Used by :func:`qgsw.qgsw` when a caller supplies a prebuilt grid. Every
        cached quantity on this object (``mask``, ``dx``, ``dy``, ``f0``,
        ``aaa``, ``bbb``, ``snu_over_dx2_int``, ...) is a pure function of these
        five arguments, so reusing a grid built from different ones silently
        changes the physics. The cost is O(ny*nx), negligible against a run of
        many time steps.
        """

        def same(a, b):
            if a is None or b is None:
                return (a is None) and (b is None)
            if a.shape != b.shape:
                return False
            return bool(
                torch.equal(torch.isnan(a), torch.isnan(b))
                and torch.equal(torch.nan_to_num(a, nan=0.0), torch.nan_to_num(b, nan=0.0))
            )

        SSH_t = as_tensor(SSH, like=self.f0)
        if tuple(SSH_t.shape) != (self.ny, self.nx):
            raise ValueError(
                f"supplied grd has shape ({self.ny}, {self.nx}) but Hi has {tuple(SSH_t.shape)}"
            )

        c_t = as_tensor(c, like=self.f0)
        if c_t.ndim == 0:
            c_t = torch.zeros_like(SSH_t) + c_t
        snu_t = None if snu is None else as_tensor(snu, like=self.f0)

        mismatches = []
        # The land mask is derived from the NaN pattern of SSH, not its values.
        if not torch.equal(torch.isnan(SSH_t), self.nan_ssh):
            mismatches.append("Hi land/NaN pattern")
        if not same(c_t, self.c):
            mismatches.append("c")
        if not same(as_tensor(lon, like=self.f0), self.lon):
            mismatches.append("lon")
        if not same(as_tensor(lat, like=self.f0), self.lat):
            mismatches.append("lat")
        if (snu_t is None) != (self.snu is None) or (
            snu_t is not None and not bool(torch.equal(snu_t, self.snu))
        ):
            mismatches.append(f"snu ({snu} vs {self.snu})")

        if mismatches:
            raise ValueError(
                "the supplied grd was built from different inputs than this call: "
                + ", ".join(mismatches)
                + ". Rebuild the grid, or drop the grd= argument."
            )


Grid = grid
