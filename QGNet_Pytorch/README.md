# QGNet_Pytorch

PyTorch port of the [`qgsw-DI`](https://github.com/redouanelg/qgsw-DI)
quasigeostrophic shallow-water model used for dynamic interpolation of SSH, with
an optional neural correction on the geostrophic velocity (the QGNet idea of
[arXiv:1911.08856](https://arxiv.org/abs/1911.08856)).

This folder lives inside the `qgsw-DI` repository; the NumPy reference model it
is validated against is one directory up.

The port is **bit-faithful**: on `ssh.nc` in float64 it agrees with the NumPy
original to `max|Δ| ≈ 5e-16` after 5 days and `1.3e-15` after 30 days
(4320 steps), and every grid attribute matches except `dx`, which differs by a
couple of ULP because the vectorised `torch.cos` and the scalar `math.cos`
round differently. See `tests/test_parity.py`.

## Module layout

| file | contents |
|---|---|
| `modgrid.py` | `grid`: metrics, Coriolis, land/boundary mask, 1-D index tables |
| `modelliptic.py` | `h2pv` (SSH → PV), `pv2h` (PV inversion), `compute_avec` |
| `moddyn.py` | `h2uv` (geostrophy), `qrhs` (PV tendency) |
| `qgsw.py` | forward propagator + neural-correction plumbing |
| `_torch_utils.py` | dtype/device helpers |

The public function and return convention match the NumPy version. Inputs may be
NumPy arrays or torch tensors; outputs are torch tensors on the inferred device.

The explicit tangent-linear and adjoint models are intentionally **not** ported:
use autograd. This is strictly better, not just more convenient — the original
`moddyn.qrhs_adj` is only an exact transpose of `qrhs_tgl` when the caller
projects input and output onto `mask == 2` (1.9e-2 relative mismatch without the
projection), and its diffusion term is inexact for variable `dx`. Autograd has
no such caveat and passes `torch.autograd.gradcheck`.

## Install

```bash
pip install -e .            # or: pip install -r requirements.txt
pip install -e '.[test]'    # adds pytest + netCDF4 for the parity suite
```

`scipy` is a hard dependency (used by the `obsspace` sampling path), not optional.

## Usage

```python
import torch
from QGNet_Pytorch import qgsw          # the submodule

Hi    = torch.as_tensor(SSH, dtype=torch.float64)   # NaN marks land
c_t   = torch.as_tensor(c,   dtype=torch.float64)
lon_t = torch.as_tensor(lon, dtype=torch.float64)
lat_t = torch.as_tensor(lat, dtype=torch.float64)

Hf, obs = qgsw.qgsw(
    Hi=Hi, c=c_t, lon=lon_t, lat=lat_t,
    tint=86400 * 5,     # negative integrates backwards
    dtout=86400,
    dt=600,
    rappel=None,
    snu=0.0,            # NOT None -- see "Gotchas"
)
```

`from QGNet_Pytorch import qgsw` gives the **module**; the function is also
exported as `qgsw_fn` for convenience.

### Precision and devices

| device | dtype | note |
|---|---|---|
| CPU | `float64` | reference; parity with NumPy at ~5e-16 |
| CPU | `float32` | ~2.7e-7 drift after one simulated day |
| CUDA | `float64` | full parity, move every input tensor to the device first |
| MPS (Apple) | `float32` only | PyTorch has no float64 on MPS; also measured *slower* than CPU at 113×170 |

Mixing float32 and float64 inputs now promotes to the highest precision and
warns, instead of silently demoting the run.

## Optional neural correction

`qgsw.qgsw(..., nn_module=net)` evaluates `net` inside the time loop. It receives
a tensor of shape `(1, 5, ny, nx)`:

| channel | content |
|---|---|
| 0 | SSH with NaN replaced by 0 |
| 1 | mask / 2 → land `0`, boundary `0.5`, interior water `1` |
| 2 | normalised `f0` |
| 3 | normalised `dx` |
| 4 | normalised `dy` |

Channels 2–4 are standardised over active interior-water cells and zeroed
elsewhere; they are constant in time and built once per run. Longitude and
latitude are not passed directly.

The module must return `(1, 2, ny, nx)`. Channel 0 is added to geostrophic `u`
and channel 1 to `v`, restricted via `torch.where` to the exact discrete
supports `qrhs` reads (`u[2:-2, 2:-1]`, `v[2:-1, 2:-2]`, interior water only), so
a NaN emitted outside the support cannot leak into the solution.

`nn_cadence_steps` controls the evaluation frequency. With `dt=600`,
`nn_cadence_steps=6` gives an hourly correction update while the dynamics still
advance every 10 minutes. The held correction keeps its autograd graph, so
gradients reach every step it covers — memory therefore grows with rollout
length regardless of cadence.

A zero-output network reproduces the pure-QG trajectory exactly (`max|Δ| = 0`).

## Gotchas inherited from the original (deliberately preserved)

* **`snu=0.` is not `snu=None`.** The masking `rq[mask <= 1] = 0` lives inside
  the `if grd.snu is not None` branch of `qrhs`. Pass `snu=0.` for the published
  configuration; `snu=None` leaves a non-zero tendency over the coastal halo
  (measured `max|Δh| = 4.5e-4 m` after one hour).
* **`pv2h` defaults to `nitr=1`**: a single optimally-scaled steepest-descent
  step. The solve is deliberately under-converged and relies on the
  time-extrapolated guess `hguess = 2hⁿ − hⁿ⁻¹`.
* **Boundary condition.** `mask == 1` rows are identity rows whose RHS is the
  guess, so boundary SSH stays frozen at the initial condition for the whole run.
* **The border ring is asymmetric**: 2 cells low, 3 cells high.
* **`q` excludes planetary vorticity**; beta enters as a separate explicit term.
* **`dy` carries the same `cos(lat)` factor as `dx`.**
* `dtout` should be an integer multiple of `dt`, otherwise output frames are not
  evenly spaced. Unlike the original, degenerate configurations (`dtout < dt`,
  `dt <= 0`) raise instead of silently returning uninitialised memory.

## Tests

```bash
pytest                                  # from this folder
QGSW_DI_PATH=/path/to/qgsw-DI pytest    # if the NumPy model is elsewhere
```

187 tests. `tests/test_parity.py` cross-checks every operator and whole
integrations (forward, backward, nudged, `snu=None`/`0`/`50`, `nitr` 1…25, with
and without land) against the NumPy model in the parent directory; those tests
skip if it cannot be found. `tests/test_regression.py` pins the float64 output
byte-for-byte, which is what makes it safe to optimise this solver.

The shared test case `ssh.nc` is the one at the repository root.

## Example notebooks

* `notebooks/validation_vs_numpy.ipynb` — reproduces the NumPy↔PyTorch parity
  figures and error metrics on `ssh.nc`.
* `notebooks/SSHforecastingQGSW_PyTorch_multistep.ipynb` — ResUNet velocity
  correction trained on NATL60 patches with a multi-day rollout curriculum.

## Notes

* Observation-space sampling follows the original `scipy.interpolate.griddata`
  path when `obsspace` is given. That path is CPU-only and non-differentiable.
  The original crashes here (`UnboundLocalError`); the port does not.
* Grid construction follows the original mask logic exactly.
