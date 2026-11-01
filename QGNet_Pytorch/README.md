# QGNet_Pytorch

PyTorch port of the [`qgsw-DI`](https://github.com/redouanelg/qgsw-DI)
quasigeostrophic shallow-water model used for dynamic interpolation of SSH, with
an optional neural correction on the geostrophic velocity (the QGNet idea of
[arXiv:1911.08856](https://arxiv.org/abs/1911.08856)).

This folder lives inside the `qgsw-DI` repository; the NumPy reference model it
is validated against is one directory up.

## Numerical parity

The port reproduces the NumPy original to **machine precision** — not bit for
bit. On `ssh.nc` in float64:

| | max\|Δ\| |
|---|---|
| `h2pv`, `h2uv`, `qrhs`, `pv2h` | ≤ 4.4e-16 |
| 5-day run (720 steps) | 5.6e-16 |
| 30-day run (4320 steps) | 1.3e-15 |

The difference grows about 3× over 4320 steps, i.e. sub-linearly: roundoff does
not amplify. See `tests/test_parity.py`.

The residual is unavoidable, because NumPy and PyTorch do not produce identical
results for the same elementary operations:

- `torch.cos` vs `math.cos` → `dx` differs on ~5% of cells by ≤ 2 ULP,
- `torch.sin` vs `numpy.sin` → `f0` differs on ~3% of cells by 1 ULP,
- `torch.dot` vs `numpy.dot` use different summation orders, which perturbs the
  conjugate-gradient line search in `pv2h` at the last-bit level.

All index tables (`indi`, `indj`, `vp1`, `vp2`, `vp2n/s/e/w`, `indp`) and the
land mask are exactly equal; `dy` happens to be exact on a regular lon/lat mesh
because its `cos(lat)` factor multiplies a zero longitude increment.

What *is* bit-exact is the internal refactoring: `tests/test_regression.py`
pins the float64 output byte-for-byte, so the performance work in `modgrid`,
`moddyn` and `modelliptic` is verified to change no bits at all.

## Module layout

| file | contents |
|---|---|
| `modgrid.py` | `grid`: metrics, Coriolis, land/boundary mask, 1-D index tables |
| `modelliptic.py` | `h2pv` (SSH → PV), `pv2h` (PV inversion), `compute_avec` |
| `moddyn.py` | `h2uv` (geostrophy), `qrhs` (PV tendency) |
| `qgsw.py` | forward propagator + neural-correction plumbing |
| `models.py` | `QGNetConvNet`, the paper's 2545-parameter correction network |
| `_torch_utils.py` | dtype/device helpers |

The public function and return convention match the NumPy version. Inputs may be
NumPy arrays or torch tensors; outputs are torch tensors on the inferred device.

The explicit tangent-linear and adjoint models are intentionally **not** ported:
use autograd. This is strictly better, not just more convenient — the original
`moddyn.qrhs_adj` is only an exact transpose of `qrhs_tgl` when the caller
projects input and output onto `mask == 2` (1.9e-2 relative mismatch without the
projection), and its diffusion term is inexact for variable `dx`. Autograd has
no such caveat and passes `torch.autograd.gradcheck`.

Gradients are finite with respect to the initial SSH, the network parameters,
`snu`, `rappel`, `c`, `lat` and `lon`, with or without land. Land (NaN in SSH)
is zero-filled before the stencils of `h2pv` and `h2uv` and the affected values
are masked, so autograd never multiplies a zero gradient by NaN. This leaves
the forward result bit-identical. To learn `c`, `lat`, `lon` or `snu`, build
the grid every iteration rather than passing a cached `grd=`: its coefficients
carry the graph, and a second `backward()` through them fails.

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
| MPS (Apple) | `float32` only | PyTorch has no float64 on MPS. A single unbatched QG run is slower than on CPU, but a batched hybrid rollout with a ResUNet is ~6× faster (0.2 s vs 1.15 s per sample-day at 128×128 on an M2) and agrees with CPU float32 to 1e-7 m |

Mixing float32 and float64 inputs promotes to the highest precision and warns,
instead of silently demoting the run. This includes NumPy arrays passed next to
tensors. NumPy-only inputs run in float64.

### Batched integration

`Hi` may be a batch `(B, ny, nx)` of independent initial states that share one
land mask (and so one grid). They are integrated together, and `SSH` comes
back as `(B, T, ny, nx)`. `Hm` may be `(ny, nx)` or `(B, ny, nx)`; `obsspace`
is not supported for a batch.

```python
SSH, _ = qgsw.qgsw(Hi=ssh_days[[0, 7, 14, 21]], c=c_t, lon=lon_t, lat=lat_t,
                   tint=3 * 86400, dtout=86400, dt=1200, snu=0.0)   # (4, 4, ny, nx)
```

Each sample gets its own line search in `pv2h`. A sample therefore agrees with
the same sample run alone to round-off (~1e-16 in float64), not bit for bit:
batched dot products sum in a different order. Unbatched runs are unchanged
and still bit-identical to before.

How much batching helps depends on the device:

* on CPU it gives about 1.5× for plain QG and little once a network
  dominates the cost;
* on a GPU (CUDA, or Apple MPS) it is what makes the device worthwhile.

## Optional neural correction

`qgsw.qgsw(..., nn_module=net)` evaluates `net` inside the time loop. It is fed
`(1, len(nn_channels), ny, nx)` and must return `(1, 2, ny, nx)`: channel 0 is
added to geostrophic `u`, channel 1 to `v`, restricted via `torch.where` to
exactly the faces `qrhs` reads for interior-water cells (`mask == 2`): a `u`
face is corrected when the cell on either side of it is interior water, and
likewise for `v`. A NaN emitted outside those faces cannot leak into the
solution.

### `nn_channels`

```python
qgsw.qgsw(..., nn_module=net, nn_channels=("ssh",))                      # default
qgsw.qgsw(..., nn_module=net, nn_channels=("ssh", "mask"))               # land present
qgsw.qgsw(..., nn_module=net, nn_channels=("ssh", "mask", "f0", "dx", "dy"))
```

| name | content |
|---|---|
| `"ssh"` | SSH with NaN replaced by 0 — the only channel that carries state |
| `"mask"` | `grd.mask / 2` → land `0`, boundary `0.5`, interior water `1` |
| `"f0"`, `"dx"`, `"dy"` | standardised over active interior-water cells, zero elsewhere |

Channels are stacked in the order given; names are validated, and if the module
exposes `in_channels` (as `nn.Conv2d` and `QGNetConvNet` do) a width mismatch is
caught before the time loop starts.

**The default is `("ssh",)` because that is what the paper specifies.**
arXiv:1911.08856 describes "a NN component that takes SSH as input", and its
stated 2545 parameters only add up with a single input channel.

Be deliberate before adding the others. They are pure functions of the grid, so
on a single fixed domain they are *constant across every sample* and act purely
as a positional encoding: that breaks the translation equivariance of the
correction, invites overfitting to the box, and destroys the property the paper
cares about — that the trained ConvNet is a portable `SSH → (du, dv)` map you
can unplug and drop into another model. `"mask"` earns its place when the domain
actually contains land; `"f0"`/`"dx"`/`"dy"` only start to matter if you train
across several regions or resolutions.

### Reference network

`QGNetConvNet` is the paper's exact architecture, 2545 parameters:

```python
from QGNet_Pytorch import QGNetConvNet

net = QGNetConvNet()           # Conv(1→16,3) BN lReLU Conv(16→16,3) BN lReLU Conv(16→2,1) × scale
assert sum(p.numel() for p in net.parameters()) == 2545
```

The trailing `scale` parameter starts at zero, so an untrained model reproduces
standard QG bit for bit and training departs from it smoothly.

`nn_cadence_steps` controls the evaluation frequency. With `dt=600`,
`nn_cadence_steps=6` gives an hourly correction update while the dynamics still
advance every 10 minutes. The held correction keeps its autograd graph, so
gradients reach every step it covers — memory therefore grows with rollout
length regardless of cadence.

A zero-output network reproduces the pure-QG trajectory exactly (`max|Δ| = 0`).

With a batched `Hi` the network receives `(B, C, ny, nx)` and must return
`(B, 2, ny, nx)`.

`nn_checkpoint=True` recomputes the network's activations in the backward
pass (`torch.utils.checkpoint`) instead of storing one set per call.

* **Memory:** on a 1-day, batch-16 rollout with the ResUNet it cuts peak
  memory from 8.7 to 2.2 GB, for ~15 % more time. That is what lets longer
  rollouts and bigger batches fit.
* **BatchNorm caveat:** each call runs twice, so BatchNorm running
  statistics are updated twice per call in `train()` mode.

`QGNetConvNet` uses BatchNorm, as the paper does. Inside the time loop it sees
a batch of one field per call, so in `train()` mode it normalises with that
field's statistics and updates its running statistics once per call, while
`eval()` uses the running statistics. The same weights give different forecasts
in the two modes; evaluate in the mode you will deploy in.

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

264 tests. `tests/test_parity.py` cross-checks every operator and whole
integrations (forward, backward, nudged, `snu=None`/`0`/`50`, `nitr` 1…25, with
and without land) against the NumPy model in the parent directory; those tests
skip if it cannot be found. `tests/test_regression.py` pins the float64 output
byte-for-byte, which is what makes it safe to optimise this solver.

The shared test case `ssh.nc` is the one at the repository root.

## Example notebooks

* `notebooks/validation_vs_numpy.ipynb` — NumPy↔PyTorch parity figures and
  error metrics on `ssh.nc`.
* `notebooks/qgnet_velocity_correction.ipynb` — the QG-Net hybrid experiment:
  trains `QGNetConvNet` to correct the geostrophic velocity and compares
  against standard QG and persistence. Runs on CPU in ~30 min. Note it uses the
  NATL60 **OSMOSIS** box shipped here, not the Gulf Stream box of the paper, so
  it demonstrates the machinery rather than replicating the published numbers.
* `notebooks/SSHforecastingQGSW_PyTorch_multistep.ipynb` — ResUNet velocity
  correction with a multi-day rollout curriculum.

## Notes

* Observation-space sampling follows the original `scipy.interpolate.griddata`
  path when `obsspace` is given. That path is CPU-only and non-differentiable.
  The 2019 original crashed when an observation fell in the first window
  (`UnboundLocalError`); the copy in this repository is patched, and the port
  does not crash.
* With `rappel`, `Hm` may be NaN on land only. NaN over water raises, where the
  original silently turns the whole field NaN.
* Grid construction follows the original mask logic exactly.
