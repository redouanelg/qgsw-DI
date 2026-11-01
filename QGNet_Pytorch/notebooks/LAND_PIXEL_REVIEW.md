# Land-Pixel Support Review — `SSHforecastingQGSW_PyTorch_multistep.ipynb`

> **Status note.** This audit was written against an earlier revision of the
> package and parts of it are now out of date:
>
> * The NaN-halo expansion is no longer a Python loop with per-cell host syncs.
>   It is a vectorised `max_pool2d` dilation (`modgrid.py`), bit-identical to
>   the loop and ~11x faster, so §2 and the "fastest practical hit" remark in
>   the bottom line no longer apply.
> * Grid construction can now be hoisted out of the loop entirely by passing
>   `grd=` to `qgsw.qgsw(...)`, which the training notebook does.
> * The CNN input channels described in §4 are under review: the paper this
>   work follows ([arXiv:1911.08856](https://arxiv.org/abs/1911.08856)) feeds
>   the correction network **SSH only**, whereas the current code feeds five
>   channels. On a land-free domain the four extra channels are constant in
>   time and act as positional encoding. See the README.
> * Line numbers quoted below are stale.
>
> The land-handling logic it describes is otherwise still accurate.

Going deeper, here's a careful audit of the land-pixel path through every stage. Short version: yes, end-to-end it works for NaN-encoded land, **including the U-Net** — but a few real (mild) issues are worth knowing.

## 1. Data loader — only one rough edge

`load_osmo_sshdegrad_dataset` (notebook §3) is hard-coded to:
- variable name `degraded_sossheig`
- sentinel `99999.0` → NaN
- coords `nav_lon/nav_lat` or `lon/lat`

If your land pixels are NaN already (or use a different fill value), bypass that loader: define `ssh_daily_np`, `ssh_daily_lon`, `ssh_daily_lat` in globals before the cell runs — the conditional at the bottom of §3 takes that path and only does the central crop.

Also note: **`modgrid.grid` treats *any* NaN as land** (`modgrid.py:54, 74`). So satellite gaps over water will be lumped with land. Pre-fill water gaps before passing data in if you want them treated as ocean.

## 2. Grid + mask (`modgrid.py`)

- `mask=2` interior, `mask=1` 2-cell halo around each interior NaN, `mask=0` on the NaN cell itself (`modgrid.py:36-74`). The boundary frame (rows 0–1, last 3 rows; cols 0–1, last 3 cols) is forced to `mask=1`, then NaN cells overwrite to 0.
- The NaN-halo expansion at `modgrid.py:58-61` is a Python loop with per-cell `.cpu().item()` syncs. Correct, but **O(n_NaN) per pair** with host–device round-trips. For 365 training pairs × thousands of land cells per patch this can dominate epoch 0 (cached after that via `pair_grid_cache`). If your land mask is time-static, this is one-time; if it varies day-to-day, it pays each pair.
- Disconnected ocean basins (e.g., island in the patch) are fine in principle: each basin gets its own boundary set, and `pv2h` uses Dirichlet from the initial guess on `vp1` (`modelliptic.py:62-64`).

## 3. Solver (`qgsw` / `moddyn` / `modelliptic`)

- `h2pv` (`modelliptic.py:24-42`) explicitly zeros `q` on `mask==0` and uses a Helmholtz-style boundary on `mask==1`. Land doesn't pollute PV.
- `pv2h` (`modelliptic.py:45-91`) only solves on `valid = mask>=1`; output is filled with NaN on land via `torch.full_like(hg, torch.nan)` — so the SSH passed to the next step has NaN reinstated on land cells.
- `qrhs` zeros `rq` on `mask<=1` only when `snu>0`; with `training_snu=0.0` (notebook default), advection still uses the boundary q values from `mask==1` cells, which is the original solver's behavior. Fine.
- `h2uv` calls `nan_to_zero` on the velocity outputs (`moddyn.py:38`), so any NaN that leaked into `u/v` near coasts becomes zero before being fed to `qrhs`.

## 4. U-Net — the part you asked about

Reading `VelocityCorrectionResUNet` (notebook §5) and `_neural_input` (`qgsw.py:52-59`) together:

### What handles land correctly

- SSH → `nan_to_num(h, nan=0.0)` before entering the net.
- The `normalize_input` method (notebook §5) computes per-image SSH mean/std using `active = mask_channel > 0.75` — that's *only* the `mask==2` interior. Land/coast contribute zero to those statistics. Then `ssh_norm = ((ssh − μ) / σ) * active`, so land/coast cells are exactly zero in the normalized SSH channel. Clean.
- `f0`, `dx`, `dy` are also normalized over `active` only (`qgsw.py:40-49`). Geometry channels are zero on land/coast.
- The **mask is fed as channel 1** (`mask/2` ∈ {0, 0.5, 1}), so the network *can* learn to gate its own response by land vs water.
- Head conv weights are zero-init (notebook §5), so the initial correction is exactly zero everywhere — no garbage on land at step 0.
- The output is gated by `_neural_u_valid_mask` / `_neural_v_valid_mask` before being added to `u/v` (`qgsw.py:97-102`). Those masks are `(mask==2) & (away from patch edges)`. **Note:** the geometric "2 in / 2 in" trim is from the patch boundary, not from the coast — so a `mask==2` cell that sits right next to land will still receive its correction. The correction simply isn't applied on land/coast cells themselves.

### Real but mild issues

1. **GroupNorm sees land.** `ResUNetBlock` uses `GroupNorm` (notebook §5) which computes mean/std over (H, W) per group — *including* the land cells (where input was forced to 0). With heavy land coverage the GN statistics get pulled toward zero, biasing every batch the same way. Not catastrophic but non-ideal. A coastal-aware norm (e.g., masked group norm) would be cleaner.
2. **Pooling/upsampling blur the coast.** `avg_pool2d` over a 2×2 block that mixes land (=0) and ocean produces a "fake" half-amplitude cell at the coarser level; bilinear upsampling smears it back. So the net learns a small halo bias near coasts. Standard CNN-on-masked-data behavior; visible mainly within ~`2^depth` cells of the coast (4 cells with `depth=2`).
3. **No explicit coast feature.** The mask channel is the only coast hint. If you wanted, adding a signed-distance-to-coast channel typically helps.
4. **Bottleneck size + land.** With `depth=2`, the bottleneck is 32×32 for a 128×128 patch. If most of those bottleneck cells correspond to mostly-land receptive fields, capacity drops. Reasonable; just check `mask==2` count.

## 5. Training loop and loss

- `finite_mse` (notebook §6) intersects `isfinite(pred) & isfinite(target)` with the per-pair `mask==2` mask. NaN targets on land are excluded automatically.
- Multi-step rollout: `ssh_daily[pi+1:pi+1+k_days]` carries NaN on land in target frames — same masking applies per day.
- `loss_mask_mode="interior"` is the right setting for land. If you switch to `"finite"` you lose the 2-cell coastal halo exclusion, and gradients from coastal noise will leak in.
- PV side loss (`loss_pv_weight > 0`) computes `h2pv` on predictions and targets; `h2pv` already zeros out non-interior PV, so it's consistent with the SSH masking. Safe.

## 6. Things to actually verify in your data

After §3 of the notebook, before training, check:

```python
m = (modgrid_torch.grid(ssh_daily[0], c_train, training_snu, lon_train, lat_train).mask == 2.0)
print("Interior cells:", int(m.sum()), "/", m.numel())
```

- If interior count is very low (<~30% of the patch), expect noisy training. Either move the crop, enlarge `training_patch_size` (and adjust `resunet_depth`), or extend `train_days`.
- If your land mask varies in time, drop `pair_grid_cache` reuse or accept the rebuild cost.
- Also confirm `ssh_daily_c` shape matches; default fills 1.5 m everywhere, only used on valid cells, so land values don't matter.

## Bottom line

The framework — solver, mask, loss, U-Net wiring — is genuinely land-aware. The U-Net handles it through (a) per-image SSH normalization restricted to interior, (b) a mask channel as input, (c) zero-init head, (d) post-hoc valid-mask gating on the output. The two real weaknesses are GroupNorm including land cells in its stats and avg-pool/bilinear smearing across coastlines — both mild, both fixable later if accuracy near coasts matters. The fastest practical hit is the Python NaN-halo loop in `modgrid.py` for many land cells; cache amortizes it after the first epoch as long as the land mask is static.
