## Learning Generalized Quasi-Geostrophic Models Using Deep Neural Numerical Models

Lguensat et al. https://arxiv.org/abs/1911.08856

`QG_PyTorch.ipynb` is the original notebook accompanying the paper: a one-layer
quasigeostrophic model written in PyTorch, adapted for GPU use and usable for
training neural networks as done in the paper. It requires double precision
(float64), is CUDA-only as written, and does not handle land pixels.

It is kept here as the supplementary material of the publication.

### Looking for a maintained PyTorch version?

See [`../QGNet_Pytorch`](../QGNet_Pytorch). It is a complete port of the NumPy
model in this repository rather than of a land-free subset:

- land/NaN masking, diffusion (`snu`), backward integration, PV nudging and
  observation sampling are all supported,
- it runs on CPU, CUDA and Apple MPS, not only CUDA,
- it agrees with the NumPy model to ~5e-16 over 5 days and ~1e-15 over 30 days,
  checked by a pytest suite,
- the neural correction hook is still there: pass `nn_module=` to
  `qgsw.qgsw(...)` to add a learned correction to the geostrophic velocity.
