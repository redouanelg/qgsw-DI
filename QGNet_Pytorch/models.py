"""Reference correction networks for the QG-Net hybrid model.

:class:`QGNetConvNet` reproduces the ConvNet of

    Lguensat, Le Sommer, Metref, Cosme, Fablet,
    "Learning Generalized Quasi-Geostrophic Models Using Deep Neural Numerical
    Models", NeurIPS 2019 Workshop on ML and the Physical Sciences.
    https://arxiv.org/abs/1911.08856

Section 3 of that paper describes it as

    "a 2-layer ConvNet with 16 {3x3} filters, Batch Normalization and leaky
     ReLU activations, the output layer is a linear layer with 2 channels that
     are added to U^g and V^g respectively to form W the new velocity fields.
     QG-Net has 2545 trainable parameters"

and, on the input,

    "a NN component that takes SSH as input and yields a deterministic
     'perturbation' of GV"

The parameter count is what pins the remaining details down::

    Conv2d(1, 16, 3, bias=False)    144
    BatchNorm2d(16)                  32
    Conv2d(16, 16, 3, bias=False)  2304
    BatchNorm2d(16)                  32
    Conv2d(16, 2, 1, bias=False)     32
    scale (scalar, init 0)            1
                                   ----
                                   2545

The trailing scalar is the paper's "to ensure a stable gradient flow at the
beginning of the optimization a scalar parameter initialized as zero is
multiplied to the ConvNet velocities", so the model starts exactly at standard
QG and departs from it only as training proceeds.
"""

from __future__ import annotations

import torch
from torch import nn

__all__ = ["QGNetConvNet"]


class QGNetConvNet(nn.Module):
    """The 2545-parameter velocity-correction ConvNet of arXiv:1911.08856.

    Maps ``(B, in_channels, ny, nx)`` to ``(B, 2, ny, nx)``: channel 0 is added
    to the geostrophic ``u`` and channel 1 to ``v``.

    Args:
        in_channels: input width. Keep the default of 1 (SSH only) to match the
            paper; it is exposed so the module can also be driven with the
            extra grid channels offered by ``qgsw.qgsw(..., nn_channels=...)``.
        hidden: number of filters in the two hidden layers (16 in the paper).
        negative_slope: leaky ReLU slope.
        zero_init_scale: start from exactly standard QG, as the paper does.

    The module advertises ``in_channels`` so ``qgsw.qgsw`` can check it against
    ``nn_channels`` before entering the time loop.
    """

    def __init__(self, in_channels: int = 1, hidden: int = 16,
                 negative_slope: float = 0.01, zero_init_scale: bool = True):
        super().__init__()
        self.in_channels = int(in_channels)
        self.hidden = int(hidden)

        self.body = nn.Sequential(
            nn.Conv2d(self.in_channels, self.hidden, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(self.hidden),
            nn.LeakyReLU(negative_slope=negative_slope),
            nn.Conv2d(self.hidden, self.hidden, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(self.hidden),
            nn.LeakyReLU(negative_slope=negative_slope),
        )
        # "the output layer is a linear layer with 2 channels"
        self.head = nn.Conv2d(self.hidden, 2, kernel_size=1, bias=False)
        self.scale = nn.Parameter(torch.zeros(1) if zero_init_scale else torch.ones(1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 4:
            raise ValueError(f"expected a 4D (B, C, ny, nx) tensor, got {tuple(x.shape)}")
        if x.shape[1] != self.in_channels:
            raise ValueError(
                f"expected {self.in_channels} input channel(s), got {x.shape[1]}"
            )
        return self.scale * self.head(self.body(x))

    def extra_repr(self) -> str:
        return f"in_channels={self.in_channels}, hidden={self.hidden}"
