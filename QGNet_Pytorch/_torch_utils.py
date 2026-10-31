"""Shared tensor helpers for the PyTorch QGSW port."""

from __future__ import annotations

import math
import warnings
from typing import Any, Optional, Tuple

import torch

#: Floating dtypes ordered from least to most precise, used to pick a common
#: dtype when the caller mixes precisions.
_DTYPE_RANK = {torch.float16: 0, torch.bfloat16: 0, torch.float32: 1, torch.float64: 2}


def infer_dtype_device(*values: Any, dtype: Optional[torch.dtype] = None, device: Optional[Any] = None) -> Tuple[torch.dtype, torch.device]:
    """Infer a floating dtype and device from existing tensors.

    When the inputs mix floating precisions the *highest* precision wins and a
    warning is emitted, so a stray ``float32`` argument cannot silently demote
    a ``float64`` run (which costs ~7 digits of parity with the NumPy model).
    All tensors must live on the same device.
    """
    found_dtype = None
    found_device = None
    mixed = False
    for value in values:
        if not isinstance(value, torch.Tensor):
            continue
        if found_device is None:
            found_device = value.device
        elif value.device != found_device:
            raise ValueError(
                f"all tensor inputs must be on the same device, got {found_device} and {value.device}"
            )
        if value.is_floating_point():
            if found_dtype is None:
                found_dtype = value.dtype
            elif value.dtype != found_dtype:
                mixed = True
                if _DTYPE_RANK.get(value.dtype, -1) > _DTYPE_RANK.get(found_dtype, -1):
                    found_dtype = value.dtype
    if mixed and dtype is None:
        warnings.warn(
            f"mixed floating dtypes in inputs; promoting everything to {found_dtype}. "
            "Pass a single dtype explicitly to silence this.",
            RuntimeWarning,
            stacklevel=3,
        )
    if dtype is None:
        dtype = found_dtype or torch.float64
    if device is None:
        device = found_device or torch.device("cpu")
    return dtype, torch.device(device)


def as_tensor(value: Any, *, dtype: Optional[torch.dtype] = None, device: Optional[Any] = None, like: Optional[torch.Tensor] = None) -> Optional[torch.Tensor]:
    """Convert arrays/scalars to tensors without losing device or dtype intent."""
    if value is None:
        return None
    if like is not None:
        if dtype is None and like.is_floating_point():
            dtype = like.dtype
        if device is None:
            device = like.device
    if isinstance(value, torch.Tensor):
        target_dtype = dtype if dtype is not None else value.dtype
        target_device = torch.device(device) if device is not None else value.device
        return value.to(dtype=target_dtype, device=target_device)
    return torch.as_tensor(value, dtype=dtype, device=device)


def sign_scalar(value: Any) -> float:
    """Return the scalar sign used by the original NumPy model."""
    if isinstance(value, torch.Tensor):
        scalar = float(value.detach().cpu().item())
    else:
        scalar = float(value)
    if scalar > 0:
        return 1.0
    if scalar < 0:
        return -1.0
    return 0.0


def scalar_float(value: Any) -> float:
    if isinstance(value, torch.Tensor):
        return float(value.detach().cpu().item())
    return float(value)


def nan_to_zero(x: torch.Tensor) -> torch.Tensor:
    """Replace NaN with 0, leaving +/-Inf untouched.

    This mirrors the original ``u[numpy.where(numpy.isnan(u))] = 0``.
    A bare ``torch.nan_to_num(x, nan=0.0)`` would additionally saturate +/-Inf
    to +/-3.4e38, which hides a blow-up (e.g. division by ``f0 == 0`` at the
    equator) instead of letting it surface, so the infinities are passed
    through explicitly.
    """
    return torch.nan_to_num(x, nan=0.0, posinf=math.inf, neginf=-math.inf)


def maybe_tensor_scalar(value: Any, like: torch.Tensor) -> Any:
    if value is None:
        return None
    return as_tensor(value, like=like)
