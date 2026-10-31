"""PyTorch port of the qgsw-DI quasigeostrophic shallow-water model.

The submodules mirror the original NumPy package::

    moddyn       geostrophic velocity (h2uv) and PV tendency (qrhs)
    modelliptic  SSH <-> PV transforms (h2pv, pv2h)
    modgrid      grid metrics, mask and index tables
    qgsw         forward propagator

Both the submodules and the main entry points are exported. Note that
``qgsw_fn`` is the *function* while ``qgsw`` is the *module*; importing the
module keeps ``from QGNet_Pytorch import qgsw; qgsw.qgsw(...)`` working the same
way as in the original package.
"""

from . import moddyn, modelliptic, modgrid, qgsw
from .modgrid import Grid
from .modgrid import grid as make_grid
from .qgsw import qgsw as qgsw_fn

__all__ = [
    "moddyn",
    "modelliptic",
    "modgrid",
    "qgsw",
    "Grid",
    "make_grid",
    "qgsw_fn",
]

__version__ = "0.2.0"
