"""
sarscov2_infectiousness: infectiousness models and related tools.
Keep top-level imports light to avoid hard optional dependencies.
"""

from .infectiousness import InfectiousnessParams, TOST, TOIT

# Export helper if present
try:
    from .infectiousness import presymptomatic_fraction  # type: ignore
except Exception:  # pragma: no cover
    presymptomatic_fraction = None  # type: ignore

__all__ = ["InfectiousnessParams", "TOST", "TOIT"]
if presymptomatic_fraction is not None:
    __all__.append("presymptomatic_fraction")

__version__ = "0.1.0"
