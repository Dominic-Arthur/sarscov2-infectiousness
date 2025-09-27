from .variable_model import InfectiousnessParams, TOST, TOIT
# presymptomatic_fraction may exist in your variable_model; import if present.
try:
    from .variable_model import presymptomatic_fraction
except Exception:  # pragma: no cover
    presymptomatic_fraction = None  # type: ignore

__all__ = ["InfectiousnessParams", "TOST", "TOIT"]
if presymptomatic_fraction is not None:
    __all__.append("presymptomatic_fraction")
