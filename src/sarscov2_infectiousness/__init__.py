"""
sarscov2_infectiousness

Variable infectiousness (E/P/I) model for SARS-CoV-2 following:
Hart WS, Maini PK, Thompson RN (2021), eLife 10:e65534.
"""

from variable_model import InfectiousnessParams, TOST, TOIT, presymptomatic_fraction

__all__ = ["InfectiousnessParams", "TOST", "TOIT", "presymptomatic_fraction"]
__version__ = "0.1.0"
