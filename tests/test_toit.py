import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOIT

def test_toit_normalizes():
    p = InfectiousnessParams()
    d = TOIT(params=p)
    y = np.linspace(0, 30, 20001)
    area = np.trapz(d.pdf(y), y)
    assert 0.98 <= area <= 1.02

def test_generation_time():
    d = TOIT(params=InfectiousnessParams())
    g = d.generation_time(1000)
    assert g.shape == (1000,)
    assert np.all(g >= 0)
