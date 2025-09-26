import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOST

def test_tost_normalizes():
    p = InfectiousnessParams()
    d = TOST(params=p)
    x = np.linspace(-10, 10, 20001)
    area = np.trapz(d.pdf(x), x)
    assert 0.98 <= area <= 1.02

def test_tost_sampling_shape():
    d = TOST(params=InfectiousnessParams())
    s = d.rvs(1000)
    assert s.shape == (1000,)
