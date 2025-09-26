import numpy as np
from sarscov2_infectiousness import InfectiousnessParams, TOIT, TOST


def test_shapes_and_support():
    p = InfectiousnessParams()
    tost = TOST(params=p)
    toit = TOIT(params=p)

    x = np.array([-5.0, 0.0, 5.0])
    y = np.array([0.0, 1.0, 5.0])

    assert np.all(tost.pdf(x) >= 0)
    assert np.all(toit.pdf(y) >= 0)
    assert np.allclose(np.trapezoid(tost.pdf(np.linspace(-10,10,5001)),
                                    np.linspace(-10,10,5001)), 1, rtol=1e-2)
    assert np.allclose(np.trapezoid(toit.pdf(np.linspace(0,30,5001)),
                                    np.linspace(0,30,5001)), 1, rtol=1e-2)
    assert (tost.rvs(1000).shape == (1000,))
    assert (toit.rvs(1000).shape == (1000,))

def test_generation_time_samples():
    toit = TOIT(params=InfectiousnessParams())
    gt = toit.generation_time(1000)
    assert np.all(gt >= 0)
    assert gt.shape == (1000,)
