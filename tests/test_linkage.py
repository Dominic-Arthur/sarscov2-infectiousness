import numpy as np
from sarscov2_infectiousness.linkage import estimate_linkage_probability

def test_linkage_shapes():
    gd = np.array([0, 1, 2], dtype=float)
    td = np.array([-5.0, 0.0, 5.0], dtype=float)
    p = estimate_linkage_probability(
        genetic_distance=gd,
        sampling_interval=td,
        num_simulations=500,  # keep test fast
        no_intermediates=3,
        intermediate_generations=(0,1),
        rng_seed=1,
    )
    assert p.shape == (3,)
    assert np.all((p >= 0) & (p <= 1))
